"""gwas_mcp — association analysis MCP (12 tools)."""
from __future__ import annotations

import json
import os
import subprocess
import sys
import time
from pathlib import Path


def _bootstrap_syspath() -> None:
    here = Path(__file__).resolve()
    for candidate in here.parents:
        src_dir = candidate / "src" / "mrbigr"
        if src_dir.is_dir():
            src_parent = str(src_dir.parent)
            if src_parent not in sys.path:
                sys.path.insert(0, src_parent)
            return


_bootstrap_syspath()

from mrbigr.core import gwas, vis  # noqa: E402

_GWAS_LMM_PROCS = {}


def _phenotype_frame(phe):
    """Return a phenotype DataFrame from a CSV path or MCP JSON object."""
    import pandas as pd

    if isinstance(phe, pd.DataFrame):
        df = phe.copy()
    elif isinstance(phe, (str, Path)):
        path = Path(phe).expanduser()
        if not path.is_file():
            raise FileNotFoundError(f"phenotype CSV not found: {path}")
        df = pd.read_csv(path, compression="infer", sep=None, engine="python", index_col=0)
    else:
        df = pd.DataFrame(phe)
        sample_col = _sample_id_column(df)
        if sample_col is not None:
            df = df.set_index(sample_col)

    if df.empty:
        raise ValueError("phenotype data is empty")
    if len(df.columns) == 0:
        raise ValueError("phenotype data must contain at least one trait column")
    df.index = df.index.astype(str)
    return df


def _sample_id_column(df):
    id_names = {
        "id",
        "iid",
        "sample",
        "sample_id",
        "sampleid",
        "genotype",
        "accession",
    }
    for col in df.columns:
        col_text = str(col).strip()
        if col_text.lower() in id_names or col_text.startswith("Unnamed:"):
            return col
    return None


def _coerce_column_list(columns):
    if columns is None:
        return None
    if isinstance(columns, str):
        text = columns.strip()
        if not text:
            return None
        try:
            parsed = json.loads(text)
        except json.JSONDecodeError:
            parsed = None
        if isinstance(parsed, list):
            values = parsed
        elif "," in text:
            values = [part.strip() for part in text.split(",")]
        else:
            values = [text]
    elif isinstance(columns, (list, tuple, set)):
        values = list(columns)
    else:
        values = [columns]
    return [str(value).strip() for value in values if str(value).strip()]


def _normalize_column_name(value):
    return "".join(ch for ch in str(value).lower() if ch.isalnum())


def _resolve_columns(df, columns):
    requested = _coerce_column_list(columns)
    if requested is None:
        return list(df.columns)

    exact = {str(col): col for col in df.columns}
    normalized = {}
    for col in df.columns:
        normalized.setdefault(_normalize_column_name(col), col)

    aliases = {
        "100gw": ["100grainweight", "100grainwt", "hundredgrainweight"],
    }

    resolved = []
    for name in requested:
        if name in exact:
            col = exact[name]
        else:
            norm = _normalize_column_name(name)
            col = normalized.get(norm)
            if col is None:
                for alias in aliases.get(norm, []):
                    col = normalized.get(alias)
                    if col is not None:
                        break
        if col is None:
            available = ", ".join(str(col) for col in df.columns)
            raise KeyError(f"phenotype column not found: {name}; available columns: {available}")
        if col not in resolved:
            resolved.append(col)
    return resolved


def _select_phenotype_frame(phe, traits=None):
    df = _phenotype_frame(phe)
    return df.loc[:, _resolve_columns(df, traits)]


def _pid_running(pid):
    try:
        os.kill(int(pid), 0)
    except (OSError, ValueError, TypeError):
        return False
    return True


def _read_job(job_file):
    try:
        return json.loads(Path(job_file).read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError):
        return None


def _write_job(job_file, data):
    job_file = Path(job_file)
    job_file.parent.mkdir(parents=True, exist_ok=True)
    job_file.write_text(json.dumps(data, ensure_ascii=False, indent=2), encoding="utf-8")


def _safe_trait_name(name):
    return str(name).replace("/", ".").replace(" ", "_")


def _as_bool(value):
    if isinstance(value, bool):
        return value
    if value is None:
        return False
    return str(value).strip().lower() in {"1", "true", "yes", "y", "on"}


def _gwas_lmm_out_names(output_name, pheno_cols):
    names = []
    for pheno_name in pheno_cols:
        pheno_name_safe = _safe_trait_name(pheno_name)
        names.append(pheno_name_safe if output_name is None else f"{output_name}_{pheno_name_safe}")
    return names


def _gwas_lmm_expected_files(output_dir, output_name, pheno_cols):
    output_dir = Path(output_dir).resolve()
    return [output_dir / f"{out_name}.assoc.txt" for out_name in _gwas_lmm_out_names(output_name, pheno_cols)]


def _gwas_lmm_output_complete(output_dir, output_name, pheno_cols):
    expected = _gwas_lmm_expected_files(output_dir, output_name, pheno_cols)
    if not expected or not all(path.is_file() and path.stat().st_size > 0 for path in expected):
        return False
    for assoc_file in expected:
        log_file = assoc_file.with_suffix("").with_suffix(".log.txt")
        if not log_file.is_file():
            return False
        try:
            log_text = log_file.read_text(errors="ignore")
        except OSError:
            return False
        if "total computation time" not in log_text:
            return False
    return True


def _gwas_lmm_job_file(output_dir, output_name):
    output_dir = Path(output_dir or "output").resolve()
    stem = output_name or "gwas_lmm"
    return output_dir / f"{stem}.mcp_job.json"


def _gwas_lmm_status(phe, geno_prefix, output_name=None, output_dir=None, auto_plot=False, cov=None, traits=None):
    import pandas as pd

    auto_plot = _as_bool(auto_plot)
    output_dir = str(Path(output_dir or "output").resolve())
    phe_df = _select_phenotype_frame(phe, traits=traits)
    pheno_cols = [str(col) for col in phe_df.columns]
    job_file = _gwas_lmm_job_file(output_dir, output_name)
    expected_files = [str(path) for path in _gwas_lmm_expected_files(output_dir, output_name, pheno_cols)]

    if _gwas_lmm_output_complete(output_dir, output_name, pheno_cols):
        if auto_plot:
            return {"status": "completed", "gwas_results": expected_files, "plots": _auto_plot(expected_files)}
        return {"status": "completed", "gwas_results": expected_files}

    job = _read_job(job_file)
    if job:
        proc = _GWAS_LMM_PROCS.get(str(job_file))
        returncode = proc.poll() if proc is not None else None
        if returncode is None and (proc is not None or _pid_running(job.get("pid"))):
            return {
                "status": "running",
                "job_id": job.get("job_id"),
                "pid": job.get("pid"),
                "output_dir": output_dir,
                "expected_files": expected_files,
                "started_at": job.get("started_at"),
            }
        if _gwas_lmm_output_complete(output_dir, output_name, pheno_cols):
            job["status"] = "completed"
            job["completed_at"] = time.time()
            job["returncode"] = 0 if returncode is None else returncode
            _write_job(job_file, job)
            if auto_plot:
                return {"status": "completed", "gwas_results": expected_files, "plots": _auto_plot(expected_files)}
            return {"status": "completed", "gwas_results": expected_files}
        job["status"] = "failed"
        job["completed_at"] = time.time()
        job["returncode"] = returncode
        _write_job(job_file, job)
        return {
            "status": "failed",
            "result": False,
            "output_dir": output_dir,
            "expected_files": expected_files,
            "returncode": returncode,
            "stdout_file": job.get("stdout_file"),
            "stderr_file": job.get("stderr_file"),
        }

    Path(output_dir).mkdir(parents=True, exist_ok=True)
    phe_file = Path(output_dir) / f"{output_name or 'gwas_lmm'}.mcp_phe.csv"
    phe_df.to_csv(phe_file)
    cov_payload = cov
    if cov is not None and not isinstance(cov, (str, Path)):
        cov_file = Path(output_dir) / f"{output_name or 'gwas_lmm'}.mcp_cov.json"
        cov_file.write_text(json.dumps(cov, ensure_ascii=False), encoding="utf-8")
        cov_payload = {"json_file": str(cov_file)}

    job = {
        "job_id": f"gwas_lmm:{output_dir}:{output_name or 'gwas_lmm'}",
        "status": "running",
        "phe_file": str(phe_file),
        "geno_prefix": geno_prefix,
        "output_name": output_name,
        "output_dir": output_dir,
        "auto_plot": bool(auto_plot),
        "cov": str(cov) if isinstance(cov, Path) else cov_payload,
        "traits": _coerce_column_list(traits),
        "pheno_cols": pheno_cols,
        "expected_files": expected_files,
        "started_at": time.time(),
    }
    _write_job(job_file, job)
    stdout_file = Path(output_dir) / f"{output_name or 'gwas_lmm'}.mcp_stdout.log"
    stderr_file = Path(output_dir) / f"{output_name or 'gwas_lmm'}.mcp_stderr.log"
    with stdout_file.open("w", encoding="utf-8") as stdout_fh, stderr_file.open("w", encoding="utf-8") as stderr_fh:
        proc = subprocess.Popen(
            [sys.executable, str(Path(__file__).resolve()), "_run_gwas_lmm_job", str(job_file)],
            stdout=stdout_fh,
            stderr=stderr_fh,
            text=True,
            start_new_session=True,
        )
    job["pid"] = proc.pid
    job["stdout_file"] = str(stdout_file)
    job["stderr_file"] = str(stderr_file)
    _GWAS_LMM_PROCS[str(job_file)] = proc
    _write_job(job_file, job)
    return {
        "status": "running",
        "job_id": job["job_id"],
        "pid": proc.pid,
        "output_dir": output_dir,
        "expected_files": expected_files,
        "started_at": job["started_at"],
    }


def _run_gwas_lmm_job(job_file):
    import pandas as pd

    job_file = Path(job_file)
    job = _read_job(job_file)
    if not job:
        raise SystemExit(f"cannot read job file: {job_file}")
    try:
        phe_df = pd.read_csv(job["phe_file"], index_col=0)
        cov = job.get("cov")
        if isinstance(cov, dict) and cov.get("json_file"):
            cov = json.loads(Path(cov["json_file"]).read_text(encoding="utf-8"))
        result_files = gwas.gwas_lmm(
            phe_df,
            job["geno_prefix"],
            output_name=job.get("output_name"),
            output_dir=job.get("output_dir"),
            cov=cov,
        )
        job["status"] = "completed" if result_files else "failed"
        job["result_files"] = result_files
        job["completed_at"] = time.time()
        job["returncode"] = 0 if result_files else 1
    except Exception as exc:  # noqa: BLE001 - persist worker failure for polling.
        job["status"] = "failed"
        job["completed_at"] = time.time()
        job["returncode"] = 1
        job["error"] = repr(exc)
    _write_job(job_file, job)
    if job.get("returncode"):
        raise SystemExit(job.get("returncode"))


def _auto_plot(result_files):
    plots = []
    for f in (result_files or []):
        if not f or not os.path.exists(f):
            continue
        stem = f.replace('.assoc.txt', '').replace('.assoc.linear', '')
        mp = vis.manhattan_plot(f, output_file=f"{stem}_manhattan.png")
        qp = vis.qq_plot(f, output_file=f"{stem}_qq.png")
        plots.append({"file": f, "manhattan": mp, "qq": qp})
    return plots


def register(mcp) -> None:  # type: ignore[no-untyped-def]
    @mcp.tool()
    def run_gwas_lm(phe, geno_prefix, output_name=None, output_dir=None, auto_plot=False, traits=None):
        """Linear Model GWAS using GEMMA. `phe` may be a CSV path or phenotype object; `traits` selects one or more phenotype columns."""
        result_files = gwas.gwas_lm(_select_phenotype_frame(phe, traits=traits), geno_prefix, output_name=output_name, output_dir=output_dir)
        if _as_bool(auto_plot) and result_files:
            return {"gwas_results": result_files, "plots": _auto_plot(result_files)}
        return result_files

    @mcp.tool()
    def run_gwas_lmm(phe, geno_prefix, output_name=None, output_dir=None, auto_plot=False, cov=None, traits=None):
        """Linear Mixed Model GWAS using GEMMA. `phe` may be a CSV path or phenotype object.

        `traits` optionally selects one or more phenotype columns by exact name,
        comma-separated names, or JSON list. `100GW` resolves to
        `100grainweight` when present.

        Long-running mode: first call starts GEMMA and returns status='running'.
        Repeating the same call polls until status='completed' or 'failed'.
        """
        return _gwas_lmm_status(
            phe, geno_prefix,
            output_name=output_name, output_dir=output_dir,
            auto_plot=auto_plot, cov=cov, traits=traits,
        )

    @mcp.tool()
    def run_gwas_plink(phe, geno_prefix, pheno_col=None, output_name="gwas", auto_plot=False):
        """GWAS using PLINK linear regression. `phe` may be a CSV path or phenotype object."""
        result_file = gwas.gwas_plink(_phenotype_frame(phe), geno_prefix, pheno_col=pheno_col, output_name=output_name)
        if _as_bool(auto_plot) and result_file:
            return {"gwas_results": [result_file], "plots": _auto_plot([result_file])}
        return result_file

    @mcp.tool()
    def add_rs_id_to_vcf(vcf_file, output_file=None):
        """Add rs IDs to VCF file if SNP IDs are missing."""
        return gwas.add_rs_id_to_vcf(vcf_file, output_file=output_file)

    @mcp.tool()
    def ensure_rs_id(geno_prefix):
        """Ensure PLINK bim file has rs IDs for SNPs."""
        return gwas.ensure_rs_id(geno_prefix)

    @mcp.tool()
    def run_simple_gwas(Y, G, cov=None):
        """Simple linear regression GWAS (pure Python, no external tools)."""
        import numpy as np
        cov_arr = np.array(cov) if cov is not None else None
        result = gwas.simple_gwas(np.array(Y), np.array(G), cov=cov_arr)
        return result.to_dict()

    @mcp.tool()
    def generate_clump(gwas_dir):
        """Generate clump input files from GWAS results."""
        return gwas.generate_clump_input(gwas_dir)

    @mcp.tool()
    def run_gwas_clump(geno_prefix, p1=0.001, p2=0.05, clump_input_dir=None, result_dir=None,
                       clump_kb=500, clump_r2=0.1):
        """GWAS result clumping using PLINK."""
        p1 = float(p1)
        p2 = float(p2)
        clump_kb = int(clump_kb)
        clump_r2 = float(clump_r2)
        return gwas.gwas_clump(
            geno_prefix, p1=p1, p2=p2,
            clump_input_dir=clump_input_dir, result_dir=result_dir,
            clump_kb=clump_kb, clump_r2=clump_r2,
        )

    @mcp.tool()
    def get_top_snps(gwas_file, n=10):
        """Get top SNPs from GWAS results."""
        return gwas.get_top_snps(gwas_file, n=n).to_dict()

    @mcp.tool()
    def calculate_lambda(pvalues):
        """Calculate genomic inflation factor (lambda)."""
        import numpy as np
        import pandas as pd

        values = pvalues
        if isinstance(pvalues, (str, Path)):
            path = Path(str(pvalues)).expanduser()
            if path.is_file():
                df = pd.read_csv(path, sep='\t')
                if 'p_wald' in df.columns:
                    values = df['p_wald'].dropna().to_numpy()
                elif 'p_score' in df.columns:
                    values = df['p_score'].dropna().to_numpy()
                else:
                    raise ValueError(f"no p-value column found in GWAS file: {path}")
            else:
                try:
                    values = json.loads(str(pvalues))
                except json.JSONDecodeError:
                    values = [float(v) for v in str(pvalues).split(',') if str(v).strip()]
        arr = np.array(values, dtype=float)
        return {"lambda": float(gwas.calculate_lambda(arr)), "n": int(arr.size)}

    @mcp.tool()
    def qq_plot_data(pvalues):
        """Generate QQ plot data from p-values."""
        import numpy as np
        return gwas.qq_plot(np.array(pvalues))

    @mcp.tool()
    def manhattan_plot_data(gwas_results):
        """Generate Manhattan plot data from GWAS results."""
        import pandas as pd
        return gwas.manhattan_plot(pd.DataFrame(gwas_results)).to_dict()


def main() -> None:
    from fastmcp import FastMCP
    mcp = FastMCP("gwas_mcp")
    register(mcp)
    mcp.run()


if __name__ == "__main__":
    if len(sys.argv) == 3 and sys.argv[1] == "_run_gwas_lmm_job":
        _run_gwas_lmm_job(sys.argv[2])
    else:
        main()
