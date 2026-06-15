"""gwas_mcp — association analysis MCP (12 tools).

Long-running gemma/plink tools (``run_gwas_lm``, ``run_gwas_lmm``,
``run_gwas_plink``) are wrapped via :mod:`mrbigr.mcp.longjob`; first call
returns ``status='running'`` instantly, repeated call polls. See
``mrbigr.mcp.longjob`` for the full poll contract and the universal
``wait_for_job`` companion tool.
"""
from __future__ import annotations

import json
import os
import sys
import tempfile
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
from mrbigr.mcp.longjob import JobSpec, register_wait_for_job, start_or_poll  # noqa: E402


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


def _safe_trait_name(name):
    return str(name).replace("/", ".").replace(" ", "_")


def _as_bool(value):
    if isinstance(value, bool):
        return value
    if value is None:
        return False
    return str(value).strip().lower() in {"1", "true", "yes", "y", "on"}


def _gwas_out_names(output_name, pheno_cols):
    names = []
    for pheno_name in pheno_cols:
        pheno_name_safe = _safe_trait_name(pheno_name)
        names.append(pheno_name_safe if output_name is None else f"{output_name}_{pheno_name_safe}")
    return names


def _gwas_expected_assoc(output_dir, output_name, pheno_cols, suffix=".assoc.txt"):
    output_dir = Path(output_dir).resolve()
    return [str(output_dir / f"{name}{suffix}") for name in _gwas_out_names(output_name, pheno_cols)]


def _gwas_expected_logs(expected_files, log_suffix=".log.txt", assoc_suffix=".assoc.txt"):
    """Map ``foo.assoc.txt`` -> ``foo.log.txt`` (GEMMA convention)."""
    out = []
    for f in expected_files:
        if f.endswith(assoc_suffix):
            out.append(f[: -len(assoc_suffix)] + log_suffix)
        else:
            out.append(str(Path(f).with_suffix("")) + log_suffix)
    return out


def _persist_phe(phe_df, output_dir, output_name, kind):
    """Save phenotype frame to a CSV the longjob worker can re-read."""
    out_dir = Path(output_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    phe_file = out_dir / f"{output_name or kind}.{kind}.mcp_phe.csv"
    phe_df.to_csv(phe_file)
    return str(phe_file)


def _persist_cov(cov, output_dir, output_name, kind):
    """If cov is structured (list/dict), serialize to JSON next to phe so the
    worker can re-read it. Returns a value safe to put into JobSpec.kwargs."""
    if cov is None or isinstance(cov, (str, Path)):
        return str(cov) if isinstance(cov, Path) else cov
    out_dir = Path(output_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    cov_file = out_dir / f"{output_name or kind}.{kind}.mcp_cov.json"
    cov_file.write_text(json.dumps(cov, ensure_ascii=False), encoding="utf-8")
    return {"__json_file__": str(cov_file)}


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


def _attach_auto_plots(result, auto_plot):
    runner_result = result.get("result")
    if isinstance(runner_result, dict):
        for key in ("sample_summary", "sample_summary_file"):
            if key in runner_result:
                result[key] = runner_result[key]
    if not (_as_bool(auto_plot) and result.get("status") == "completed"):
        return result
    if isinstance(runner_result, dict) and runner_result.get("plots"):
        result["plots"] = runner_result["plots"]
        return result
    result["plots"] = _auto_plot(result.get("result_files") or [])
    return result


def _clump_dirs(clump_input_dir=None, result_dir=None):
    input_dir_text = clump_input_dir if clump_input_dir is not None else "./clump_input"
    input_dir = Path(input_dir_text).expanduser().resolve()
    if result_dir is None:
        if input_dir_text != "./clump_input":
            result_path = input_dir.parent / "clump_result"
        elif os.access(Path.cwd(), os.W_OK):
            result_path = Path("./clump_result").resolve()
        else:
            result_path = Path(tempfile.gettempdir()) / "mrbigr2_mcp_jobs" / "gwas_clump" / "clump_result"
    else:
        result_path = Path(result_dir).expanduser().resolve()
    return input_dir, result_path


def _clump_expected_files(clump_input_dir, result_dir):
    expected = [str(result_dir / "gwas_clump.mcp_result.json")]
    if clump_input_dir.is_dir():
        for path in sorted(p for p in clump_input_dir.glob("*") if p.is_file()):
            phe_name = path.name.split(".")[0]
            expected.append(str(result_dir / f"{phe_name}.clumped"))
    return expected


def _read_json_file(path):
    try:
        return json.loads(Path(path).read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError):
        return {}


def register(mcp) -> None:  # type: ignore[no-untyped-def]
    @mcp.tool()
    def run_gwas_lm(phe, geno_prefix, output_name=None, output_dir=None, auto_plot=False, traits=None):
        """Linear Model GWAS using GEMMA. `phe` may be a CSV path or phenotype object; `traits` selects one or more phenotype columns.

        Long-running: first call returns status='running' instantly; re-call
        with the same arguments to poll, or call wait_for_job(job_file). Do
        not bypass via shell gemma — the server already runs it detached.
        """
        output_dir = str(Path(output_dir or "output").resolve())
        phe_df = _select_phenotype_frame(phe, traits=traits)
        pheno_cols = [str(c) for c in phe_df.columns]
        expected = _gwas_expected_assoc(output_dir, output_name, pheno_cols)
        log_files = _gwas_expected_logs(expected)
        phe_csv = _persist_phe(phe_df, output_dir, output_name, "gwas_lm")
        spec = JobSpec(
            kind="gwas_lm",
            runner={
                "type": "python",
                "module": "mrbigr.mcp.runners",
                "function": "gwas_lm_runner",
                "kwargs": {
                    "phe_csv": phe_csv,
                    "geno_prefix": geno_prefix,
                    "output_name": output_name,
                    "output_dir": output_dir,
                    "auto_plot": _as_bool(auto_plot),
                },
            },
            output_dir=output_dir,
            output_name=output_name or "gwas_lm",
            expected_files=expected,
            completion_marker={
                "type": "log_contains",
                "pattern": "total computation time",
                "log_files": log_files,
            },
            eta_seconds=300 * max(1, len(pheno_cols)),
        )
        result = start_or_poll(spec)
        return _attach_auto_plots(result, auto_plot)

    @mcp.tool()
    def run_gwas_lmm(phe, geno_prefix, output_name=None, output_dir=None, auto_plot=False, cov=None, traits=None, kinship_file=None):
        """Linear Mixed Model GWAS using GEMMA. `phe` may be a CSV path or phenotype object.

        `traits` optionally selects one or more phenotype columns by exact name,
        comma-separated names, or JSON list. `100GW` resolves to
        `100grainweight` when present.

        `kinship_file` optionally points to a precomputed GEMMA `.cXX.txt`
        matrix. When provided, the LMM reuses it instead of generating a new
        kinship file inside `output_dir`.

        Long-running: first call starts GEMMA and returns status='running'.
        Re-call with the same arguments to poll, or use wait_for_job(job_file).
        Do not bypass to shell gemma — same dedup key prevents accidental
        parallel runs.
        """
        output_dir = str(Path(output_dir or "output").resolve())
        phe_df = _select_phenotype_frame(phe, traits=traits)
        pheno_cols = [str(c) for c in phe_df.columns]
        expected = _gwas_expected_assoc(output_dir, output_name, pheno_cols)
        log_files = _gwas_expected_logs(expected)
        phe_csv = _persist_phe(phe_df, output_dir, output_name, "gwas_lmm")
        cov_payload = _persist_cov(cov, output_dir, output_name, "gwas_lmm")
        if kinship_file is not None:
            kinship_path = Path(kinship_file).expanduser().resolve()
            if not kinship_path.is_file():
                raise FileNotFoundError(f"kinship file not found: {kinship_path}")
            kinship_file = str(kinship_path)
        spec = JobSpec(
            kind="gwas_lmm",
            runner={
                "type": "python",
                "module": "mrbigr.mcp.runners",
                "function": "gwas_lmm_runner",
                "kwargs": {
                    "phe_csv": phe_csv,
                    "geno_prefix": geno_prefix,
                    "output_name": output_name,
                    "output_dir": output_dir,
                    "cov": cov_payload,
                    "kinship_file": kinship_file,
                    "auto_plot": _as_bool(auto_plot),
                },
            },
            output_dir=output_dir,
            output_name=output_name or "gwas_lmm",
            expected_files=expected,
            completion_marker={
                "type": "log_contains",
                "pattern": "total computation time",
                "log_files": log_files,
            },
            eta_seconds=900 * max(1, len(pheno_cols)),
            extras={"kinship_file": kinship_file} if kinship_file else {},
        )
        result = start_or_poll(spec)
        return _attach_auto_plots(result, auto_plot)

    @mcp.tool()
    def run_gwas_plink(phe, geno_prefix, pheno_col=None, output_name="gwas", auto_plot=False):
        """GWAS using PLINK linear regression. `phe` may be a CSV path or phenotype object.

        Long-running: first call returns status='running' instantly; re-call to
        poll or use wait_for_job(job_file). Do not run plink in shell — the
        server already runs it detached and dedups by output_name.
        """
        phe_df = _phenotype_frame(phe)
        # PLINK gwas writes into the current working dir of the runner; we
        # therefore use the output_name's parent (default cwd) as output_dir.
        out_path = Path(output_name)
        if out_path.is_absolute() or out_path.parent != Path("."):
            output_dir = str(out_path.parent.resolve())
            out_stem = out_path.name
        else:
            output_dir = str(Path.cwd().resolve())
            out_stem = output_name
        expected = [str(Path(output_dir) / f"{out_stem}.assoc.linear")]
        phe_csv = _persist_phe(phe_df, output_dir, out_stem, "gwas_plink")
        spec = JobSpec(
            kind="gwas_plink",
            runner={
                "type": "python",
                "module": "mrbigr.mcp.runners",
                "function": "gwas_plink_runner",
                "kwargs": {
                    "phe_csv": phe_csv,
                    "geno_prefix": geno_prefix,
                    "pheno_col": pheno_col,
                    "output_name": str(Path(output_dir) / out_stem),
                    "auto_plot": _as_bool(auto_plot),
                },
            },
            output_dir=output_dir,
            output_name=out_stem,
            expected_files=expected,
            completion_marker={"type": "file_nonempty"},
            eta_seconds=600,
        )
        result = start_or_poll(spec)
        return _attach_auto_plots(result, auto_plot)

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
        from mrbigr.core._argjson import maybe_json_loads
        Y = maybe_json_loads(Y)
        G = maybe_json_loads(G)
        cov = maybe_json_loads(cov)
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
        """GWAS result clumping using PLINK.

        Long-running: first call returns status='running' instantly; re-call to
        poll or use wait_for_job(job_file). Completion requires the result JSON
        plus each expected PLINK .clumped file for the current clump_input_dir.
        """
        p1 = float(p1)
        p2 = float(p2)
        clump_kb = int(clump_kb)
        clump_r2 = float(clump_r2)
        clump_dir, out_dir = _clump_dirs(clump_input_dir=clump_input_dir, result_dir=result_dir)
        expected = _clump_expected_files(clump_dir, out_dir)
        spec = JobSpec(
            kind="gwas_clump",
            runner={
                "type": "python",
                "module": "mrbigr.mcp.runners",
                "function": "gwas_clump_runner",
                "kwargs": {
                    "geno_prefix": geno_prefix,
                    "p1": p1,
                    "p2": p2,
                    "clump_input_dir": str(clump_dir),
                    "result_dir": str(out_dir),
                    "clump_kb": clump_kb,
                    "clump_r2": clump_r2,
                    "result_json": str(out_dir / "gwas_clump.mcp_result.json"),
                },
            },
            output_dir=str(out_dir),
            output_name="gwas_clump",
            expected_files=expected,
            completion_marker={"type": "file_nonempty"},
            eta_seconds=600,
        )
        response = start_or_poll(spec)
        result_json = out_dir / "gwas_clump.mcp_result.json"
        if response.get("status") == "completed":
            payload = _read_json_file(result_json)
            response["result"] = payload.get("result_dir", str(out_dir))
            response["result_dir"] = payload.get("result_dir", str(out_dir))
            response["clumped_files"] = payload.get("clumped_files", [])
        return response

    @mcp.tool()
    def get_top_snps(gwas_file, n=10):
        """Get top SNPs from GWAS results."""
        return gwas.get_top_snps(gwas_file, n=int(n)).to_dict()

    @mcp.tool()
    def calculate_lambda(pvalues):
        """Calculate genomic inflation factor (lambda)."""
        import numpy as np
        import pandas as pd

        values = pvalues
        if isinstance(pvalues, (str, Path)):
            path = Path(str(pvalues)).expanduser()
            if path.is_file():
                # Auto-detect separator (tab/comma/whitespace) so .assoc.txt,
                # .csv and .tsv all work; recognize the p-value column under any
                # of the names used across the toolbox (consistent with
                # qq_plot_data / manhattan_plot_data / get_top_snps).
                df = pd.read_csv(path, sep=None, engine='python')
                candidates = ['p_wald', 'p_score', 'p_lrt', 'pvalue', 'p_value',
                              'p-value', 'p', 'p_bolt_lmm', 'p_bolt_lmm_inf']
                lower = {str(c).lower(): c for c in df.columns}
                pcol = next((lower[c] for c in candidates if c in lower), None)
                if pcol is None:
                    raise ValueError(
                        f"no p-value column found in GWAS file: {path}. "
                        f"Columns present: {list(df.columns)}. Pass a file with one of "
                        f"{candidates}, or pass an array/comma-list of p-values."
                    )
                values = df[pcol].dropna().to_numpy()
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
        from mrbigr.core._argjson import maybe_json_loads
        return gwas.qq_plot(np.array(maybe_json_loads(pvalues)))

    @mcp.tool()
    def manhattan_plot_data(gwas_results):
        """Generate Manhattan plot data from GWAS results."""
        from mrbigr.core._argjson import coerce_table
        return gwas.manhattan_plot(coerce_table(gwas_results)).to_dict()

    register_wait_for_job(mcp)


def main() -> None:
    from fastmcp import FastMCP
    mcp = FastMCP("gwas_mcp")
    register(mcp)
    mcp.run()


if __name__ == "__main__":
    main()
