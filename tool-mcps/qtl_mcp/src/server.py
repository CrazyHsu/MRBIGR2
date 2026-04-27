"""qtl_mcp — QTL analysis MCP (11 tools)."""
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

from mrbigr.core import qtl  # noqa: E402

_QTL_EXTRACT_PROCS = {}


def _pid_running(pid):
    try:
        os.kill(int(pid), 0)
    except (OSError, ValueError, TypeError):
        return False
    return True


def _read_job(job_file):
    try:
        return json.loads(job_file.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError):
        return None


def _write_json(path, data):
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(data, ensure_ascii=False, indent=2), encoding="utf-8")


def _qtl_extract_job_file(output_prefix):
    return Path(f"{output_prefix}.mcp_job.json")


def _qtl_extract_result_file(output_prefix):
    return Path(f"{output_prefix}.mcp_result.json")


def _qtl_extract_expected_prefixes(qtl_data, output_prefix):
    return [f"{output_prefix}_qtl_{idx}" for idx, _ in qtl_data.iterrows()]


def _qtl_extract_output_complete(expected_prefixes):
    for prefix in expected_prefixes:
        expected = [Path(f"{prefix}{suffix}") for suffix in (".bed", ".bim", ".fam")]
        if not all(path.is_file() and path.stat().st_size > 0 for path in expected):
            return False
        log_file = Path(f"{prefix}.log")
        if not log_file.is_file():
            return False
        try:
            log_text = log_file.read_text(errors="ignore")
        except OSError:
            return False
        if "End time:" not in log_text:
            return False
    return True


def _clean_partial_qtl_extract_outputs(output_prefix, expected_prefixes):
    suffixes = (".bed", ".bim", ".fam", ".log", ".nosex")
    for prefix in expected_prefixes:
        for suffix in suffixes:
            path = Path(f"{prefix}{suffix}")
            try:
                if path.exists():
                    path.unlink()
            except OSError:
                pass
    for suffix in (".mcp_job.json", ".mcp_result.json", ".mcp_qtl_regions.json", ".mcp_stdout.log", ".mcp_stderr.log"):
        path = Path(f"{output_prefix}{suffix}")
        try:
            if path.exists():
                path.unlink()
        except OSError:
            pass


def _qtl_extract_completed(output_prefix, expected_prefixes):
    result_file = _qtl_extract_result_file(output_prefix)
    outputs = expected_prefixes
    if result_file.is_file():
        try:
            outputs = json.loads(result_file.read_text(encoding="utf-8")).get("outputs", outputs)
        except (OSError, json.JSONDecodeError):
            outputs = expected_prefixes
    return {
        "status": "completed",
        "result": outputs,
        "output_prefix": output_prefix,
        "files": outputs,
    }


def _qtl_extract_status(geno_prefix, qtl_data, output_prefix):
    output_prefix = str(Path(output_prefix).resolve())
    expected_prefixes = _qtl_extract_expected_prefixes(qtl_data, output_prefix)
    job_file = _qtl_extract_job_file(output_prefix)
    stdout_path = Path(f"{output_prefix}.mcp_stdout.log")
    stderr_path = Path(f"{output_prefix}.mcp_stderr.log")

    if _qtl_extract_output_complete(expected_prefixes):
        return _qtl_extract_completed(output_prefix, expected_prefixes)

    job = _read_job(job_file)
    if job:
        proc = _QTL_EXTRACT_PROCS.get(str(job_file))
        returncode = proc.poll() if proc is not None else None
        if returncode is None and (proc is not None or _pid_running(job.get("pid"))):
            return {
                "status": "running",
                "job_id": job.get("job_id"),
                "pid": job.get("pid"),
                "output_prefix": output_prefix,
                "expected_files": expected_prefixes,
                "started_at": job.get("started_at"),
            }
        if _qtl_extract_output_complete(expected_prefixes):
            job["status"] = "completed"
            job["completed_at"] = time.time()
            job["returncode"] = 0 if returncode is None else returncode
            _write_json(job_file, job)
            return _qtl_extract_completed(output_prefix, expected_prefixes)
        job["status"] = "failed"
        job["completed_at"] = time.time()
        job["returncode"] = returncode
        _write_json(job_file, job)
        return {
            "status": "failed",
            "result": False,
            "output_prefix": output_prefix,
            "returncode": returncode,
            "stdout_file": str(stdout_path),
            "stderr_file": str(stderr_path),
        }

    _clean_partial_qtl_extract_outputs(output_prefix, expected_prefixes)
    Path(output_prefix).parent.mkdir(parents=True, exist_ok=True)
    qtl_regions_file = Path(f"{output_prefix}.mcp_qtl_regions.json")
    result_file = _qtl_extract_result_file(output_prefix)
    qtl_data.to_json(qtl_regions_file, orient="split", force_ascii=False)
    helper = (
        "import json, sys; "
        "from pathlib import Path; "
        "import pandas as pd; "
        "from mrbigr.core import qtl; "
        "qtl_data = pd.read_json(sys.argv[2], orient='split'); "
        "outputs = qtl.extract_qtl_genotypes(sys.argv[1], qtl_data, sys.argv[3]); "
        "Path(sys.argv[4]).write_text(json.dumps({'outputs': outputs}, ensure_ascii=False), encoding='utf-8')"
    )
    env = os.environ.copy()
    sys_path = [path for path in sys.path if path]
    if env.get("PYTHONPATH"):
        sys_path.append(env["PYTHONPATH"])
    env["PYTHONPATH"] = os.pathsep.join(sys_path)
    stdout_fh = stdout_path.open("w", encoding="utf-8")
    stderr_fh = stderr_path.open("w", encoding="utf-8")
    proc = subprocess.Popen(
        [sys.executable, "-c", helper, geno_prefix, str(qtl_regions_file), output_prefix, str(result_file)],
        stdout=stdout_fh,
        stderr=stderr_fh,
        text=True,
        start_new_session=True,
        env=env,
    )
    stdout_fh.close()
    stderr_fh.close()
    job = {
        "job_id": f"extract_qtl_genotypes:{output_prefix}",
        "status": "running",
        "pid": proc.pid,
        "input_prefix": geno_prefix,
        "output_prefix": output_prefix,
        "expected_files": expected_prefixes,
        "started_at": time.time(),
    }
    _QTL_EXTRACT_PROCS[str(job_file)] = proc
    _write_json(job_file, job)
    return {
        "status": "running",
        "job_id": job["job_id"],
        "pid": proc.pid,
        "output_prefix": output_prefix,
        "expected_files": expected_prefixes,
        "started_at": job["started_at"],
    }


def register(mcp) -> None:  # type: ignore[no-untyped-def]
    @mcp.tool()
    def detect_qtl_regions(gwas_file, p1=1e-7, p2=1e-5, p2n=5, window=500000):
        """Detect QTL regions from GWAS results."""
        p1 = float(p1)
        p2 = float(p2)
        p2n = int(p2n)
        window = int(window)
        result = qtl.detect_qtl(gwas_file, p1=p1, p2=p2, p2n=p2n, window=window)
        if result is None:
            return None
        output_file = None
        if isinstance(gwas_file, (str, Path)):
            gwas_path = Path(gwas_file).expanduser()
            if gwas_path.is_file():
                if gwas_path.suffix == ".txt" and gwas_path.name.endswith(".assoc.txt"):
                    output_file = str(gwas_path.with_suffix("").with_suffix(".qtl.csv"))
                else:
                    output_file = str(gwas_path.with_suffix(".qtl.csv"))
                Path(output_file).parent.mkdir(parents=True, exist_ok=True)
                result.to_csv(output_file, index=False)
        return {
            "qtl_regions": result.to_dict(),
            "output_file": output_file,
            "n_regions": int(len(result)),
        }

    @mcp.tool()
    def get_lead_snp(gwas_file, region_chr, region_start, region_end):
        """Get lead SNP in a region."""
        region_start = int(region_start)
        region_end = int(region_end)
        result = qtl.get_lead_snp(gwas_file, region_chr, region_start, region_end)
        return result.to_dict() if result is not None else None

    @mcp.tool()
    def identify_peak_snps(gwas_dir, p_threshold=1e-5, max_peaks=50):
        """Identify peak SNPs from GWAS results."""
        p_threshold = float(p_threshold)
        max_peaks = int(max_peaks)
        result = qtl.identify_peak_snps(gwas_dir, p_threshold=p_threshold, max_peaks=max_peaks)
        return result.to_dict() if result is not None else None

    @mcp.tool()
    def extract_qtl_genotypes(geno_prefix, qtl_regions, output_prefix):
        """Extract genotypes for QTL regions.

        Long-running mode: first call starts a background PLINK extraction and
        returns status='running'. Repeating the same call polls until completed.
        """
        import pandas as pd
        if isinstance(qtl_regions, dict):
            qtl_data = pd.DataFrame(qtl_regions)
        elif isinstance(qtl_regions, (str, Path)):
            qtl_data = pd.read_csv(qtl_regions)
        elif isinstance(qtl_regions, pd.DataFrame):
            qtl_data = qtl_regions
        else:
            qtl_data = pd.DataFrame(qtl_regions)
        return _qtl_extract_status(geno_prefix, qtl_data, output_prefix)

    @mcp.tool()
    def calculate_qtl_haplo(geno_prefix: str, qtl_snp_list: list[str] | str, output_prefix: str):
        """Calculate haplotypes for SNPs in a QTL region."""
        result = qtl.calculate_qtl_haplo(geno_prefix, qtl_snp_list, output_prefix)
        return result.to_dict() if result is not None else None

    @mcp.tool()
    def plot_qtl_region(gwas_file, chr_val, start, end, output_file=None,
                        highlight_snps=None, significance=5e-8, suggest=1e-5):
        """Plot GWAS results for QTL region."""
        start = int(start)
        end = int(end)
        significance = float(significance)
        suggest = float(suggest)
        return qtl.plot_qtl_region(
            gwas_file, chr_val, start, end,
            output_file=output_file,
            highlight_snps=highlight_snps,
            significance=significance,
            suggest=suggest,
        )

    @mcp.tool()
    def qtl_summary(gwas_dir, output_prefix='qtl_summary'):
        """Generate summary of all QTL regions."""
        result = qtl.qtl_summary(gwas_dir, output_prefix=output_prefix)
        return result.to_dict() if result is not None else None

    @mcp.tool()
    def map_qtl_to_genes(qtl_df, annotation_file):
        """Map QTL regions to genes."""
        import pandas as pd
        qtl_data = pd.DataFrame(qtl_df) if isinstance(qtl_df, dict) else qtl_df
        result = qtl.map_qtl_to_genes(qtl_data, annotation_file)
        return result.to_dict() if result is not None else None

    @mcp.tool()
    def qtl_enrichment_test(qtl_genes, background_genes, gene_sets):
        """Test enrichment of QTL genes against gene sets."""
        result = qtl.qtl_enrichment_test(qtl_genes, background_genes, gene_sets)
        return result.to_dict() if result is not None else None

    @mcp.tool()
    def get_qtl_overlap(qtl1, qtl2):
        """Find overlapping QTL intervals between two QTL sets."""
        import pandas as pd
        qtl1_data = pd.DataFrame(qtl1) if isinstance(qtl1, dict) else qtl1
        qtl2_data = pd.DataFrame(qtl2) if isinstance(qtl2, dict) else qtl2
        result = qtl.get_qtl_overlap(qtl1_data, qtl2_data)
        return result.to_dict() if result is not None else None

    @mcp.tool()
    def export_qtl_bed(qtl_df, output_file):
        """Export QTL as BED file."""
        import pandas as pd
        qtl_data = pd.DataFrame(qtl_df) if isinstance(qtl_df, dict) else qtl_df
        return qtl.export_qtl_bed(qtl_data, output_file)


def main() -> None:
    from fastmcp import FastMCP
    mcp = FastMCP("qtl_mcp")
    register(mcp)
    mcp.run()


if __name__ == "__main__":
    main()
