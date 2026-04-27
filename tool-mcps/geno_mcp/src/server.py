"""geno_mcp — genotype processing MCP (12 tools).

Two entry points:
  * ``register(mcp)`` — attach the 12 tools to an existing FastMCP
    instance. Used by src/server.py (the legacy aggregator).
  * ``python tool-mcps/geno_mcp/src/server.py`` — standalone stdio server.
"""
from __future__ import annotations

import os
import json
import subprocess
import sys
import time
from pathlib import Path


def _bootstrap_syspath() -> None:
    """Let this file be run directly without `pip install -e .`.

    When executed as ``python tool-mcps/geno_mcp/src/server.py`` there's
    no guarantee the mrbigr package is importable, so we walk up from
    this file and put <repo>/src on sys.path.
    """
    here = Path(__file__).resolve()
    for candidate in here.parents:
        src_dir = candidate / "src" / "mrbigr"
        if src_dir.is_dir():
            src_parent = str(src_dir.parent)
            if src_parent not in sys.path:
                sys.path.insert(0, src_parent)
            return


_bootstrap_syspath()

from mrbigr.core import geno  # noqa: E402

_SNP_QC_PROCS = {}
_KINSHIP_PROCS = {}


def _snp_qc_job_file(output_prefix):
    return Path(f"{output_prefix}.mcp_job.json")


def _snp_qc_expected_files(output_prefix):
    return [Path(f"{output_prefix}{suffix}") for suffix in (".bed", ".bim", ".fam")]


def _snp_qc_output_complete(output_prefix):
    expected = _snp_qc_expected_files(output_prefix)
    if not all(path.is_file() and path.stat().st_size > 0 for path in expected):
        return False
    log_file = Path(f"{output_prefix}.log")
    if not log_file.is_file():
        return False
    try:
        log_text = log_file.read_text(errors="ignore")
    except OSError:
        return False
    return "End time:" in log_text and "--make-bed" in log_text


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


def _write_job(job_file, data):
    job_file.parent.mkdir(parents=True, exist_ok=True)
    job_file.write_text(json.dumps(data, ensure_ascii=False, indent=2), encoding="utf-8")


def _kinship_job_file(output_prefix):
    return Path(f"{output_prefix}.mcp_job.json")


def _kinship_expected_file(output_prefix):
    return Path(f"{output_prefix}.cXX.txt")


def _kinship_output_complete(output_prefix):
    kinship_file = _kinship_expected_file(output_prefix)
    if not (kinship_file.is_file() and kinship_file.stat().st_size > 0):
        return False
    log_file = Path(f"{output_prefix}.log.txt")
    if not log_file.is_file():
        return False
    try:
        log_text = log_file.read_text(errors="ignore")
    except OSError:
        return False
    return "GEMMA Version" in log_text and "Computation Time" in log_text


def _clean_partial_kinship_outputs(output_prefix):
    for suffix in (".cXX.txt", ".log.txt", ".mcp_stdout.log", ".mcp_stderr.log"):
        path = Path(f"{output_prefix}{suffix}")
        try:
            if path.exists():
                path.unlink()
        except OSError:
            pass


def _kinship_status(input_prefix, output_prefix):
    output_prefix = str(Path(output_prefix).resolve())
    job_file = _kinship_job_file(output_prefix)
    kinship_file = _kinship_expected_file(output_prefix)
    log_file = Path(f"{output_prefix}.log.txt")
    stdout_path = Path(f"{output_prefix}.mcp_stdout.log")
    stderr_path = Path(f"{output_prefix}.mcp_stderr.log")

    if _kinship_output_complete(output_prefix):
        return {
            "status": "completed",
            "result": True,
            "kinship_file": str(kinship_file),
            "output_prefix": output_prefix,
            "log_file": str(log_file),
        }

    job = _read_job(job_file)
    if job:
        proc = _KINSHIP_PROCS.get(str(job_file))
        returncode = proc.poll() if proc is not None else None
        if returncode is None and (proc is not None or _pid_running(job.get("pid"))):
            return {
                "status": "running",
                "job_id": job.get("job_id"),
                "pid": job.get("pid"),
                "output_prefix": output_prefix,
                "started_at": job.get("started_at"),
            }
        if _kinship_output_complete(output_prefix):
            job["status"] = "completed"
            job["completed_at"] = time.time()
            job["returncode"] = 0 if returncode is None else returncode
            _write_job(job_file, job)
            return {
                "status": "completed",
                "result": True,
                "kinship_file": str(kinship_file),
                "output_prefix": output_prefix,
                "log_file": str(log_file),
            }
        job["status"] = "failed"
        job["completed_at"] = time.time()
        job["returncode"] = returncode
        _write_job(job_file, job)
        return {
            "status": "failed",
            "result": False,
            "kinship_file": None,
            "output_prefix": output_prefix,
            "returncode": returncode,
            "log_file": str(log_file),
            "stdout_file": str(stdout_path),
            "stderr_file": str(stderr_path),
        }

    _clean_partial_kinship_outputs(output_prefix)
    output_path = Path(output_prefix)
    output_dir = output_path.parent
    output_stem = output_path.name
    output_dir.mkdir(parents=True, exist_ok=True)
    cmd = [
        str(geno.SCRIPT_DIR / "utils" / "gemma.linux"),
        "-bfile", input_prefix,
        "-gk", "1",
        "-outdir", str(output_dir),
        "-o", output_stem,
    ]
    stdout_fh = stdout_path.open("w", encoding="utf-8")
    stderr_fh = stderr_path.open("w", encoding="utf-8")
    proc = subprocess.Popen(
        cmd,
        stdout=stdout_fh,
        stderr=stderr_fh,
        text=True,
        start_new_session=True,
    )
    stdout_fh.close()
    stderr_fh.close()
    job = {
        "job_id": f"kinship:{output_prefix}",
        "status": "running",
        "pid": proc.pid,
        "cmd": cmd,
        "input_prefix": input_prefix,
        "output_prefix": output_prefix,
        "started_at": time.time(),
    }
    _KINSHIP_PROCS[str(job_file)] = proc
    _write_job(job_file, job)
    return {
        "status": "running",
        "job_id": job["job_id"],
        "pid": proc.pid,
        "output_prefix": output_prefix,
        "started_at": job["started_at"],
    }


def _clean_partial_snp_qc_outputs(output_prefix):
    for suffix in (".bed", ".bim", ".fam", ".log", ".nosex", ".mcp_stdout.log", ".mcp_stderr.log"):
        path = Path(f"{output_prefix}{suffix}")
        try:
            if path.exists():
                path.unlink()
        except OSError:
            pass


def _snp_qc_status(input_prefix, output_prefix, maf, missing_rate, mind):
    job_file = _snp_qc_job_file(output_prefix)
    if _snp_qc_output_complete(output_prefix):
        return {
            "status": "completed",
            "result": True,
            "output_prefix": output_prefix,
            "files": [str(path) for path in _snp_qc_expected_files(output_prefix)],
        }

    job = _read_job(job_file)
    if job:
        proc = _SNP_QC_PROCS.get(str(job_file))
        returncode = proc.poll() if proc is not None else None
        if returncode is None and (proc is not None or _pid_running(job.get("pid"))):
            return {
                "status": "running",
                "job_id": job.get("job_id"),
                "pid": job.get("pid"),
                "output_prefix": output_prefix,
                "started_at": job.get("started_at"),
            }
        if _snp_qc_output_complete(output_prefix):
            job["status"] = "completed"
            job["completed_at"] = time.time()
            job["returncode"] = 0 if returncode is None else returncode
            _write_job(job_file, job)
            return {
                "status": "completed",
                "result": True,
                "output_prefix": output_prefix,
                "files": [str(path) for path in _snp_qc_expected_files(output_prefix)],
            }
        job["status"] = "failed"
        job["completed_at"] = time.time()
        job["returncode"] = returncode
        _write_job(job_file, job)
        return {
            "status": "failed",
            "result": False,
            "output_prefix": output_prefix,
            "returncode": returncode,
            "log_file": f"{output_prefix}.log",
            "stderr_file": f"{output_prefix}.mcp_stderr.log",
        }

    _clean_partial_snp_qc_outputs(output_prefix)
    output_path = Path(output_prefix)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    cmd = [
        geno.PLINK_BIN,
        "--bfile", input_prefix,
        "--out", output_prefix,
        "--maf", str(maf),
        "--geno", str(missing_rate),
        "--mind", str(mind),
        "--make-bed",
    ]
    stdout_path = Path(f"{output_prefix}.mcp_stdout.log")
    stderr_path = Path(f"{output_prefix}.mcp_stderr.log")
    stdout_fh = stdout_path.open("w", encoding="utf-8")
    stderr_fh = stderr_path.open("w", encoding="utf-8")
    proc = subprocess.Popen(
        cmd,
        stdout=stdout_fh,
        stderr=stderr_fh,
        text=True,
        start_new_session=True,
    )
    stdout_fh.close()
    stderr_fh.close()
    job = {
        "job_id": f"snp_qc:{output_prefix}",
        "status": "running",
        "pid": proc.pid,
        "cmd": cmd,
        "input_prefix": input_prefix,
        "output_prefix": output_prefix,
        "maf": maf,
        "missing_rate": missing_rate,
        "mind": mind,
        "started_at": time.time(),
    }
    _SNP_QC_PROCS[str(job_file)] = proc
    _write_job(job_file, job)
    return {
        "status": "running",
        "job_id": job["job_id"],
        "pid": proc.pid,
        "output_prefix": output_prefix,
        "started_at": job["started_at"],
    }


def register(mcp) -> None:  # type: ignore[no-untyped-def]
    """Attach all geno_mcp tools to the given FastMCP instance."""

    @mcp.tool()
    def run_snp_qc(input_prefix, output_prefix, maf=0.05, missing_rate=0.2, mind=0.2):
        """SNP quality control using PLINK.

        Long-running mode: first call starts the PLINK QC job and returns
        status='running'. Repeating the same call polls the job until it returns
        status='completed' or status='failed'.
        """
        return _snp_qc_status(input_prefix, output_prefix, maf, missing_rate, mind)

    @mcp.tool()
    def subset_genotype(input_prefix, output_prefix, chromosomes=None, proportion=None, seed=42):
        """Subset PLINK genotype data by chromosome or random proportion."""
        return geno.subset_plink(input_prefix, output_prefix, chromosomes=chromosomes, proportion=proportion, seed=seed)

    @mcp.tool()
    def run_genotype_pca(input_prefix, n_components=10, max_snps=20000):
        """PCA analysis for genotype data."""
        n_components = int(n_components)
        max_snps = int(max_snps) if max_snps is not None else None
        pc_df, var_ratio = geno.calculate_pca(input_prefix, n_components=n_components, max_snps=max_snps)
        input_path = Path(os.path.abspath(input_prefix))
        output_dir = input_path.parent
        if not os.access(output_dir, os.W_OK):
            output_dir = Path("/tmp/mrbigr_pca")
            output_dir.mkdir(parents=True, exist_ok=True)
        output_file = str(output_dir / f"{input_path.name}.pca.csv")
        output_df = pc_df.copy()
        if "sample" in output_df.columns:
            output_df["sample"] = output_df["sample"].astype(str)
        else:
            output_df.index = output_df.index.astype(str)
            output_df.index.name = "sample"
            output_df = output_df.reset_index()
        output_df.to_csv(output_file, index=False)
        return {
            "pc_data": pc_df.to_dict(),
            "variance_explained": var_ratio.tolist(),
            "output_file": output_file,
            "rows": int(pc_df.shape[0]),
            "columns": ["sample"] + [str(col) for col in pc_df.columns],
        }

    @mcp.tool()
    def run_calculate_ibd(input_prefix, output_prefix):
        """Calculate Identity by Descent (IBD) matrix."""
        return geno.calculate_ibd(input_prefix, output_prefix)

    @mcp.tool()
    def convert_vcf(vcf_file, output_prefix):
        """Convert VCF to PLINK format."""
        return geno.vcf_to_plink(vcf_file, output_prefix)

    @mcp.tool()
    def convert_hapmap(hapmap_file, output_prefix):
        """Convert HapMap format to PLINK."""
        return geno.hapmap_to_plink(hapmap_file, output_prefix)

    @mcp.tool()
    def run_plink_to_vcf(bed_prefix, output_prefix):
        """Convert PLINK to VCF format."""
        return geno.plink_to_vcf(bed_prefix, output_prefix)

    @mcp.tool()
    def calculate_kinship(input_prefix, output_prefix):
        """Calculate kinship/relatedness matrix using GEMMA.

        Long-running mode: first call starts GEMMA and returns status='running'.
        Repeating the same call polls until status='completed' or 'failed'.
        """
        return _kinship_status(input_prefix, output_prefix)

    @mcp.tool()
    def impute_genotype(input_prefix, output_prefix, method='mean'):
        """Impute missing genotype values."""
        return geno.snp_impute(input_prefix, output_prefix, method=method)

    @mcp.tool()
    def run_snp_pruning(input_prefix, output_prefix, window=50, shift=5, r2=0.5, maf=0.05):
        """LD-based SNP pruning using PLINK."""
        return geno.snp_pruning(input_prefix, output_prefix, window=window, shift=shift, r2=r2, maf=maf)

    @mcp.tool()
    def run_snp_clumping(input_prefix, output_prefix, r2=0.5, maf=0.05, window_kb=250):
        """LD-based SNP clumping."""
        return geno.snp_clumping(input_prefix, output_prefix, r2=r2, maf=maf, window_kb=window_kb)

    @mcp.tool()
    def get_snp_statistics(input_prefix):
        """Get basic SNP statistics."""
        return geno.get_snp_stats(input_prefix)


def main() -> None:
    from fastmcp import FastMCP
    mcp = FastMCP("geno_mcp")
    register(mcp)
    mcp.run()


if __name__ == "__main__":
    main()
