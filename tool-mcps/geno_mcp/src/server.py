"""geno_mcp — genotype processing MCP (12 tools).

All long-running tools (gemma kinship, plink QC/IBD/pruning/clumping,
VCF/HapMap conversions, mean-imputation) are wrapped via
:mod:`mrbigr.mcp.longjob`: first call returns ``status='running'`` instantly,
re-call with the same args to poll, or use the universal ``wait_for_job``
companion tool.
"""
from __future__ import annotations

import os
import sys
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
from mrbigr.mcp.longjob import JobSpec, register_wait_for_job, start_or_poll  # noqa: E402


def _split_prefix(output_prefix: str) -> tuple[str, str, str]:
    """Return ``(absolute_prefix, output_dir, output_stem)``."""
    abs_prefix = str(Path(output_prefix).resolve())
    out_path = Path(abs_prefix)
    return abs_prefix, str(out_path.parent), out_path.name


def _core_call_spec(
    *,
    kind: str,
    function: str,
    kwargs: dict,
    output_prefix: str,
    expected_suffixes: list[str],
    log_pattern: str | None = None,
    log_suffix: str = ".log",
    eta_seconds: int = 300,
) -> JobSpec:
    """Build a JobSpec that calls ``mrbigr.core.geno.<function>(**kwargs)`` via
    the longjob worker.

    Caller passes the un-resolved ``output_prefix`` for naming, and the
    suffixes that should appear when the job completes successfully (e.g.
    ``[".bed", ".bim", ".fam"]``).
    """
    abs_prefix, output_dir, out_stem = _split_prefix(output_prefix)
    expected = [f"{abs_prefix}{suf}" for suf in expected_suffixes]
    if log_pattern:
        marker: dict = {
            "type": "log_contains",
            "pattern": log_pattern,
            "log_files": [f"{abs_prefix}{log_suffix}"],
        }
    else:
        marker = {"type": "file_nonempty"}
    return JobSpec(
        kind=kind,
        runner={
            "type": "python",
            "module": "mrbigr.mcp.runners",
            "function": "core_call_runner",
            "kwargs": {
                "module": "mrbigr.core.geno",
                "function": function,
                "kwargs": dict(kwargs),
            },
        },
        output_dir=output_dir,
        output_name=out_stem,
        expected_files=expected,
        completion_marker=marker,
        eta_seconds=eta_seconds,
    )


def register(mcp) -> None:  # type: ignore[no-untyped-def]
    """Attach all geno_mcp tools to the given FastMCP instance."""

    @mcp.tool()
    def run_snp_qc(input_prefix, output_prefix, maf=0.05, missing_rate=0.2, mind=0.2):
        """SNP quality control using PLINK.

        Long-running: first call returns status='running' instantly; re-call
        with the same arguments to poll, or use wait_for_job(job_file).
        """
        return start_or_poll(_core_call_spec(
            kind="snp_qc",
            function="snp_qc",
            kwargs={
                "input_prefix": input_prefix,
                "output_prefix": output_prefix,
                "maf": float(maf),
                "missing_rate": float(missing_rate),
                "mind": float(mind),
            },
            output_prefix=output_prefix,
            expected_suffixes=[".bed", ".bim", ".fam"],
            log_pattern="--make-bed",
            log_suffix=".log",
            eta_seconds=300,
        ))

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
            "variance_explained": var_ratio.tolist(),
            "output_file": output_file,
            "rows": int(pc_df.shape[0]),
            "columns": [str(col) for col in output_df.columns],
        }

    @mcp.tool()
    def run_calculate_ibd(input_prefix, output_prefix):
        """Calculate Identity by Descent (IBD) matrix using PLINK.

        Long-running: first call returns status='running'; re-call to poll or
        use wait_for_job(job_file).
        """
        return start_or_poll(_core_call_spec(
            kind="ibd",
            function="calculate_ibd",
            kwargs={"input_prefix": input_prefix, "output_prefix": output_prefix},
            output_prefix=output_prefix,
            expected_suffixes=[".genome"],
            log_pattern=None,
            eta_seconds=600,
        ))

    @mcp.tool()
    def convert_vcf(vcf_file, output_prefix):
        """Convert VCF to PLINK format using PLINK.

        Long-running: first call returns status='running'; re-call to poll or
        use wait_for_job(job_file).
        """
        return start_or_poll(_core_call_spec(
            kind="vcf_to_plink",
            function="vcf_to_plink",
            kwargs={"vcf_file": vcf_file, "output_prefix": output_prefix},
            output_prefix=output_prefix,
            expected_suffixes=[".bed", ".bim", ".fam"],
            log_pattern=None,
            eta_seconds=600,
        ))

    @mcp.tool()
    def convert_hapmap(hapmap_file, output_prefix):
        """Convert HapMap format to PLINK.

        Long-running: first call returns status='running'; re-call to poll or
        use wait_for_job(job_file).
        """
        return start_or_poll(_core_call_spec(
            kind="hapmap_to_plink",
            function="hapmap_to_plink",
            kwargs={"hapmap_file": hapmap_file, "output_prefix": output_prefix},
            output_prefix=output_prefix,
            expected_suffixes=[".bed", ".bim", ".fam"],
            log_pattern=None,
            eta_seconds=600,
        ))

    @mcp.tool()
    def run_plink_to_vcf(bed_prefix, output_prefix):
        """Convert PLINK to VCF format using PLINK.

        Long-running: first call returns status='running'; re-call to poll or
        use wait_for_job(job_file).
        """
        return start_or_poll(_core_call_spec(
            kind="plink_to_vcf",
            function="plink_to_vcf",
            kwargs={"bed_prefix": bed_prefix, "output_prefix": output_prefix},
            output_prefix=output_prefix,
            expected_suffixes=[".vcf"],
            log_pattern=None,
            eta_seconds=600,
        ))

    @mcp.tool()
    def calculate_kinship(input_prefix, output_prefix):
        """Calculate kinship/relatedness matrix using GEMMA.

        Long-running: first call returns status='running'; re-call with the
        same arguments to poll, or use wait_for_job(job_file). Do not bypass
        to shell gemma — same dedup key prevents accidental parallel runs.
        """
        return start_or_poll(_core_call_spec(
            kind="kinship",
            function="calculate_kinship",
            kwargs={"input_prefix": input_prefix, "output_prefix": output_prefix},
            output_prefix=output_prefix,
            expected_suffixes=[".cXX.txt"],
            log_pattern="Computation Time",
            log_suffix=".log.txt",
            eta_seconds=600,
        ))

    @mcp.tool()
    def impute_genotype(input_prefix, output_prefix, method='mean'):
        """Impute missing genotype values.

        Long-running: first call returns status='running'; re-call to poll or
        use wait_for_job(job_file).
        """
        return start_or_poll(_core_call_spec(
            kind="impute",
            function="snp_impute",
            kwargs={"input_prefix": input_prefix, "output_prefix": output_prefix, "method": method},
            output_prefix=output_prefix,
            expected_suffixes=[".bed", ".bim", ".fam"],
            log_pattern=None,
            eta_seconds=900,
        ))

    @mcp.tool()
    def run_snp_pruning(input_prefix, output_prefix, window=50, shift=5, r2=0.5, maf=0.05):
        """LD-based SNP pruning using PLINK.

        Produces ``{output_prefix}_pruned.{bed,bim,fam}``. Long-running: first
        call returns status='running'; re-call to poll or use wait_for_job.
        """
        # snp_pruning writes its results under "{output_prefix}_pruned*"; use
        # that as the resolved prefix for completion detection so the job_file
        # also lives next to it.
        pruned_prefix = f"{output_prefix}_pruned"
        return start_or_poll(_core_call_spec(
            kind="snp_pruning",
            function="snp_pruning",
            kwargs={
                "input_prefix": input_prefix,
                "output_prefix": output_prefix,
                "window": int(window),
                "shift": int(shift),
                "r2": float(r2),
                "maf": float(maf),
            },
            output_prefix=pruned_prefix,
            expected_suffixes=[".bed", ".bim", ".fam"],
            log_pattern=None,
            eta_seconds=600,
        ))

    @mcp.tool()
    def run_snp_clumping(input_prefix, output_prefix, r2=0.5, maf=0.05, window_kb=250):
        """LD-based SNP clumping (pure Python).

        Long-running: first call returns status='running'; re-call to poll or
        use wait_for_job(job_file).
        """
        return start_or_poll(_core_call_spec(
            kind="snp_clumping",
            function="snp_clumping",
            kwargs={
                "input_prefix": input_prefix,
                "output_prefix": output_prefix,
                "r2": float(r2),
                "maf": float(maf),
                "window_kb": int(window_kb),
            },
            output_prefix=output_prefix,
            expected_suffixes=[".bed", ".bim", ".fam"],
            log_pattern=None,
            eta_seconds=600,
        ))

    @mcp.tool()
    def get_snp_statistics(input_prefix):
        """Get basic SNP statistics."""
        return geno.get_snp_stats(input_prefix)

    register_wait_for_job(mcp)


def main() -> None:
    from fastmcp import FastMCP
    mcp = FastMCP("geno_mcp")
    register(mcp)
    mcp.run()


if __name__ == "__main__":
    main()
