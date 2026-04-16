"""geno_mcp — genotype processing MCP (12 tools).

Two entry points:
  * ``register(mcp)`` — attach the 12 tools to an existing FastMCP
    instance. Used by src/server.py (the legacy aggregator).
  * ``python tool-mcps/geno_mcp/src/server.py`` — standalone stdio server.
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


def register(mcp) -> None:  # type: ignore[no-untyped-def]
    """Attach all geno_mcp tools to the given FastMCP instance."""

    @mcp.tool()
    def run_snp_qc(input_prefix, output_prefix, maf=0.05, missing_rate=0.2, mind=0.2):
        """SNP quality control using PLINK."""
        return geno.snp_qc(input_prefix, output_prefix, maf=maf, missing_rate=missing_rate, mind=mind)

    @mcp.tool()
    def subset_genotype(input_prefix, output_prefix, chromosomes=None, proportion=None, seed=42):
        """Subset PLINK genotype data by chromosome or random proportion."""
        return geno.subset_plink(input_prefix, output_prefix, chromosomes=chromosomes, proportion=proportion, seed=seed)

    @mcp.tool()
    def run_genotype_pca(input_prefix, n_components=10, max_snps=20000):
        """PCA analysis for genotype data."""
        pc_df, var_ratio = geno.calculate_pca(input_prefix, n_components=n_components, max_snps=max_snps)
        return {"pc_data": pc_df.to_dict(), "variance_explained": var_ratio.tolist()}

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
        """Calculate kinship/relatedness matrix using GEMMA."""
        return geno.calculate_kinship(input_prefix, output_prefix)

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
