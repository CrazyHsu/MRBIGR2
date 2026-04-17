"""peak_mcp — peak & haplotype visualization MCP (4 tools)."""
from __future__ import annotations

import sys
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

from mrbigr.core import peak  # noqa: E402


def register(mcp) -> None:  # type: ignore[no-untyped-def]
    @mcp.tool()
    def plot_qtl_boxplot(pheno_file, geno_prefix, qtl_df):
        """Generate boxplots for QTL regions."""
        return peak.plot_qtl_boxplot(pheno_file, geno_prefix, qtl_df)

    @mcp.tool()
    def plot_grouped_boxplot(data_dict, output_file=None, test_method='t-test'):
        """Generate grouped boxplot."""
        return peak.plot_grouped_boxplot(data_dict, output_file=output_file, test_method=test_method)

    @mcp.tool()
    def haplotype_test(pheno_df, geno_df, snp_id, test_method='t-test'):
        """Perform statistical test for haplotype effect."""
        result = peak.haplotype_test(pheno_df, geno_df, snp_id, test_method=test_method)
        return result.to_dict() if result is not None else None

    @mcp.tool()
    def multi_trait_qtl_plot(gwas_dir, qtl_file, output_prefix):
        """Generate multi-trait Manhattan plot with QTL regions."""
        return peak.multi_trait_qtl_plot(gwas_dir, qtl_file, output_prefix)


def main() -> None:
    from fastmcp import FastMCP
    mcp = FastMCP("peak_mcp")
    register(mcp)
    mcp.run()


if __name__ == "__main__":
    main()
