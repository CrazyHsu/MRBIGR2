"""vis_mcp — visualization MCP (9 tools)."""
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

from mrbigr.core import vis  # noqa: E402


def register(mcp) -> None:  # type: ignore[no-untyped-def]
    @mcp.tool()
    def plot_manhattan(gwas_file, output_file=None, significance=5e-8, suggest=1e-5):
        """Generate Manhattan plot from GWAS results."""
        significance = float(significance)
        suggest = float(suggest)
        return vis.manhattan_plot(gwas_file, output_file=output_file,
                                  significance=significance, suggest=suggest)

    @mcp.tool()
    def plot_qq(gwas_file, output_file=None):
        """Generate Q-Q plot from GWAS results."""
        return vis.qq_plot(gwas_file, output_file=output_file)

    @mcp.tool()
    def plot_pca(pca_file, output_file=None, pc_x='PC1', pc_y='PC2'):
        """Generate PCA scatter plot."""
        return vis.pca_plot(pca_file, output_file=output_file, pc_x=pc_x, pc_y=pc_y)

    @mcp.tool()
    def plot_tsne(tsne_file, output_file=None):
        """Generate t-SNE scatter plot."""
        return vis.tsne_plot(tsne_file, output_file=output_file)

    @mcp.tool()
    def plot_ld_heatmap(geno_prefix, output_file=None, max_snps=500):
        """Generate LD heatmap from genotype data."""
        max_snps = int(max_snps)
        return vis.ld_heatmap(geno_prefix, output_file=output_file, max_snps=max_snps)

    @mcp.tool()
    def plot_phenotype_hist(pheno_file, output_file=None):
        """Generate histogram for phenotype distribution."""
        return vis.phenotype_hist(pheno_file, output_file=output_file)

    @mcp.tool()
    def plot_phenotype_boxplot(pheno_file, output_file=None):
        """Generate boxplot for phenotypes."""
        return vis.phenotype_boxplot(pheno_file, output_file=output_file)

    @mcp.tool()
    def plot_phenotype_correlation(pheno_file, output_file=None, method='pearson'):
        """Generate phenotype correlation heatmap."""
        return vis.phenotype_correlation(pheno_file, output_file=output_file, method=method)

    @mcp.tool()
    def gwas_summary(gwas_dir, output_prefix='gwas_summary'):
        """Generate comprehensive GWAS summary plots (Manhattan + QQ)."""
        return vis.gwas_summary_plot(gwas_dir, output_prefix=output_prefix)


def main() -> None:
    from fastmcp import FastMCP
    mcp = FastMCP("vis_mcp")
    register(mcp)
    mcp.run()


if __name__ == "__main__":
    main()
