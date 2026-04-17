"""qtl_mcp — QTL analysis MCP (11 tools)."""
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

from mrbigr.core import qtl  # noqa: E402


def register(mcp) -> None:  # type: ignore[no-untyped-def]
    @mcp.tool()
    def detect_qtl_regions(gwas_file, p1=1e-7, p2=1e-5, p2n=5, window=500000):
        """Detect QTL regions from GWAS results."""
        result = qtl.detect_qtl(gwas_file, p1=p1, p2=p2, p2n=p2n, window=window)
        return result.to_dict() if result is not None else None

    @mcp.tool()
    def get_lead_snp(gwas_file, region_chr, region_start, region_end):
        """Get lead SNP in a region."""
        result = qtl.get_lead_snp(gwas_file, region_chr, region_start, region_end)
        return result.to_dict() if result is not None else None

    @mcp.tool()
    def identify_peak_snps(gwas_dir, p_threshold=1e-5, max_peaks=50):
        """Identify peak SNPs from GWAS results."""
        result = qtl.identify_peak_snps(gwas_dir, p_threshold=p_threshold, max_peaks=max_peaks)
        return result.to_dict() if result is not None else None

    @mcp.tool()
    def extract_qtl_genotypes(geno_prefix, qtl_regions, output_prefix):
        """Extract genotypes for QTL regions."""
        import pandas as pd
        qtl_data = pd.DataFrame(qtl_regions) if isinstance(qtl_regions, dict) else qtl_regions
        return qtl.extract_qtl_genotypes(geno_prefix, qtl_data, output_prefix)

    @mcp.tool()
    def calculate_qtl_haplo(geno_prefix, qtl_snp_list, output_prefix):
        """Calculate haplotypes for SNPs in a QTL region."""
        result = qtl.calculate_qtl_haplo(geno_prefix, qtl_snp_list, output_prefix)
        return result.to_dict() if result is not None else None

    @mcp.tool()
    def plot_qtl_region(gwas_file, chr_val, start, end, output_file=None):
        """Plot GWAS results for QTL region."""
        return qtl.plot_qtl_region(gwas_file, chr_val, start, end, output_file=output_file)

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
