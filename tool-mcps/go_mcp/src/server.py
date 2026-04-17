"""go_mcp — GO / KEGG / GSEA enrichment MCP (9 tools)."""
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

from mrbigr.core import go  # noqa: E402


def register(mcp) -> None:  # type: ignore[no-untyped-def]
    @mcp.tool()
    def run_go_enrichment(gene_list, organism='Mouse', pvalue_cutoff=0.05, qvalue_cutoff=0.05,
                          mode='local', gene_sets_file=None, background_genes=None,
                          go_obo_file=None, obo_cache_dir=None, auto_download_obo=True,
                          drop_unmapped_terms=True):
        """GO enrichment analysis (BP, MF, CC)."""
        import pandas as pd
        genes = pd.DataFrame({'genes': gene_list}) if isinstance(gene_list, list) else gene_list
        result = go.go_enrich(
            genes, organism=organism, pvalue_cutoff=pvalue_cutoff,
            qvalue_cutoff=qvalue_cutoff, mode=mode, gene_sets_file=gene_sets_file,
            background_genes=background_genes, go_obo_file=go_obo_file,
            obo_cache_dir=obo_cache_dir, auto_download_obo=auto_download_obo,
            drop_unmapped_terms=drop_unmapped_terms,
        )
        return result.to_dict() if result is not None else None

    @mcp.tool()
    def run_gsea_analysis(gene_rank, gene_symbol_col='gene', score_col='score', organism='Mouse',
                          mode='local', gene_sets_file=None):
        """Gene Set Enrichment Analysis (GSEA)."""
        import pandas as pd
        df = pd.DataFrame(gene_rank) if isinstance(gene_rank, dict) else gene_rank
        result = go.gsea_enrich(
            None, gene_rank=df, gene_symbol_col=gene_symbol_col, score_col=score_col,
            organism=organism, mode=mode, gene_sets_file=gene_sets_file,
        )
        return result.to_dict() if result is not None else None

    @mcp.tool()
    def plot_go_enrichment(enrich_result, output_prefix, plot_types=None, top_n=20, figsize=None,
                           split_ontology=True, label_col='Term',
                           bar_x='GeneRatio', bar_color='Adjusted P-value',
                           dot_x='GeneRatio', dot_color='Adjusted P-value',
                           dot_size='Count', sort_by=None, ascending=None):
        """Generate GO enrichment plots."""
        import pandas as pd
        df = pd.DataFrame(enrich_result) if isinstance(enrich_result, dict) else enrich_result
        if plot_types is None:
            plot_types = ['barplot', 'dotplot']
        if figsize is None:
            figsize = (8, 6)
        return go.go_plot(
            df, output_prefix, plot_types=plot_types, top_n=top_n, figsize=tuple(figsize),
            split_ontology=split_ontology, label_col=label_col,
            bar_x=bar_x, bar_color=bar_color, dot_x=dot_x, dot_color=dot_color,
            dot_size=dot_size, sort_by=sort_by, ascending=ascending,
        )

    @mcp.tool()
    def plot_gsea_results(gsea_result, output_prefix, top_n=20, figsize=None, format='png',
                          plot_mode='curve', terms=None, curve_terms=1, trace_terms=3, rank_metric=None):
        """Generate GSEA plots (curve, trace, NES summary)."""
        import pandas as pd
        df = pd.DataFrame(gsea_result) if isinstance(gsea_result, dict) else gsea_result
        if figsize is None:
            figsize = (8, 6)
        return go.gsea_plot(
            df, output_prefix, top_n=top_n, figsize=tuple(figsize), format=format,
            plot_mode=plot_mode, terms=terms, curve_terms=curve_terms,
            trace_terms=trace_terms, rank_metric=rank_metric,
        )

    @mcp.tool()
    def run_kegg_enrichment(gene_list, organism='mmu', pvalue_cutoff=0.05):
        """KEGG pathway enrichment analysis."""
        import pandas as pd
        genes = pd.DataFrame({'genes': gene_list}) if isinstance(gene_list, list) else gene_list
        result = go.kegg_enrich(genes, organism=organism, pvalue_cutoff=pvalue_cutoff)
        return result.to_dict() if result is not None else None

    @mcp.tool()
    def get_enrichr_libraries(organism='Mouse'):
        """Get available gene set libraries."""
        return go.enrichr_library_list(organism=organism)

    @mcp.tool()
    def extract_go_from_gtf(gtf_file, output_file=None):
        """Extract GO annotations from GTF file."""
        result = go.go_annotation_from_gtf(gtf_file, output_file=output_file)
        return result.to_dict() if result is not None else None

    @mcp.tool()
    def simplify_go_results(enrich_result, similarity_cutoff=0.7):
        """Simplify GO terms by removing redundancy."""
        import pandas as pd
        df = pd.DataFrame(enrich_result) if isinstance(enrich_result, dict) else enrich_result
        result = go.simplify_go_terms(df, similarity_cutoff=similarity_cutoff)
        return result.to_dict() if result is not None else None

    @mcp.tool()
    def export_go_report(enrich_result, output_file, format='excel'):
        """Export GO enrichment results."""
        import pandas as pd
        df = pd.DataFrame(enrich_result) if isinstance(enrich_result, dict) else enrich_result
        return go.export_go_report(df, output_file, format=format)


def main() -> None:
    from fastmcp import FastMCP
    mcp = FastMCP("go_mcp")
    register(mcp)
    mcp.run()


if __name__ == "__main__":
    main()
