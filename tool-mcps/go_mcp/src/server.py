"""go_mcp — GO / KEGG / GSEA enrichment MCP (9 tools)."""
from __future__ import annotations

import hashlib
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

from mrbigr.core import go  # noqa: E402
from mrbigr.mcp.longjob import JobSpec, register_wait_for_job, start_or_poll  # noqa: E402


def _stable_hash(value) -> str:
    try:
        payload = json.dumps(value, ensure_ascii=False, sort_keys=True, default=str)
    except TypeError:
        payload = repr(value)
    return hashlib.sha256(payload.encode("utf-8")).hexdigest()[:12]


def _job_root() -> Path:
    cwd = Path.cwd()
    if os.access(cwd, os.W_OK):
        return cwd / ".mrbigr_mcp_jobs"
    return Path(tempfile.gettempdir()) / "mrbigr2_mcp_jobs"


def _default_output(kind: str, payload, suffix: str = ".csv") -> Path:
    job_dir = (_job_root() / kind / _stable_hash(payload)).resolve()
    return job_dir / f"result{suffix}"


def _as_bool(value) -> bool:
    if isinstance(value, bool):
        return value
    if value is None:
        return False
    return str(value).strip().lower() in {"1", "true", "yes", "y", "on"}


def _attach_table_result(response: dict, output_file: Path, *, sep=",") -> dict:
    response["output_file"] = str(output_file)
    if response.get("status") == "completed" and output_file.is_file():
        import pandas as pd
        try:
            response["result"] = pd.read_csv(output_file, sep=sep).to_dict()
        except pd.errors.EmptyDataError:
            response["result"] = {}
    return response


def register(mcp) -> None:  # type: ignore[no-untyped-def]
    @mcp.tool()
    def run_go_enrichment(gene_list, organism='Mouse', pvalue_cutoff=0.05, qvalue_cutoff=0.05,
                          mode='local', gene_sets_file=None, background_genes=None,
                          go_obo_file=None, obo_cache_dir=None, auto_download_obo=True,
                          drop_unmapped_terms=True, output_file=None):
        """GO enrichment analysis (BP, MF, CC)."""
        payload = {
            "gene_list": gene_list,
            "organism": organism,
            "pvalue_cutoff": pvalue_cutoff,
            "qvalue_cutoff": qvalue_cutoff,
            "mode": mode,
            "gene_sets_file": gene_sets_file,
            "background_genes": background_genes,
            "go_obo_file": go_obo_file,
            "obo_cache_dir": obo_cache_dir,
            "auto_download_obo": auto_download_obo,
            "drop_unmapped_terms": drop_unmapped_terms,
        }
        out_file = Path(output_file).expanduser().resolve() if output_file else _default_output("run_go_enrichment", payload)
        spec = JobSpec(
            kind="run_go_enrichment",
            runner={
                "type": "python",
                "module": "mrbigr.mcp.runners",
                "function": "go_enrich_runner",
                "kwargs": {
                    **payload,
                    "pvalue_cutoff": float(pvalue_cutoff),
                    "qvalue_cutoff": float(qvalue_cutoff),
                    "auto_download_obo": _as_bool(auto_download_obo),
                    "drop_unmapped_terms": _as_bool(drop_unmapped_terms),
                    "output_file": str(out_file),
                },
            },
            output_dir=str(out_file.parent),
            output_name=out_file.stem,
            expected_files=[str(out_file)],
            completion_marker={"type": "file_nonempty"},
            eta_seconds=600,
        )
        return _attach_table_result(start_or_poll(spec), out_file)

    @mcp.tool()
    def run_gsea_analysis(gene_rank, gene_symbol_col='gene', score_col='score', organism='Mouse',
                          mode='local', gene_sets_file=None, output_file=None):
        """Gene Set Enrichment Analysis (GSEA)."""
        payload = {
            "gene_rank": gene_rank,
            "gene_symbol_col": gene_symbol_col,
            "score_col": score_col,
            "organism": organism,
            "mode": mode,
            "gene_sets_file": gene_sets_file,
        }
        out_file = Path(output_file).expanduser().resolve() if output_file else _default_output("run_gsea_analysis", payload)
        spec = JobSpec(
            kind="run_gsea_analysis",
            runner={
                "type": "python",
                "module": "mrbigr.mcp.runners",
                "function": "gsea_enrich_runner",
                "kwargs": {
                    **payload,
                    "output_file": str(out_file),
                },
            },
            output_dir=str(out_file.parent),
            output_name=out_file.stem,
            expected_files=[str(out_file)],
            completion_marker={"type": "file_nonempty"},
            eta_seconds=600,
        )
        return _attach_table_result(start_or_poll(spec), out_file)

    @mcp.tool()
    def plot_go_enrichment(enrich_result, output_prefix, plot_types=None, top_n=20, figsize=None,
                           split_ontology=True, label_col='Term',
                           bar_x='GeneRatio', bar_color='Adjusted P-value',
                           dot_x='GeneRatio', dot_color='Adjusted P-value',
                           dot_size='Count', sort_by=None, ascending=None):
        """Generate GO enrichment plots."""
        from mrbigr.core._argjson import coerce_table
        df = coerce_table(enrich_result)
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
        from mrbigr.core._argjson import coerce_table
        df = coerce_table(gsea_result)
        if figsize is None:
            figsize = (8, 6)
        return go.gsea_plot(
            df, output_prefix, top_n=top_n, figsize=tuple(figsize), format=format,
            plot_mode=plot_mode, terms=terms, curve_terms=curve_terms,
            trace_terms=trace_terms, rank_metric=rank_metric,
        )

    @mcp.tool()
    def run_kegg_enrichment(gene_list, organism='mmu', pvalue_cutoff=0.05, output_file=None):
        """KEGG pathway enrichment analysis."""
        payload = {
            "gene_list": gene_list,
            "organism": organism,
            "pvalue_cutoff": pvalue_cutoff,
        }
        out_file = Path(output_file).expanduser().resolve() if output_file else _default_output("run_kegg_enrichment", payload)
        spec = JobSpec(
            kind="run_kegg_enrichment",
            runner={
                "type": "python",
                "module": "mrbigr.mcp.runners",
                "function": "kegg_enrich_runner",
                "kwargs": {
                    **payload,
                    "pvalue_cutoff": float(pvalue_cutoff),
                    "output_file": str(out_file),
                },
            },
            output_dir=str(out_file.parent),
            output_name=out_file.stem,
            expected_files=[str(out_file)],
            completion_marker={"type": "file_nonempty"},
            eta_seconds=600,
        )
        return _attach_table_result(start_or_poll(spec), out_file)

    @mcp.tool()
    def get_enrichr_libraries(organism='Mouse'):
        """Get available gene set libraries."""
        return go.enrichr_library_list(organism=organism)

    @mcp.tool()
    def extract_go_from_gtf(gtf_file, output_file=None):
        """Extract GO annotations from GTF file."""
        out_file = (
            Path(output_file).expanduser().resolve()
            if output_file
            else _default_output("extract_go_from_gtf", {"gtf_file": gtf_file}, suffix=".tsv")
        )
        spec = JobSpec(
            kind="extract_go_from_gtf",
            runner={
                "type": "python",
                "module": "mrbigr.mcp.runners",
                "function": "core_call_runner",
                "kwargs": {
                    "module": "mrbigr.core.go",
                    "function": "go_annotation_from_gtf",
                    "kwargs": {
                        "gtf_file": gtf_file,
                        "output_file": str(out_file),
                    },
                },
            },
            output_dir=str(out_file.parent),
            output_name=out_file.stem,
            expected_files=[str(out_file)],
            completion_marker={"type": "file_nonempty"},
            eta_seconds=300,
        )
        return _attach_table_result(start_or_poll(spec), out_file, sep="\t")

    @mcp.tool()
    def simplify_go_results(enrich_result, similarity_cutoff=0.7):
        """Simplify GO terms by removing redundancy."""
        from mrbigr.core._argjson import coerce_table
        df = coerce_table(enrich_result)
        result = go.simplify_go_terms(df, similarity_cutoff=similarity_cutoff)
        return result.to_dict() if result is not None else None

    @mcp.tool()
    def export_go_report(enrich_result, output_file, format='excel'):
        """Export GO enrichment results."""
        from mrbigr.core._argjson import coerce_table
        df = coerce_table(enrich_result)
        return go.export_go_report(df, output_file, format=format)

    register_wait_for_job(mcp)


def main() -> None:
    from fastmcp import FastMCP
    mcp = FastMCP("go_mcp")
    register(mcp)
    mcp.run()


if __name__ == "__main__":
    main()
