"""anno_mcp — annotation MCP (6 tools)."""
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

from mrbigr.core import anno  # noqa: E402


def register(mcp) -> None:  # type: ignore[no-untyped-def]
    @mcp.tool()
    def parse_gtf(gtf_file, output_file=None, return_records=False):
        """Parse GTF annotation file.

        Default return is a compact summary. Set return_records=true only for
        small annotations where a full MCP JSON table is practical.
        """
        result = anno.parse_gtf(gtf_file)
        if result is None:
            return None
        if output_file:
            Path(output_file).parent.mkdir(parents=True, exist_ok=True)
            result.to_csv(output_file, index=False)
        if str(return_records).lower() in {"1", "true", "yes"}:
            return result.to_dict()
        gene_count = 0
        if "feature" in result.columns:
            gene_count = int((result["feature"] == "gene").sum())
        return {
            "gtf_file": gtf_file,
            "output_file": output_file,
            "n_records": int(len(result)),
            "n_genes": gene_count,
            "columns": list(result.columns),
        }

    @mcp.tool()
    def annotate_snps(snp_df, gtf_file, window=5000):
        """Annotate SNPs using GTF file."""
        window = int(window)
        result = anno.annotate_snps_simple(snp_df, gtf_file, window=window)
        return result.to_dict() if result is not None else None

    @mcp.tool()
    def qtl_annotation(qtl_file, annotation_file):
        """Annotate QTL regions with genes."""
        result = anno.qtl_annotation(qtl_file, annotation_file)
        return result.to_dict() if result is not None else None

    @mcp.tool()
    def create_annotation_db(gtf_file, output_prefix):
        """Create local annotation database from GTF."""
        return anno.create_annotation_db(gtf_file, output_prefix)

    @mcp.tool()
    def predict_variant_effect(snp_chr=None, snp_pos=None, ref_allele=None, alt_allele=None,
                               gtf_file=None, simple_ref=None, simple_alt=None, region_type=None):
        """Predict variant functional effects for a single variant."""
        if snp_pos is not None:
            snp_pos = int(snp_pos)
        return anno.predict_variant_effect(
            snp_chr=snp_chr, snp_pos=snp_pos,
            ref_allele=ref_allele, alt_allele=alt_allele,
            gtf_file=gtf_file, simple_ref=simple_ref,
            simple_alt=simple_alt, region_type=region_type,
        )

    @mcp.tool()
    def get_genes_in_region(chr_name, start, end, annotation_file, mode='overlap'):
        """Get all genes in a genomic region."""
        start = int(start)
        end = int(end)
        result = anno.get_genes_in_region(chr_name, start, end, annotation_file, mode=mode)
        return result.to_dict() if result is not None else None


def main() -> None:
    from fastmcp import FastMCP
    mcp = FastMCP("anno_mcp")
    register(mcp)
    mcp.run()


if __name__ == "__main__":
    main()
