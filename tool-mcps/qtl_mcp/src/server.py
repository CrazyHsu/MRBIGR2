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
from mrbigr.mcp.longjob import JobSpec, register_wait_for_job, start_or_poll  # noqa: E402


def _qtl_extract_result_file(output_prefix):
    return Path(f"{output_prefix}.mcp_result.json")


def _qtl_extract_expected_prefixes(qtl_data, output_prefix):
    return [f"{output_prefix}_qtl_{idx}" for idx, _ in qtl_data.iterrows()]


def _qtl_extract_expected_files(qtl_data, output_prefix):
    expected = [str(_qtl_extract_result_file(output_prefix))]
    for prefix in _qtl_extract_expected_prefixes(qtl_data, output_prefix):
        expected.extend(str(Path(f"{prefix}{suffix}")) for suffix in (".bed", ".bim", ".fam"))
    return expected


def _qtl_extract_log_files(qtl_data, output_prefix):
    return [str(Path(f"{prefix}.log")) for prefix in _qtl_extract_expected_prefixes(qtl_data, output_prefix)]


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
        output_prefix = str(Path(output_prefix).expanduser().resolve())
        Path(output_prefix).parent.mkdir(parents=True, exist_ok=True)
        qtl_csv = Path(f"{output_prefix}.mcp_qtl_regions.csv")
        qtl_data.to_csv(qtl_csv, index=False)
        expected = _qtl_extract_expected_files(qtl_data, output_prefix)
        spec = JobSpec(
            kind="extract_qtl_genotypes",
            runner={
                "type": "python",
                "module": "mrbigr.mcp.runners",
                "function": "qtl_extract_runner",
                "kwargs": {
                    "geno_prefix": geno_prefix,
                    "qtl_regions_csv": str(qtl_csv),
                    "output_prefix": output_prefix,
                    "result_json": str(_qtl_extract_result_file(output_prefix)),
                },
            },
            output_dir=str(Path(output_prefix).parent),
            output_name=Path(output_prefix).name,
            expected_files=expected,
            completion_marker={
                "type": "log_contains",
                "pattern": "End time:",
                "log_files": _qtl_extract_log_files(qtl_data, output_prefix),
            },
            eta_seconds=600,
        )
        result = start_or_poll(spec)
        if result.get("status") == "completed":
            prefixes = _qtl_extract_expected_prefixes(qtl_data, output_prefix)
            result["result"] = prefixes
            result["files"] = prefixes
            result["output_prefix"] = output_prefix
        return result

    @mcp.tool()
    def calculate_qtl_haplo(geno_prefix: str, qtl_snp_list: list[str] | str, output_prefix: str):
        """Calculate haplotypes for SNPs in a QTL region.

        Long-running mode: first call starts PLINK dosage extraction and
        repeated calls poll the stable haplotype CSV.
        """
        output_prefix = str(Path(output_prefix).expanduser().resolve())
        output_csv = f"{output_prefix}_haplo.csv"
        spec = JobSpec(
            kind="calculate_qtl_haplo",
            runner={
                "type": "python",
                "module": "mrbigr.mcp.runners",
                "function": "qtl_haplo_runner",
                "kwargs": {
                    "geno_prefix": geno_prefix,
                    "qtl_snp_list": qtl_snp_list,
                    "output_prefix": output_prefix,
                    "output_csv": output_csv,
                },
            },
            output_dir=str(Path(output_prefix).parent),
            output_name=Path(output_prefix).name,
            expected_files=[output_csv],
            completion_marker={"type": "file_nonempty"},
            eta_seconds=300,
        )
        result = start_or_poll(spec)
        if result.get("status") == "completed":
            result["output_file"] = output_csv
        return result

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
    def map_qtl_to_genes(qtl_df, annotation_file, output_file=None):
        """Map QTL regions to genes.

        If output_file is provided, also persist the mapped table as CSV for
        downstream workflow handoff.
        """
        import pandas as pd
        qtl_data = pd.DataFrame(qtl_df) if isinstance(qtl_df, dict) else qtl_df
        result = qtl.map_qtl_to_genes(qtl_data, annotation_file)
        if result is None:
            return None
        if output_file:
            Path(output_file).parent.mkdir(parents=True, exist_ok=True)
            result.to_csv(output_file, index=False)
        return result.to_dict()

    @mcp.tool()
    def qtl_enrichment_test(qtl_genes, background_genes, gene_sets):
        """Test enrichment of QTL genes against gene sets."""
        from mrbigr.core._argjson import maybe_json_loads
        result = qtl.qtl_enrichment_test(
            maybe_json_loads(qtl_genes), maybe_json_loads(background_genes), maybe_json_loads(gene_sets))
        return result.to_dict() if result is not None else None

    @mcp.tool()
    def get_qtl_overlap(qtl1, qtl2):
        """Find overlapping QTL intervals between two QTL sets."""
        result = qtl.get_qtl_overlap(qtl._coerce_qtl_df(qtl1), qtl._coerce_qtl_df(qtl2))
        return result.to_dict() if result is not None else None

    @mcp.tool()
    def export_qtl_bed(qtl_df, output_file):
        """Export QTL as BED file."""
        import pandas as pd
        qtl_data = pd.DataFrame(qtl_df) if isinstance(qtl_df, dict) else qtl_df
        return qtl.export_qtl_bed(qtl_data, output_file)

    register_wait_for_job(mcp)


def main() -> None:
    from fastmcp import FastMCP
    mcp = FastMCP("qtl_mcp")
    register(mcp)
    mcp.run()


if __name__ == "__main__":
    main()
