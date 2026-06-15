"""anno_mcp — annotation MCP (6 tools)."""
from __future__ import annotations

import hashlib
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

from mrbigr.core import anno  # noqa: E402
from mrbigr.mcp.longjob import JobSpec, register_wait_for_job, start_or_poll  # noqa: E402


def _hash_text(value) -> str:
    return hashlib.sha256(str(value).encode("utf-8")).hexdigest()[:12]


def _job_root() -> Path:
    cwd = Path.cwd()
    if os.access(cwd, os.W_OK):
        return cwd / ".mrbigr_mcp_jobs"
    return Path(tempfile.gettempdir()) / "mrbigr2_mcp_jobs"


def _default_parse_gtf_output(gtf_file) -> Path:
    job_dir = (_job_root() / "parse_gtf" / _hash_text(Path(gtf_file).expanduser())).resolve()
    return job_dir / "parsed_gtf.csv"


def _qtl_annotation_paths(qtl_file, annotation_file):
    """Resolve (qtl_csv, output_csv, output_dir, output_name) for qtl_annotation.

    When qtl_file is a path, the annotated CSV lands next to it as
    ``<stem>.annotated.csv``. When it's an in-memory dict/DataFrame, the
    inputs are serialized into a content-hashed job dir under ``_job_root()``.
    """
    if isinstance(qtl_file, (str, Path)):
        qtl_path = Path(qtl_file).expanduser().resolve()
        out_csv = qtl_path.parent / f"{qtl_path.stem}.annotated.csv"
        return str(qtl_path), str(out_csv), str(out_csv.parent), out_csv.stem

    import pandas as pd

    job_dir = (_job_root() / "qtl_annotation" / _hash_text(f"{qtl_file!r}|{annotation_file!r}")).resolve()
    job_dir.mkdir(parents=True, exist_ok=True)
    qtl_df = pd.DataFrame(qtl_file) if isinstance(qtl_file, dict) else qtl_file
    if not isinstance(qtl_df, pd.DataFrame):
        qtl_df = pd.DataFrame(qtl_df)
    qtl_csv = job_dir / "qtl.csv"
    qtl_df.to_csv(qtl_csv, index=False)
    out_csv = job_dir / "qtl_annotated.csv"
    return str(qtl_csv), str(out_csv), str(job_dir), "qtl_annotation"


def register(mcp) -> None:  # type: ignore[no-untyped-def]
    @mcp.tool()
    def parse_gtf(gtf_file, output_file=None, return_records=False):
        """Parse GTF annotation file.

        Default return is a compact summary. Set return_records=true only for
        small annotations where a full MCP JSON table is practical. Long-running
        mode writes parsed records to output_file and returns a job status.
        """
        out_file = Path(output_file).expanduser().resolve() if output_file else _default_parse_gtf_output(gtf_file)
        spec = JobSpec(
            kind="parse_gtf",
            runner={
                "type": "python",
                "module": "mrbigr.mcp.runners",
                "function": "parse_gtf_runner",
                "kwargs": {
                    "gtf_file": gtf_file,
                    "output_file": str(out_file),
                },
            },
            output_dir=str(out_file.parent),
            output_name=out_file.stem,
            expected_files=[str(out_file)],
            completion_marker={"type": "file_nonempty"},
            eta_seconds=300,
        )
        response = start_or_poll(spec)
        response["gtf_file"] = gtf_file
        response["output_file"] = str(out_file)
        if response.get("status") == "completed":
            import pandas as pd
            df = pd.read_csv(out_file)
            gene_count = int((df["feature"] == "gene").sum()) if "feature" in df.columns else 0
            response.update({
                "n_records": int(len(df)),
                "n_genes": gene_count,
                "columns": list(df.columns),
            })
            if str(return_records).lower() in {"1", "true", "yes"}:
                response["records"] = df.to_dict()
        return response

    @mcp.tool()
    def annotate_snps(snp_df, gtf_file, window=5000):
        """Annotate SNPs using GTF file."""
        window = int(window)
        result = anno.annotate_snps_simple(snp_df, gtf_file, window=window)
        return result.to_dict() if result is not None else None

    @mcp.tool()
    def qtl_annotation(qtl_file, annotation_file):
        """Annotate QTL regions with genes.

        Long-running: first call returns status='running'; re-call with the
        same arguments to poll, or use wait_for_job(job_file). Result is
        written to ``<qtl_file_stem>.annotated.csv`` (or to a hashed job dir
        when qtl_file is in-memory).
        """
        qtl_csv, out_csv, output_dir, output_name = _qtl_annotation_paths(qtl_file, annotation_file)
        spec = JobSpec(
            kind="qtl_annotation",
            runner={
                "type": "python",
                "module": "mrbigr.mcp.runners",
                "function": "qtl_annotation_runner",
                "kwargs": {
                    "qtl_file": qtl_csv,
                    "annotation_file": str(annotation_file),
                    "output_csv": out_csv,
                },
            },
            output_dir=output_dir,
            output_name=output_name,
            expected_files=[out_csv],
            completion_marker={"type": "file_nonempty"},
            eta_seconds=300,
        )
        response = start_or_poll(spec)
        response["output_file"] = out_csv
        if response.get("status") == "completed" and Path(out_csv).is_file():
            import pandas as pd
            response["result"] = pd.read_csv(out_csv).to_dict()
        return response

    @mcp.tool()
    def create_annotation_db(gtf_file, output_prefix):
        """Create local annotation database from GTF."""
        output_prefix = str(Path(output_prefix).expanduser().resolve())
        db_file = f"{output_prefix}_db.pkl"
        spec = JobSpec(
            kind="create_annotation_db",
            runner={
                "type": "python",
                "module": "mrbigr.mcp.runners",
                "function": "core_call_runner",
                "kwargs": {
                    "module": "mrbigr.core.anno",
                    "function": "create_annotation_db",
                    "kwargs": {
                        "gtf_file": gtf_file,
                        "output_prefix": output_prefix,
                    },
                },
            },
            output_dir=str(Path(output_prefix).parent),
            output_name=Path(output_prefix).name,
            expected_files=[db_file],
            completion_marker={"type": "file_nonempty"},
            eta_seconds=300,
        )
        response = start_or_poll(spec)
        if response.get("status") == "completed":
            response["output_file"] = db_file
        return response

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

    register_wait_for_job(mcp)


def main() -> None:
    from fastmcp import FastMCP
    mcp = FastMCP("anno_mcp")
    register(mcp)
    mcp.run()


if __name__ == "__main__":
    main()
