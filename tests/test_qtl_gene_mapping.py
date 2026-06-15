"""Regression tests for QTL-to-gene mapping."""
from __future__ import annotations

import importlib.util
import sys
from pathlib import Path

import pandas as pd

REPO_ROOT = Path(__file__).resolve().parent.parent
REPO_SRC = REPO_ROOT / "src"
QTL_SERVER = REPO_ROOT / "tool-mcps" / "qtl_mcp" / "src" / "server.py"
if str(REPO_SRC) not in sys.path:
    sys.path.insert(0, str(REPO_SRC))


class FakeMCP:
    def __init__(self) -> None:
        self.tools = {}

    def tool(self):
        def decorate(fn):
            self.tools[fn.__name__] = fn
            return fn

        return decorate


def _load_qtl_server():
    spec = importlib.util.spec_from_file_location("qtl_mcp_test_server", QTL_SERVER)
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


def test_map_qtl_to_genes_counts_unique_gene_rows_from_parsed_gtf_csv(tmp_path: Path) -> None:
    from mrbigr.core import qtl

    annotation = tmp_path / "parsed_gtf.csv"
    annotation.write_text(
        "chr,feature,start,end,strand,gene_id,gene_name,gene_type\n"
        "5,gene,100,200,+,geneA,,protein_coding\n"
        "5,transcript,100,200,+,geneA,,protein_coding\n"
        "5,exon,120,150,+,geneA,,protein_coding\n"
        "5,gene,300,400,-,geneB,,protein_coding\n",
        encoding="utf-8",
    )
    qtl_df = pd.DataFrame({
        "CHR": [5],
        "qtl_start": [90],
        "qtl_end": [210],
        "SNP": ["chr5.s_150"],
        "P": [1e-7],
    })

    result = qtl.map_qtl_to_genes(qtl_df, str(annotation))

    assert result.loc[0, "n_genes"] == 1
    assert result.loc[0, "genes"] == "geneA"


def test_map_qtl_to_genes_mcp_can_persist_csv(tmp_path: Path) -> None:
    server = _load_qtl_server()
    fake_mcp = FakeMCP()
    server.register(fake_mcp)

    annotation = tmp_path / "parsed_gtf.csv"
    annotation.write_text(
        "chr,feature,start,end,strand,gene_id,gene_name,gene_type\n"
        "5,gene,100,200,+,geneA,,protein_coding\n",
        encoding="utf-8",
    )
    output_file = tmp_path / "qtl_genes.csv"

    result = fake_mcp.tools["map_qtl_to_genes"](
        {
            "CHR": {"0": 5},
            "qtl_start": {"0": 90},
            "qtl_end": {"0": 210},
            "SNP": {"0": "chr5.s_150"},
            "P": {"0": 1e-7},
        },
        str(annotation),
        output_file=str(output_file),
    )

    assert output_file.is_file()
    persisted = pd.read_csv(output_file)
    assert result["genes"][0] == "geneA"
    assert persisted.loc[0, "genes"] == "geneA"
