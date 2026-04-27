"""Input coercion tests for the phenotype MCP wrapper."""
from __future__ import annotations

import importlib.util
import sys
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parent.parent
SERVER_PY = REPO_ROOT / "tool-mcps" / "pheno_mcp" / "src" / "server.py"


def _load_pheno_server():
    spec = importlib.util.spec_from_file_location("pheno_mcp_test_server", SERVER_PY)
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


class FakeMCP:
    def __init__(self) -> None:
        self.tools = {}

    def tool(self):
        def decorate(fn):
            self.tools[fn.__name__] = fn
            return fn

        return decorate


def test_phenotype_frame_reads_csv_path_with_first_column_as_index(tmp_path: Path) -> None:
    server = _load_pheno_server()
    phe_csv = tmp_path / "pheno.csv"
    phe_csv.write_text("ID,trait_a,trait_b\ns1,1.0,2.0\ns2,3.0,4.0\n", encoding="utf-8")

    df = server._phenotype_frame(str(phe_csv))

    assert list(df.index) == ["s1", "s2"]
    assert list(df.columns) == ["trait_a", "trait_b"]


def test_phenotype_frame_sets_explicit_sample_id_column() -> None:
    server = _load_pheno_server()

    df = server._phenotype_frame({"ID": ["s1", "s2"], "trait": [1.0, 2.0]})

    assert list(df.index) == ["s1", "s2"]
    assert list(df.columns) == ["trait"]


def test_filter_missing_accepts_csv_path(tmp_path: Path) -> None:
    server = _load_pheno_server()
    phe_csv = tmp_path / "pheno.csv"
    phe_csv.write_text("ID,trait\ns1,1.0\ns2,2.0\n", encoding="utf-8")
    fake_mcp = FakeMCP()
    server.register(fake_mcp)
    captured = {}
    original = server.pheno.missing_filter

    def fake_missing_filter(d, missing_ratio):
        captured["d"] = d
        captured["missing_ratio"] = missing_ratio
        return d

    try:
        server.pheno.missing_filter = fake_missing_filter
        result = fake_mcp.tools["filter_missing"](str(phe_csv), 0.1)
    finally:
        server.pheno.missing_filter = original

    assert result == {"trait": {"s1": 1.0, "s2": 2.0}}
    assert list(captured["d"].index) == ["s1", "s2"]
    assert list(captured["d"].columns) == ["trait"]
    assert captured["missing_ratio"] == 0.1


def test_select_phenotype_columns_resolves_100gw_alias(tmp_path: Path) -> None:
    server = _load_pheno_server()
    phe_csv = tmp_path / "pheno.csv"
    phe_csv.write_text(
        "ID,Plantheight,100grainweight\ns1,10.0,1.0\ns2,20.0,2.0\n",
        encoding="utf-8",
    )
    fake_mcp = FakeMCP()
    server.register(fake_mcp)

    result = fake_mcp.tools["select_phenotype_columns"](str(phe_csv), "100GW")

    assert result == {"100grainweight": {"s1": 1.0, "s2": 2.0}}


def test_select_phenotype_columns_can_write_csv(tmp_path: Path) -> None:
    server = _load_pheno_server()
    phe_csv = tmp_path / "pheno.csv"
    out_csv = tmp_path / "selected.csv"
    phe_csv.write_text("ID,trait_a,trait_b\ns1,1.0,2.0\ns2,3.0,4.0\n", encoding="utf-8")
    fake_mcp = FakeMCP()
    server.register(fake_mcp)

    result = fake_mcp.tools["select_phenotype_columns"](str(phe_csv), '["trait_b"]', str(out_csv))

    assert result == {"output_file": str(out_csv), "rows": 2, "columns": ["trait_b"]}
    assert out_csv.read_text(encoding="utf-8") == "ID,trait_b\ns1,2.0\ns2,4.0\n"
