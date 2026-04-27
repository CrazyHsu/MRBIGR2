"""Input coercion tests for the GWAS MCP wrapper."""
from __future__ import annotations

import importlib.util
import sys
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parent.parent
SERVER_PY = REPO_ROOT / "tool-mcps" / "gwas_mcp" / "src" / "server.py"


def _load_gwas_server():
    spec = importlib.util.spec_from_file_location("gwas_mcp_test_server", SERVER_PY)
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
    server = _load_gwas_server()
    phe_csv = tmp_path / "pheno.csv"
    phe_csv.write_text("ID,trait_a,trait_b\ns1,1.0,2.0\ns2,3.0,4.0\n", encoding="utf-8")

    df = server._phenotype_frame(str(phe_csv))

    assert list(df.index) == ["s1", "s2"]
    assert list(df.columns) == ["trait_a", "trait_b"]


def test_phenotype_frame_sets_explicit_sample_id_column() -> None:
    server = _load_gwas_server()

    df = server._phenotype_frame({"ID": ["s1", "s2"], "trait": [1.0, 2.0]})

    assert list(df.index) == ["s1", "s2"]
    assert list(df.columns) == ["trait"]


def test_run_gwas_lmm_accepts_csv_path_and_traits(tmp_path: Path) -> None:
    server = _load_gwas_server()
    phe_csv = tmp_path / "pheno.csv"
    out_dir = tmp_path / "gwas"
    out_dir.mkdir()
    result_file = out_dir / "100grainweight.assoc.txt"
    result_file.write_text("rs\tp_wald\nsnp1\t0.1\n", encoding="utf-8")
    (out_dir / "100grainweight.log.txt").write_text("total computation time\n", encoding="utf-8")
    phe_csv.write_text(
        "ID,Plantheight,100grainweight\ns1,10.0,1.0\ns2,20.0,2.0\n",
        encoding="utf-8",
    )
    fake_mcp = FakeMCP()
    server.register(fake_mcp)

    result = fake_mcp.tools["run_gwas_lmm"](
        str(phe_csv),
        "geno_prefix",
        output_dir=str(out_dir),
        traits="100GW",
    )

    assert result == {"status": "completed", "gwas_results": [str(result_file)]}


def test_select_phenotype_frame_resolves_100gw_alias(tmp_path: Path) -> None:
    server = _load_gwas_server()
    phe_csv = tmp_path / "pheno.csv"
    phe_csv.write_text(
        "ID,Plantheight,100grainweight\ns1,10.0,1.0\ns2,20.0,2.0\n",
        encoding="utf-8",
    )

    df = server._select_phenotype_frame(str(phe_csv), traits="100GW")

    assert list(df.index) == ["s1", "s2"]
    assert list(df.columns) == ["100grainweight"]
