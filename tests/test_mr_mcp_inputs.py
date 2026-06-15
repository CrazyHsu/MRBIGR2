"""Input coercion tests for the MR MCP wrapper."""
from __future__ import annotations

import importlib.util
import sys
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parent.parent
SERVER_PY = REPO_ROOT / "tool-mcps" / "mr_mcp" / "src" / "server.py"
REPO_SRC = REPO_ROOT / "src"
if str(REPO_SRC) not in sys.path:
    sys.path.insert(0, str(REPO_SRC))


def _load_mr_server():
    spec = importlib.util.spec_from_file_location("mr_mcp_test_server", SERVER_PY)
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


def test_format_qtl_for_mr_infers_missing_phe_name(tmp_path: Path) -> None:
    qtl_csv = tmp_path / "100grainweight_lmm_100grainweight.qtl.csv"
    out_csv = tmp_path / "formatted.csv"
    qtl_csv.write_text(
        "CHR,qtl_start,qtl_end,SNP,P,qtl_length\n"
        "5,56668617,57003917,chr5.s_56943468,9.44252e-7,335300\n",
        encoding="utf-8",
    )
    server = _load_mr_server()
    fake_mcp = FakeMCP()
    server.register(fake_mcp)

    result = fake_mcp.tools["format_qtl_for_mr"](str(qtl_csv), output_file=str(out_csv))

    assert result["phe_name"][0] == "100grainweight"
    assert out_csv.is_file()


def test_persist_gwas_filters_plink_clumped_instruments(tmp_path: Path) -> None:
    server = _load_mr_server()
    gwas_csv = tmp_path / "assoc.txt"
    clumped = tmp_path / "trait.clumped"
    out_dir = tmp_path / "mr"
    gwas_csv.write_text(
        "rs\tbeta\tse\tp_wald\n"
        "snp1\t0.1\t0.01\t0.001\n"
        "snp2\t0.2\t0.02\t0.002\n"
        "snp3\t0.3\t0.03\t0.003\n",
        encoding="utf-8",
    )
    clumped.write_text(
        " CHR    F     SNP       BP        P    TOTAL\n"
        "   1    1     snp2      20   0.002        3\n"
        "   1    1     snp3      30   0.003        2\n",
        encoding="utf-8",
    )

    result = server._persist_gwas(str(gwas_csv), out_dir, "exposure", instrument_file=str(clumped))

    assert Path(result).read_text(encoding="utf-8").splitlines() == [
        "rs,beta,se,pvalue",
        "snp2,0.2,0.02,0.002",
        "snp3,0.3,0.03,0.003",
    ]


def test_run_mr_analysis_accepts_instrument_file(monkeypatch, tmp_path: Path) -> None:
    server = _load_mr_server()
    exposure = tmp_path / "exposure.tsv"
    outcome = tmp_path / "outcome.tsv"
    clumped = tmp_path / "trait.clumped"
    exposure.write_text(
        "rs\tbeta\tse\nsnp1\t0.1\t0.01\nsnp2\t0.2\t0.02\n",
        encoding="utf-8",
    )
    outcome.write_text(
        "rs\tbeta\tse\nsnp1\t0.3\t0.03\nsnp2\t0.4\t0.04\n",
        encoding="utf-8",
    )
    clumped.write_text("CHR SNP BP P\n1 snp2 20 0.002\n", encoding="utf-8")
    captured = {}

    def fake_start_or_poll(spec):
        captured["spec"] = spec
        return {"status": "running", "kind": spec.kind}

    monkeypatch.setattr(server, "start_or_poll", fake_start_or_poll)
    fake_mcp = FakeMCP()
    server.register(fake_mcp)

    result = fake_mcp.tools["run_mr_analysis"](
        str(exposure),
        str(outcome),
        output_dir=str(tmp_path / "mr"),
        output_name="filtered",
        instrument_file=str(clumped),
    )

    spec = captured["spec"]
    assert result["status"] == "running"
    assert spec.extras["instrument_file"] == str(clumped)
    assert Path(spec.runner["kwargs"]["exposure_gwas_csv"]).read_text(encoding="utf-8").splitlines() == [
        "rs,beta,se",
        "snp2,0.2,0.02",
    ]
    assert Path(spec.runner["kwargs"]["outcome_gwas_csv"]).read_text(encoding="utf-8").splitlines() == [
        "rs,beta,se",
        "snp2,0.4,0.04",
    ]
