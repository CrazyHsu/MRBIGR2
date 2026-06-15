"""Longjob wrapper tests for phase2 MCP servers."""
from __future__ import annotations

import importlib.util
import sys
from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parent.parent


class FakeMCP:
    def __init__(self) -> None:
        self.tools = {}

    def tool(self):
        def decorate(fn):
            self.tools[fn.__name__] = fn
            return fn

        return decorate


def _load_server(name: str):
    server_py = REPO_ROOT / "tool-mcps" / name / "src" / "server.py"
    spec = importlib.util.spec_from_file_location(f"{name}_phase2_test_server", server_py)
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


def _capture_start_or_poll(monkeypatch, module):
    captured = {}

    def fake_start(spec):
        captured["spec"] = spec
        return {
            "status": "running",
            "kind": spec.kind,
            "job_file": str(Path(spec.output_dir) / f"{spec.output_name}.{spec.kind}.mcp_job.json"),
            "expected_files": list(spec.expected_files),
        }

    monkeypatch.setattr(module, "start_or_poll", fake_start)
    return captured


def test_gwas_clump_uses_longjob_and_expected_clumped_files(tmp_path: Path, monkeypatch) -> None:
    server = _load_server("gwas_mcp")
    captured = _capture_start_or_poll(monkeypatch, server)
    fake_mcp = FakeMCP()
    server.register(fake_mcp)
    clump_input = tmp_path / "clump_input"
    result_dir = tmp_path / "clump_result"
    clump_input.mkdir()
    (clump_input / "trait_a.assoc").write_text("SNP\tP\ns1\t0.1\n", encoding="utf-8")

    result = fake_mcp.tools["run_gwas_clump"](
        "geno",
        clump_input_dir=str(clump_input),
        result_dir=str(result_dir),
    )

    spec = captured["spec"]
    assert result["status"] == "running"
    assert spec.kind == "gwas_clump"
    assert spec.runner["function"] == "gwas_clump_runner"
    assert str(result_dir / "gwas_clump.mcp_result.json") in spec.expected_files
    assert str(result_dir / "trait_a.clumped") in spec.expected_files


def test_qtl_phase2_tools_use_longjob(tmp_path: Path, monkeypatch) -> None:
    server = _load_server("qtl_mcp")
    captured = _capture_start_or_poll(monkeypatch, server)
    fake_mcp = FakeMCP()
    server.register(fake_mcp)

    prefix = tmp_path / "qtl_extract"
    fake_mcp.tools["extract_qtl_genotypes"](
        "geno",
        {"CHR": [1], "qtl_start": [10], "qtl_end": [20]},
        str(prefix),
    )
    spec = captured["spec"]
    assert spec.kind == "extract_qtl_genotypes"
    assert spec.runner["function"] == "qtl_extract_runner"
    assert str(tmp_path / "qtl_extract_qtl_0.bed") in spec.expected_files
    assert str(tmp_path / "qtl_extract.mcp_result.json") in spec.expected_files

    fake_mcp.tools["calculate_qtl_haplo"]("geno", ["snp1", "snp2"], str(tmp_path / "hap"))
    spec = captured["spec"]
    assert spec.kind == "calculate_qtl_haplo"
    assert spec.runner["function"] == "qtl_haplo_runner"
    assert spec.expected_files == [str(tmp_path / "hap_haplo.csv")]


def test_net_anno_peak_go_servers_register_longjob_tools(tmp_path: Path, monkeypatch) -> None:
    net_server = _load_server("net_mcp")
    net_capture = _capture_start_or_poll(monkeypatch, net_server)
    net_mcp = FakeMCP()
    net_server.register(net_mcp)
    net_mcp.tools["module_identify"]({"row": ["a"], "col": ["b"], "weight": [0.9]})
    assert net_capture["spec"].kind == "module_identify"
    assert net_capture["spec"].runner["function"] == "net_module_identify_runner"
    assert "wait_for_job" in net_mcp.tools

    anno_server = _load_server("anno_mcp")
    anno_capture = _capture_start_or_poll(monkeypatch, anno_server)
    anno_mcp = FakeMCP()
    anno_server.register(anno_mcp)
    anno_mcp.tools["parse_gtf"]("genes.gtf", output_file=str(tmp_path / "genes.csv"))
    assert anno_capture["spec"].kind == "parse_gtf"
    assert anno_capture["spec"].runner["function"] == "parse_gtf_runner"
    anno_mcp.tools["create_annotation_db"]("genes.gtf", str(tmp_path / "anno"))
    assert anno_capture["spec"].kind == "create_annotation_db"
    assert str(tmp_path / "anno_db.pkl") in anno_capture["spec"].expected_files
    assert "wait_for_job" in anno_mcp.tools

    peak_server = _load_server("peak_mcp")
    peak_capture = _capture_start_or_poll(monkeypatch, peak_server)
    peak_mcp = FakeMCP()
    peak_server.register(peak_mcp)
    peak_mcp.tools["plot_qtl_boxplot"](
        "pheno.csv",
        "geno",
        {"SNP": ["snp1"], "phe_name": ["trait"]},
        output_dir=str(tmp_path / "boxplot"),
    )
    assert peak_capture["spec"].kind == "plot_qtl_boxplot"
    assert peak_capture["spec"].runner["function"] == "peak_boxplot_runner"
    peak_mcp.tools["multi_trait_qtl_plot"]("gwas", "qtl.csv", str(tmp_path / "multi"))
    assert peak_capture["spec"].kind == "multi_trait_qtl_plot"
    assert str(tmp_path / "multi_multi_trait.png") in peak_capture["spec"].expected_files
    assert "wait_for_job" in peak_mcp.tools

    go_server = _load_server("go_mcp")
    go_capture = _capture_start_or_poll(monkeypatch, go_server)
    go_mcp = FakeMCP()
    go_server.register(go_mcp)
    go_mcp.tools["run_go_enrichment"](["gene1"], gene_sets_file="sets.gmt")
    assert go_capture["spec"].kind == "run_go_enrichment"
    assert go_capture["spec"].runner["function"] == "go_enrich_runner"
    go_mcp.tools["run_gsea_analysis"]({"gene": ["gene1"], "score": [1.0]}, gene_sets_file="sets.gmt")
    assert go_capture["spec"].kind == "run_gsea_analysis"
    assert go_capture["spec"].runner["function"] == "gsea_enrich_runner"
    go_mcp.tools["run_kegg_enrichment"](["gene1"])
    assert go_capture["spec"].kind == "run_kegg_enrichment"
    assert go_capture["spec"].runner["function"] == "kegg_enrich_runner"
    go_mcp.tools["extract_go_from_gtf"]("genes.gtf", output_file=str(tmp_path / "go.tsv"))
    assert go_capture["spec"].kind == "extract_go_from_gtf"
    assert str(tmp_path / "go.tsv") in go_capture["spec"].expected_files
    assert "wait_for_job" in go_mcp.tools
