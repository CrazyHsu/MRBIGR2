"""Input coercion tests for the GWAS MCP wrapper."""
from __future__ import annotations

import importlib.util
import sys
from pathlib import Path

import pandas as pd

REPO_ROOT = Path(__file__).resolve().parent.parent
SERVER_PY = REPO_ROOT / "tool-mcps" / "gwas_mcp" / "src" / "server.py"
REPO_SRC = REPO_ROOT / "src"
if str(REPO_SRC) not in sys.path:
    sys.path.insert(0, str(REPO_SRC))


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

    assert result["status"] == "completed"
    assert result["kind"] == "gwas_lmm"
    assert result["result_files"] == [str(result_file)]


def test_run_gwas_lmm_wires_external_kinship_file(monkeypatch, tmp_path: Path) -> None:
    server = _load_gwas_server()
    phe_csv = tmp_path / "pheno.csv"
    out_dir = tmp_path / "gwas"
    kinship = tmp_path / "kinship" / "chr_HAMP_qc.cXX.txt"
    kinship.parent.mkdir()
    kinship.write_text("1\n", encoding="utf-8")
    phe_csv.write_text("ID,100grainweight\ns1,1.0\ns2,2.0\n", encoding="utf-8")
    captured = {}

    def fake_start_or_poll(spec):
        captured["spec"] = spec
        return {"status": "running", "kind": spec.kind}

    monkeypatch.setattr(server, "start_or_poll", fake_start_or_poll)
    fake_mcp = FakeMCP()
    server.register(fake_mcp)

    result = fake_mcp.tools["run_gwas_lmm"](
        str(phe_csv),
        "geno_prefix",
        output_dir=str(out_dir),
        traits="100GW",
        kinship_file=str(kinship),
    )

    spec = captured["spec"]
    assert result["status"] == "running"
    assert spec.runner["kwargs"]["kinship_file"] == str(kinship.resolve())
    assert spec.extras["kinship_file"] == str(kinship.resolve())


def test_run_gwas_lmm_promotes_sample_summary(monkeypatch, tmp_path: Path) -> None:
    server = _load_gwas_server()
    phe_csv = tmp_path / "pheno.csv"
    out_dir = tmp_path / "gwas"
    phe_csv.write_text("ID,100grainweight\ns1,1.0\ns2,2.0\n", encoding="utf-8")
    sample_summary = {"traits": {"100grainweight": {"matched_nonmissing": 2}}}

    def fake_start_or_poll(spec):
        return {
            "status": "completed",
            "kind": spec.kind,
            "result_files": spec.expected_files,
            "result": {"result_files": spec.expected_files, "sample_summary": sample_summary},
        }

    monkeypatch.setattr(server, "start_or_poll", fake_start_or_poll)
    fake_mcp = FakeMCP()
    server.register(fake_mcp)

    result = fake_mcp.tools["run_gwas_lmm"](
        str(phe_csv),
        "geno_prefix",
        output_dir=str(out_dir),
        traits="100GW",
    )

    assert result["sample_summary"] == sample_summary


def test_core_gwas_lmm_reuses_external_kinship(monkeypatch, tmp_path: Path) -> None:
    from mrbigr.core import gwas

    geno_prefix = tmp_path / "chr_HAMP_qc"
    for ext in [".bed", ".bim"]:
        (tmp_path / f"chr_HAMP_qc{ext}").write_text("", encoding="utf-8")
    (tmp_path / "chr_HAMP_qc.fam").write_text("0 s1 0 0 0 -9\n", encoding="utf-8")
    out_dir = tmp_path / "gwas"
    kinship = tmp_path / "kinship" / "chr_HAMP_qc.cXX.txt"
    kinship.parent.mkdir()
    kinship.write_text("1\n", encoding="utf-8")
    phe = pd.DataFrame({"100grainweight": [1.0]}, index=["s1"])
    kinship_commands = []
    gwas_commands = []

    def fake_subprocess_run(cmd, *args, **kwargs):
        kinship_commands.append(cmd)

        class Result:
            returncode = 0

        return Result()

    def fake_run_shell_commands(cmds, num_threads, label):
        gwas_commands.extend(cmds)
        return [0 for _ in cmds]

    monkeypatch.setattr(gwas.subprocess, "run", fake_subprocess_run)
    monkeypatch.setattr(gwas, "_run_shell_commands", fake_run_shell_commands)

    results = gwas.gwas_lmm(
        phe,
        str(geno_prefix),
        output_name="100gw",
        output_dir=str(out_dir),
        kinship_file=str(kinship),
    )

    assert kinship_commands == []
    assert results == [str(out_dir.resolve() / "100gw_100grainweight.assoc.txt")]
    assert f"-k {kinship.resolve()}" in gwas_commands[0]
    assert not (out_dir / "chr_HAMP_qc.cXX.txt").exists()


def test_core_gwas_clump_preserves_mcp_job_files(tmp_path: Path) -> None:
    from mrbigr.core import gwas

    clump_input = tmp_path / "clump_input"
    result_dir = tmp_path / "clump_result"
    clump_input.mkdir()
    result_dir.mkdir()
    job_file = result_dir / "gwas_clump.gwas_clump.mcp_job.json"
    stdout_file = result_dir / "gwas_clump.gwas_clump.mcp_stdout.log"
    stale_output = result_dir / "old.clumped"
    job_file.write_text('{"status":"running"}\n', encoding="utf-8")
    stdout_file.write_text("", encoding="utf-8")
    stale_output.write_text("stale\n", encoding="utf-8")

    result = gwas.gwas_clump(
        str(tmp_path / "geno"),
        clump_input_dir=str(clump_input),
        result_dir=str(result_dir),
    )

    assert Path(result) == result_dir
    assert job_file.read_text(encoding="utf-8") == '{"status":"running"}\n'
    assert stdout_file.is_file()
    assert not stale_output.exists()


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


def test_get_top_snps_accepts_string_n(tmp_path: Path) -> None:
    server = _load_gwas_server()
    gwas_file = tmp_path / "assoc.txt"
    gwas_file.write_text(
        "rs\tp_wald\nsnp1\t0.5\nsnp2\t0.001\nsnp3\t0.01\n",
        encoding="utf-8",
    )
    fake_mcp = FakeMCP()
    server.register(fake_mcp)

    result = fake_mcp.tools["get_top_snps"](str(gwas_file), n="2")

    assert list(result["rs"].values()) == ["snp2", "snp3"]


def test_gwas_lmm_runner_auto_plot_generates_plots(monkeypatch, tmp_path: Path) -> None:
    from mrbigr.core import gwas, vis
    from mrbigr.mcp import runners

    phe_csv = tmp_path / "pheno.csv"
    phe_csv.write_text("ID,trait\ns1,1.0\ns2,\ns3,3.0\n", encoding="utf-8")
    geno_prefix = tmp_path / "geno"
    (tmp_path / "geno.fam").write_text(
        "0 s1 0 0 0 -9\n0 s2 0 0 0 -9\n0 s4 0 0 0 -9\n",
        encoding="utf-8",
    )
    assoc_file = tmp_path / "trait.assoc.txt"
    assoc_file.write_text("rs\tp_wald\nsnp1\t0.1\n", encoding="utf-8")

    def fake_gwas_lmm(*args, **kwargs):
        return [str(assoc_file)]

    def fake_manhattan_plot(gwas_file, output_file=None):
        Path(output_file).write_text("manhattan\n", encoding="utf-8")
        return output_file

    def fake_qq_plot(gwas_file, output_file=None):
        Path(output_file).write_text("qq\n", encoding="utf-8")
        return output_file

    monkeypatch.setattr(gwas, "gwas_lmm", fake_gwas_lmm)
    monkeypatch.setattr(vis, "manhattan_plot", fake_manhattan_plot)
    monkeypatch.setattr(vis, "qq_plot", fake_qq_plot)

    result = runners.gwas_lmm_runner(
        str(phe_csv),
        str(geno_prefix),
        output_name="trait",
        output_dir=str(tmp_path),
        auto_plot="true",
    )

    assert result["result_files"] == [str(assoc_file)]
    assert Path(result["plots"][0]["manhattan"]).is_file()
    assert Path(result["plots"][0]["qq"]).is_file()
    assert result["sample_summary"]["phenotype_rows"] == 3
    assert result["sample_summary"]["fam_rows"] == 3
    assert result["sample_summary"]["traits"]["trait"]["phenotype_nonmissing"] == 2
    assert result["sample_summary"]["traits"]["trait"]["matched_nonmissing"] == 1
    assert Path(result["sample_summary_file"]).is_file()


def test_calculate_lambda_reads_csv_with_generic_p_column(tmp_path: Path) -> None:
    """Comma-separated GWAS file with a generic 'P' column is accepted
    (auto-detected separator + broadened p-value column recognition)."""
    server = _load_gwas_server()
    gwas_csv = tmp_path / "gwas.csv"
    gwas_csv.write_text(
        "SNP,P\nrs1,0.9\nrs2,0.5\nrs3,0.1\nrs4,0.05\nrs5,0.01\nrs6,0.001\n",
        encoding="utf-8",
    )
    fake_mcp = FakeMCP()
    server.register(fake_mcp)

    result = fake_mcp.tools["calculate_lambda"](str(gwas_csv))

    assert result["n"] == 6
    assert isinstance(result["lambda"], float)


def test_calculate_lambda_reads_tab_assoc_p_wald(tmp_path: Path) -> None:
    """Tab-separated .assoc.txt with p_wald still works (no regression)."""
    server = _load_gwas_server()
    assoc = tmp_path / "lm.assoc.txt"
    assoc.write_text("rs\tp_wald\nsnp1\t0.5\nsnp2\t0.1\nsnp3\t0.01\n", encoding="utf-8")
    fake_mcp = FakeMCP()
    server.register(fake_mcp)

    result = fake_mcp.tools["calculate_lambda"](str(assoc))

    assert result["n"] == 3


def test_calculate_lambda_accepts_inline_array_and_json_string() -> None:
    """Both a real list and a JSON-string array (how this MCP client delivers
    inline args) decode to the same result."""
    server = _load_gwas_server()
    fake_mcp = FakeMCP()
    server.register(fake_mcp)

    as_list = fake_mcp.tools["calculate_lambda"]([0.9, 0.5, 0.1, 0.05, 0.01])
    as_json = fake_mcp.tools["calculate_lambda"]("[0.9, 0.5, 0.1, 0.05, 0.01]")

    assert as_list["n"] == 5
    assert as_json["n"] == 5
    assert as_list["lambda"] == as_json["lambda"]
