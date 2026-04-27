"""Tests for deploying MRBIGR2 workflow skills into Agent Skills directories."""
from __future__ import annotations

import os
import subprocess
import sys
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parent.parent
SRC = REPO_ROOT / "src"
sys.path.insert(0, str(SRC))

from mrbigr.skill.deployer import deploy_all_skills, deploy_skill, resolve_targets, uninstall_all_skills  # noqa: E402


ACTIVE_SKILLS = {
    "gwas-pipeline",
    "qtl-to-target",
    "causal-network",
    "functional-enrichment",
    "reproduce-v1-case",
}


def test_deploy_skill_uses_agent_skills_directory_shape(tmp_path: Path) -> None:
    target = deploy_skill("gwas_pipeline", target_dir=tmp_path)

    assert target == tmp_path / "gwas-pipeline" / "SKILL.md"
    assert target.is_file()
    assert "name: gwas-pipeline" in target.read_text(encoding="utf-8")


def test_deploy_skill_accepts_hyphen_skill_name(tmp_path: Path) -> None:
    target = deploy_skill("qtl-to-target", target_dir=tmp_path)

    assert target == tmp_path / "qtl-to-target" / "SKILL.md"
    assert target.is_file()


def test_deploy_all_skips_deprecated_by_default(tmp_path: Path) -> None:
    targets = deploy_all_skills(target_dir=tmp_path)
    installed = {path.parent.name for path in targets}

    assert installed == ACTIVE_SKILLS
    assert not (tmp_path / "mr-analysis-deprecated").exists()
    assert not (tmp_path / "qtl-mapping-deprecated").exists()


def test_resolve_targets_defaults_to_claude() -> None:
    targets = resolve_targets()

    assert len(targets) == 1
    assert targets[0].name == "claude"
    assert targets[0].path == Path.home() / ".claude" / "skills"


def test_resolve_targets_all_includes_default_client_targets(monkeypatch) -> None:
    monkeypatch.setenv("CODEX_HOME", "/tmp/mrbigr-codex-home")

    targets = resolve_targets(targets=["all"])
    by_name = {target.name: target.path for target in targets}

    assert set(by_name) == {"claude", "codex", "gemini", "opencode"}
    assert by_name["codex"] == Path("/tmp/mrbigr-codex-home/skills")
    assert by_name["gemini"] == Path.home() / ".gemini" / "skills"


def test_resolve_targets_agents_is_explicit() -> None:
    targets = resolve_targets(targets=["agents"])

    assert len(targets) == 1
    assert targets[0].name == "agents"
    assert targets[0].path == Path.home() / ".agents" / "skills"


def test_resolve_targets_dedupes_custom_dir() -> None:
    claude_dir = Path.home() / ".claude" / "skills"

    targets = resolve_targets(targets=["claude"], target_dirs=[claude_dir])

    assert len(targets) == 1
    assert targets[0].name == "claude"


def test_resolve_targets_rejects_unknown_alias() -> None:
    try:
        resolve_targets(targets=["unknown-agent"])
    except ValueError as e:
        assert "unknown skill target" in str(e)
        assert "--target-dir" in str(e)
    else:
        raise AssertionError("expected unknown target to fail")


def test_deployer_import_does_not_load_heavy_core_package() -> None:
    script = (
        "import sys; "
        "sys.path.insert(0, 'src'); "
        "import mrbigr.skill.deployer; "
        "print('mrbigr.core' in sys.modules)"
    )
    result = subprocess.run(
        [sys.executable, "-c", script],
        cwd=REPO_ROOT,
        capture_output=True,
        text=True,
        timeout=30,
    )

    assert result.returncode == 0, result.stderr
    assert result.stdout.strip() == "False"


def test_skill_cli_install_all_supports_multiple_target_dirs(tmp_path: Path) -> None:
    env = os.environ.copy()
    env["PYTHONPATH"] = str(SRC)
    env["MRBIGR_ROOT"] = str(REPO_ROOT)
    target_a = tmp_path / "agent-a"
    target_b = tmp_path / "agent-b"
    result = subprocess.run(
        [
            sys.executable,
            "-m",
            "mrbigr.skill_cli",
            "install-all",
            "--target-dir",
            str(target_a),
            "--target-dir",
            str(target_b),
        ],
        cwd=REPO_ROOT,
        env=env,
        capture_output=True,
        text=True,
        timeout=30,
    )

    assert result.returncode == 0, result.stderr
    assert {path.name for path in target_a.iterdir()} == ACTIVE_SKILLS
    assert {path.name for path in target_b.iterdir()} == ACTIVE_SKILLS


def test_uninstall_all_removes_active_skills_only(tmp_path: Path) -> None:
    deploy_all_skills(target_dir=tmp_path)
    third_party = tmp_path / "third-party" / "SKILL.md"
    third_party.parent.mkdir()
    third_party.write_text("---\nname: third-party\n---\n", encoding="utf-8")

    removed = uninstall_all_skills(target_dir=tmp_path)

    assert {info.skill_name for info, was_removed in removed if was_removed} == ACTIVE_SKILLS
    assert third_party.is_file()
    for skill_name in ACTIVE_SKILLS:
        assert not (tmp_path / skill_name).exists()


def test_skill_cli_list_defaults_to_default_client_targets() -> None:
    env = os.environ.copy()
    env["PYTHONPATH"] = str(SRC)
    env["MRBIGR_ROOT"] = str(REPO_ROOT)
    result = subprocess.run(
        [sys.executable, "-m", "mrbigr.skill_cli", "list"],
        cwd=REPO_ROOT,
        env=env,
        capture_output=True,
        text=True,
        timeout=30,
    )

    assert result.returncode == 0, result.stderr
    assert "  claude:" in result.stdout
    assert "  codex:" in result.stdout
    assert "  gemini:" in result.stdout
    assert "  opencode:" in result.stdout
    assert "  agents:" not in result.stdout


def test_skill_cli_uninstall_all_supports_one_target_dir(tmp_path: Path) -> None:
    deploy_all_skills(target_dir=tmp_path)
    env = os.environ.copy()
    env["PYTHONPATH"] = str(SRC)
    env["MRBIGR_ROOT"] = str(REPO_ROOT)
    result = subprocess.run(
        [
            sys.executable,
            "-m",
            "mrbigr.skill_cli",
            "uninstall-all",
            "--target-dir",
            str(tmp_path),
        ],
        cwd=REPO_ROOT,
        env=env,
        capture_output=True,
        text=True,
        timeout=30,
    )

    assert result.returncode == 0, result.stderr
    assert "gwas-pipeline" in result.stdout
    assert not any((tmp_path / skill_name).exists() for skill_name in ACTIVE_SKILLS)
