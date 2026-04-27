"""Tests for MCP runtime environment generation."""
from __future__ import annotations

import os
import sys
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parent.parent
SRC = REPO_ROOT / "src"
sys.path.insert(0, str(SRC))

from mrbigr.mcp.mcp import MCP  # noqa: E402
from mrbigr.core import paths  # noqa: E402


def test_resolved_env_prepends_bundled_library_dir(tmp_path: Path) -> None:
    bundled_lib = tmp_path / "utils" / "lib"
    bundled_lib.mkdir(parents=True)
    mcp = MCP(
        name="gwas_mcp",
        path="tool-mcps/gwas_mcp",
        env_vars={
            "MRBIGR_ROOT": "${REPO_ROOT}",
            "LD_LIBRARY_PATH": "/opt/custom",
        },
    )

    env = mcp.resolved_env(tmp_path)

    assert env["MRBIGR_ROOT"] == str(tmp_path)
    assert env["LD_LIBRARY_PATH"] == os.pathsep.join([str(bundled_lib), "/opt/custom"])


def test_resolved_env_does_not_duplicate_bundled_library_dir(tmp_path: Path) -> None:
    bundled_lib = tmp_path / "utils" / "lib"
    bundled_lib.mkdir(parents=True)
    mcp = MCP(
        name="gwas_mcp",
        path="tool-mcps/gwas_mcp",
        env_vars={"LD_LIBRARY_PATH": str(bundled_lib)},
    )

    env = mcp.resolved_env(tmp_path)

    assert env["LD_LIBRARY_PATH"] == str(bundled_lib)


def test_ensure_bundled_runtime_libs_updates_process_env(tmp_path: Path, monkeypatch) -> None:
    bundled_lib = tmp_path / "utils" / "lib"
    bundled_lib.mkdir(parents=True)
    monkeypatch.setenv("MRBIGR_ROOT", str(tmp_path))
    monkeypatch.setenv("LD_LIBRARY_PATH", "/opt/custom")

    value = paths.ensure_bundled_runtime_libs()

    assert value == os.pathsep.join([str(bundled_lib), "/opt/custom"])
    assert os.environ["LD_LIBRARY_PATH"] == value
