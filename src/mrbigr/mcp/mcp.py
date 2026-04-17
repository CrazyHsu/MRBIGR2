"""MCP dataclass — one entry per row of mcps.yaml.

Each MCP has a four-state lifecycle mirroring ProteinMCP:

    NOT_INSTALLED  server.py missing on disk AND not registered with Claude
    INSTALLED      server.py present          but not yet registered
    REGISTERED     registered with Claude    but server.py missing (broken)
    BOTH           present AND registered    — ready to use

"installed" here means "the tool-mcp's server.py file exists at the path
declared in mcps.yaml". "registered" means `claude mcp list` reports this
name. Registration is delegated to the Claude Code CLI (`claude mcp add`).
"""
from __future__ import annotations

import os
import shutil
import subprocess
import sys
from dataclasses import dataclass, field
from enum import Enum
from pathlib import Path
from typing import Any

try:
    import termios
except ImportError:  # non-POSIX (Windows) — terminal corruption is moot there
    termios = None  # type: ignore[assignment]


def _run_claude(cmd: list[str], **kwargs: Any) -> subprocess.CompletedProcess[Any]:
    """Run `claude ...` with the parent's terminal protected.

    Claude Code's Node.js CLI puts the controlling TTY into raw mode via
    process.stdin.setRawMode(true). When it exits without restoring termios,
    the parent shell is left with ICRNL/ICANON/ECHO/ISIG cleared — typed
    characters stop echoing, Enter shows ^M, Ctrl-C shows ^C. We defend
    against this by:
      1. detaching the child from our stdin (so it can't read /dev/tty
         through inherited fd 0); and
      2. snapshotting and restoring our own termios across the call, so
         even if the child touched /dev/tty directly, the parent shell
         comes back exactly as the user left it.
    """
    saved = None
    if termios is not None and sys.stdin.isatty():
        try:
            saved = termios.tcgetattr(sys.stdin.fileno())
        except (termios.error, OSError):
            saved = None
    kwargs.setdefault("stdin", subprocess.DEVNULL)
    try:
        return subprocess.run(cmd, **kwargs)
    finally:
        if saved is not None:
            try:
                termios.tcsetattr(sys.stdin.fileno(), termios.TCSADRAIN, saved)
            except (termios.error, OSError):
                pass


class Status(str, Enum):
    NOT_INSTALLED = "NOT_INSTALLED"
    INSTALLED = "INSTALLED"
    REGISTERED = "REGISTERED"
    BOTH = "BOTH"


@dataclass
class MCP:
    name: str
    path: str
    runtime: str = "python"
    entry: str = "src/server.py"
    description: str = ""
    category: str = ""
    tool_count: int = 0
    tags: list[str] = field(default_factory=list)
    env_vars: dict[str, str] = field(default_factory=dict)
    dependencies: dict[str, Any] = field(default_factory=dict)

    def server_path(self, root: Path) -> Path:
        return root / self.path / self.entry

    def is_installed(self, root: Path) -> bool:
        return self.server_path(root).is_file()

    def is_registered(self) -> bool:
        claude = shutil.which("claude")
        if claude is None:
            return False
        try:
            out = _run_claude(
                [claude, "mcp", "list"],
                capture_output=True,
                text=True,
                timeout=10,
            )
        except (subprocess.TimeoutExpired, OSError):
            return False
        return self.name in out.stdout

    def status(self, root: Path) -> Status:
        installed = self.is_installed(root)
        registered = self.is_registered()
        if installed and registered:
            return Status.BOTH
        if installed:
            return Status.INSTALLED
        if registered:
            return Status.REGISTERED
        return Status.NOT_INSTALLED

    def resolved_env(self, root: Path) -> dict[str, str]:
        """Substitute ${REPO_ROOT} in env_vars."""
        out = {}
        for k, v in self.env_vars.items():
            out[k] = v.replace("${REPO_ROOT}", str(root)) if isinstance(v, str) else v
        return out

    def register(self, root: Path, python_bin: str | None = None) -> None:
        """Register this MCP with Claude Code via `claude mcp add`."""
        if not self.is_installed(root):
            raise FileNotFoundError(f"server.py missing: {self.server_path(root)}")
        claude = shutil.which("claude")
        if claude is None:
            raise RuntimeError("`claude` CLI not found on PATH — install Claude Code first")
        python_bin = python_bin or shutil.which("python") or "python"
        env_flags = []
        for k, v in self.resolved_env(root).items():
            env_flags += ["--env", f"{k}={v}"]
        cmd = [claude, "mcp", "add", self.name, *env_flags,
               "--", python_bin, str(self.server_path(root))]
        _run_claude(cmd, check=True)

    def unregister(self) -> None:
        claude = shutil.which("claude")
        if claude is None:
            return
        _run_claude([claude, "mcp", "remove", self.name], check=False)
