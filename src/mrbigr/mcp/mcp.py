"""MCP dataclass — one entry per row of mcps.yaml.

Each MCP has a four-state lifecycle mirroring ProteinMCP:

    NOT_INSTALLED  server.py missing on disk AND not registered with the client
    INSTALLED      server.py present          but not yet registered
    REGISTERED     registered with the client but server.py missing (broken)
    BOTH           present AND registered    — ready to use

"installed" here means "the tool-mcp's server.py file exists at the path
declared in mcps.yaml". "registered" means the selected MCP client reports this
name. Direct registration is delegated to supported client CLIs.
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

DEFAULT_MCP_CLIENT = "claude"
DIRECT_MCP_CLIENTS = ("claude", "codex", "gemini")
EXPORT_CONFIG_FORMATS = ("mcpservers", "gemini", "opencode")


def _run_mcp_client(cmd: list[str], **kwargs: Any) -> subprocess.CompletedProcess[Any]:
    """Run an MCP client CLI with the parent's terminal protected.

    Claude Code's Node.js CLI puts the controlling TTY into raw mode via
    process.stdin.setRawMode(true). When it exits without restoring termios,
    the parent shell is left with ICRNL/ICANON/ECHO/ISIG cleared — typed
    characters stop echoing, Enter shows ^M, Ctrl-C shows ^C.

    We defend by snapshotting and restoring our own termios across the
    call. subprocess.run blocks until the child exits, so the restore
    always completes before control returns to the user's shell — any
    intermediate raw-mode state is invisible to the user.

    We deliberately do NOT redirect the child's stdin to /dev/null: that
    breaks some MCP client commands, which appear to probe fd 0 at
    startup and behave unpredictably (empty output, silent exit) when it
    is not the terminal. The termios save/restore is sufficient on its
    own.
    """
    saved = None
    if termios is not None and sys.stdin.isatty():
        try:
            saved = termios.tcgetattr(sys.stdin.fileno())
        except (termios.error, OSError):
            saved = None
    try:
        return subprocess.run(cmd, **kwargs)
    finally:
        if saved is not None:
            try:
                termios.tcsetattr(sys.stdin.fileno(), termios.TCSADRAIN, saved)
            except (termios.error, OSError):
                pass


def normalize_mcp_client(client: str | None = None) -> str:
    normalized = (client or DEFAULT_MCP_CLIENT).lower()
    if normalized not in DIRECT_MCP_CLIENTS:
        known = ", ".join(DIRECT_MCP_CLIENTS)
        raise ValueError(f"unsupported direct MCP client: {client}. Known direct clients: {known}")
    return normalized


def mcp_client_list_output(client: str | None = None) -> str:
    client = normalize_mcp_client(client)
    cli = shutil.which(client)
    if cli is None:
        return ""
    try:
        out = _run_mcp_client(
            [cli, "mcp", "list"],
            capture_output=True,
            text=True,
            timeout=10,
        )
    except (subprocess.TimeoutExpired, OSError):
        return ""
    return f"{out.stdout}\n{out.stderr or ''}"


def status_from_flags(installed: bool, registered: bool) -> "Status":
    if installed and registered:
        return Status.BOTH
    if installed:
        return Status.INSTALLED
    if registered:
        return Status.REGISTERED
    return Status.NOT_INSTALLED


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

    def is_registered(self, client: str | None = None, list_output: str | None = None) -> bool:
        output = mcp_client_list_output(client) if list_output is None else list_output
        return self.name in output

    def status(self, root: Path, client: str | None = None, list_output: str | None = None) -> Status:
        installed = self.is_installed(root)
        registered = self.is_registered(client=client, list_output=list_output)
        return status_from_flags(installed, registered)

    def resolved_env(self, root: Path) -> dict[str, str]:
        """Substitute ${REPO_ROOT} in env_vars and inject runtime defaults.

        FASTMCP_CHECK_FOR_UPDATES=off disables FastMCP's on-startup call to
        PyPI. The check is pointless for an MCP server (no interactive user
        to see the notice) and it hard-fails the process on networks where
        httpx cannot construct a proxy transport — e.g. shells with
        `all_proxy=socks5://...` but without `httpx[socks]` installed, or
        offline environments. Can be overridden by setting it explicitly in
        mcps.yaml's env_vars.
        """
        defaults: dict[str, str] = {"FASTMCP_CHECK_FOR_UPDATES": "off"}
        out = dict(defaults)
        for k, v in self.env_vars.items():
            out[k] = v.replace("${REPO_ROOT}", str(root)) if isinstance(v, str) else v
        return out

    def launch_command(self, root: Path, python_bin: str | None = None) -> list[str]:
        python_bin = python_bin or shutil.which("python") or "python"
        return [python_bin, str(self.server_path(root))]

    def mcpservers_config(self, root: Path, python_bin: str | None = None) -> dict[str, Any]:
        command = self.launch_command(root, python_bin=python_bin)
        return {
            "command": command[0],
            "args": command[1:],
            "env": self.resolved_env(root),
        }

    def opencode_config(self, root: Path, python_bin: str | None = None) -> dict[str, Any]:
        return {
            "type": "local",
            "command": self.launch_command(root, python_bin=python_bin),
            "environment": self.resolved_env(root),
            "enabled": True,
        }

    def register(self, root: Path, python_bin: str | None = None, client: str | None = None) -> None:
        """Register this MCP with a supported MCP client.

        For clients that support scopes, user scope makes the MCP visible from
        any working directory. Local/project scope ties registration to the cwd
        where the add command ran, which would break sessions launched from
        anywhere else.

        Idempotent: `mcp remove` is always attempted first with
        `check=False`, so the "not registered" error on a fresh install
        is swallowed and re-installs cleanly overwrite stale paths /
        env vars (e.g. after switching conda envs or moving the repo).
        We do NOT gate this on `is_registered()` — under subprocess
        wrappers that check has proven unreliable in practice, and
        always-remove is simpler and correct.
        """
        if not self.is_installed(root):
            raise FileNotFoundError(f"server.py missing: {self.server_path(root)}")
        client = normalize_mcp_client(client)
        cli = shutil.which(client)
        if cli is None:
            raise RuntimeError(f"`{client}` CLI not found on PATH — install {client} first or export config instead")
        if client == "gemini":
            gemini_home = Path(os.environ.get("HOME", str(Path.home()))) / ".gemini"
            gemini_home.mkdir(parents=True, exist_ok=True)
            projects_json = gemini_home / "projects.json"
            if not projects_json.exists():
                projects_json.write_text('{"projects":{}}\n', encoding="utf-8")
        python_bin = python_bin or shutil.which("python") or "python"

        remove_cmd = [cli, "mcp", "remove"]
        if client in {"claude", "gemini"}:
            remove_cmd += ["--scope", "user"]
        remove_cmd.append(self.name)
        _run_mcp_client(remove_cmd, check=False)

        env_flags = []
        for k, v in self.resolved_env(root).items():
            env_flags += ["--env", f"{k}={v}"]

        add_cmd = [cli, "mcp", "add"]
        if client == "claude":
            add_cmd += ["--scope", "user", self.name, *env_flags, "--", python_bin, str(self.server_path(root))]
        elif client == "codex":
            add_cmd += [self.name, *env_flags, "--", python_bin, str(self.server_path(root))]
        elif client == "gemini":
            gemini_env_flags = []
            for k, v in self.resolved_env(root).items():
                gemini_env_flags += ["--env", f"{k}={v}"]
            add_cmd += ["--scope", "user", *gemini_env_flags, self.name, python_bin, str(self.server_path(root))]
        else:
            raise RuntimeError(f"unsupported MCP client: {client}")
        try:
            _run_mcp_client(add_cmd, check=True)
        except subprocess.CalledProcessError as e:
            raise RuntimeError(f"`{client} mcp add` failed with exit code {e.returncode}") from e

    def unregister(self, client: str | None = None) -> None:
        client = normalize_mcp_client(client)
        cli = shutil.which(client)
        if cli is None:
            return
        cmd = [cli, "mcp", "remove"]
        if client in {"claude", "gemini"}:
            cmd += ["--scope", "user"]
        cmd.append(self.name)
        _run_mcp_client(cmd, check=False)
