"""`mrbigr` CLI — meta-orchestrator for MRBIGR2 MCPs.

Subcommands mirror `pmcp` from ProteinMCP:

    mrbigr list                     show all registered MCPs and their status
    mrbigr status <name>            show a single MCP's lifecycle state
    mrbigr install <name>           register the MCP with a supported MCP client
    mrbigr uninstall <name>         unregister from a supported MCP client
    mrbigr install-all              register every MCP whose server.py exists
    mrbigr export-config [name|all] emit MCP client configuration JSON
"""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

from .mcp.manager import MCPManager
from .mcp.mcp import (
    DEFAULT_MCP_CLIENT,
    DIRECT_MCP_CLIENTS,
    EXPORT_CONFIG_FORMATS,
    Status,
    mcp_client_list_output,
)


def _format_row(name: str, status: Status, tool_count: int, description: str) -> str:
    status_col = f"{status.value:<14}"
    tool_col = f"({tool_count:>2} tools)" if tool_count else "         "
    return f"  {name:<12}  {status_col}  {tool_col}  {description}"


def cmd_list(args: argparse.Namespace) -> int:
    mgr = MCPManager()
    try:
        clients = _clients(args)
    except ValueError as e:
        print(e, file=sys.stderr)
        return 1
    mcps = mgr.list()
    if not mcps:
        print(f"No MCPs found in {mgr.registry_path}")
        return 0
    print(f"MRBIGR2 registry: {mgr.registry_path}  ({len(mcps)} MCPs)")
    for client in clients:
        list_output = mcp_client_list_output(client)
        print()
        print(f"Client: {client}")
        print(_format_row("NAME", Status.NOT_INSTALLED, 0, "DESCRIPTION").replace(
            Status.NOT_INSTALLED.value.ljust(14), "STATUS".ljust(14)
        ).replace("(0 tools)", "         "))
        print("  " + "-" * 76)
        for m in mcps:
            print(_format_row(m.name, m.status(mgr.root, client=client, list_output=list_output), m.tool_count, m.description))
    return 0


def cmd_status(args: argparse.Namespace) -> int:
    mgr = MCPManager()
    try:
        clients = _clients(args)
    except ValueError as e:
        print(e, file=sys.stderr)
        return 1
    try:
        mcp = mgr.get(args.name)
    except KeyError as e:
        print(e, file=sys.stderr)
        return 1
    for client in clients:
        list_output = mcp_client_list_output(client)
        status = mcp.status(mgr.root, client=client, list_output=list_output)
        print(f"{mcp.name} [{client}]: {status.value}")
        print(f"  path:        {mgr.root / mcp.path}")
        print(f"  entry:       {mcp.entry}")
        print(f"  tool_count:  {mcp.tool_count}")
        print(f"  category:    {mcp.category}")
        print(f"  description: {mcp.description}")
        print(f"  server.py:   {'present' if mcp.is_installed(mgr.root) else 'MISSING'}")
        print(f"  registered:  {mcp.is_registered(client=client, list_output=list_output)}")
    return 0


def cmd_install(args: argparse.Namespace) -> int:
    mgr = MCPManager()
    try:
        clients = _clients(args)
    except ValueError as e:
        print(e, file=sys.stderr)
        return 1
    rc = 0
    for client in clients:
        try:
            new_status = mgr.install(args.name, python_bin=args.python, client=client)
        except KeyError as e:
            print(e, file=sys.stderr)
            return 1
        except FileNotFoundError as e:
            print(f"Install failed [{client}]: {e}", file=sys.stderr)
            rc = 2
            continue
        except RuntimeError as e:
            print(f"Install failed [{client}]: {e}", file=sys.stderr)
            rc = 3
            continue
        print(f"{args.name} [{client}]: {new_status.value}")
    return rc


def cmd_uninstall(args: argparse.Namespace) -> int:
    mgr = MCPManager()
    try:
        clients = _clients(args)
    except ValueError as e:
        print(e, file=sys.stderr)
        return 1
    rc = 0
    for client in clients:
        try:
            new_status = mgr.uninstall(args.name, client=client)
        except KeyError as e:
            print(e, file=sys.stderr)
            return 1
        except RuntimeError as e:
            print(f"Uninstall failed [{client}]: {e}", file=sys.stderr)
            rc = 3
            continue
        print(f"{args.name} [{client}]: {new_status.value}")
    return rc


def cmd_install_all(args: argparse.Namespace) -> int:
    mgr = MCPManager()
    try:
        clients = _clients(args)
    except ValueError as e:
        print(e, file=sys.stderr)
        return 1
    rc = 0
    for client in clients:
        for mcp in mgr.list():
            if not mcp.is_installed(mgr.root):
                print(f"skip {mcp.name} [{client}]: server.py not yet created")
                continue
            try:
                new_status = mgr.install(mcp.name, python_bin=args.python, client=client)
                print(f"{mcp.name} [{client}]: {new_status.value}")
            except Exception as e:  # noqa: BLE001 — surface any failure, continue others
                print(f"{mcp.name} [{client}]: FAILED ({e})", file=sys.stderr)
                rc = 4
    return rc


def cmd_export_config(args: argparse.Namespace) -> int:
    mgr = MCPManager()
    try:
        selected = mgr.list() if args.name == "all" else [mgr.get(args.name)]
    except KeyError as e:
        print(e, file=sys.stderr)
        return 1

    if args.format in {"mcpservers", "gemini"}:
        data = {
            "mcpServers": {
                mcp.name: mcp.mcpservers_config(mgr.root, python_bin=args.python)
                for mcp in selected
            }
        }
    elif args.format == "opencode":
        data = {
            "$schema": "https://opencode.ai/config.json",
            "mcp": {
                mcp.name: mcp.opencode_config(mgr.root, python_bin=args.python)
                for mcp in selected
            },
        }
    else:
        print(f"Unknown export format: {args.format}", file=sys.stderr)
        return 1

    text = json.dumps(data, indent=2, ensure_ascii=False) + "\n"
    if args.output:
        output = Path(args.output).expanduser()
        output.parent.mkdir(parents=True, exist_ok=True)
        output.write_text(text, encoding="utf-8")
        print(f"wrote {args.format} config: {output}")
    else:
        print(text, end="")
    return 0


def _clients(args: argparse.Namespace) -> list[str]:
    requested = args.client or [DEFAULT_MCP_CLIENT]
    clients: list[str] = []
    for value in requested:
        name = value.lower()
        if name == "all":
            clients.extend(DIRECT_MCP_CLIENTS)
        elif name in DIRECT_MCP_CLIENTS:
            clients.append(name)
        elif name == "opencode":
            raise ValueError(
                "OpenCode's `opencode mcp add` is interactive; use "
                "`mrbigr export-config --format opencode` and merge the generated mcp block."
            )
        else:
            known = ", ".join((*DIRECT_MCP_CLIENTS, "all"))
            raise ValueError(f"unknown MCP client: {value}. Known direct clients: {known}")
    return _dedupe(clients)


def _dedupe(values: list[str]) -> list[str]:
    seen: set[str] = set()
    out: list[str] = []
    for value in values:
        if value in seen:
            continue
        seen.add(value)
        out.append(value)
    return out


def _add_client_arg(parser: argparse.ArgumentParser) -> None:
    known = ", ".join((*DIRECT_MCP_CLIENTS, "all"))
    parser.add_argument(
        "--client",
        action="append",
        default=None,
        help=f"MCP client for direct registration/status ({known}); repeatable",
    )


_ORCHESTRATOR_COMMANDS = {"list", "status", "install", "uninstall", "install-all", "export-config"}


def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(prog="mrbigr", description="MRBIGR2 MCP orchestrator")
    sub = p.add_subparsers(dest="command", required=True)

    sp = sub.add_parser("list", help="list all MCPs and their status")
    _add_client_arg(sp)
    sp.set_defaults(func=cmd_list)

    sp = sub.add_parser("status", help="show status of one MCP")
    sp.add_argument("name")
    _add_client_arg(sp)
    sp.set_defaults(func=cmd_status)

    sp = sub.add_parser("install", help="register an MCP with a supported MCP client")
    sp.add_argument("name")
    sp.add_argument("--python", default=None, help="python binary to launch the server with")
    _add_client_arg(sp)
    sp.set_defaults(func=cmd_install)

    sp = sub.add_parser("uninstall", help="unregister an MCP")
    sp.add_argument("name")
    _add_client_arg(sp)
    sp.set_defaults(func=cmd_uninstall)

    sp = sub.add_parser("install-all", help="install every MCP whose server.py exists")
    sp.add_argument("--python", default=None)
    _add_client_arg(sp)
    sp.set_defaults(func=cmd_install_all)

    sp = sub.add_parser("export-config", help="export MCP client configuration JSON")
    sp.add_argument("name", nargs="?", default="all", help="MCP name or 'all' (default: all)")
    sp.add_argument("--format", choices=EXPORT_CONFIG_FORMATS, default="mcpservers")
    sp.add_argument("--python", default=None, help="python binary to launch the server with")
    sp.add_argument("-o", "--output", default=None, help="write JSON to this path instead of stdout")
    sp.set_defaults(func=cmd_export_config)

    return p


def _delegate_to_legacy_cli(argv: list[str]) -> int:
    """Forward non-orchestrator commands to the legacy file-based CLI."""
    cli_py = Path(__file__).resolve().parent.parent / "mrbigr_cli.py"
    if not cli_py.is_file():
        print(f"Unknown command: {argv[0]}", file=sys.stderr)
        print(f"Valid orchestrator commands: {', '.join(sorted(_ORCHESTRATOR_COMMANDS))}")
        print("Legacy CLI not found at expected location.")
        return 1
    import subprocess
    return subprocess.call([sys.executable, str(cli_py)] + argv)


def main(argv: list[str] | None = None) -> int:
    if argv is None:
        argv = sys.argv[1:]
    if not argv or argv[0] in ("-h", "--help"):
        build_parser().parse_args(argv)
        return 0
    if argv[0] not in _ORCHESTRATOR_COMMANDS:
        return _delegate_to_legacy_cli(argv)
    parser = build_parser()
    args = parser.parse_args(argv)
    return args.func(args)


if __name__ == "__main__":
    raise SystemExit(main())
