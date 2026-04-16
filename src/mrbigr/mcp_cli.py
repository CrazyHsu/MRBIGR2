"""`mrbigr` CLI — meta-orchestrator for MRBIGR2 MCPs.

Subcommands mirror `pmcp` from ProteinMCP:

    mrbigr list                     show all registered MCPs and their status
    mrbigr status <name>            show a single MCP's lifecycle state
    mrbigr install <name>           register the MCP with Claude Code
    mrbigr uninstall <name>         unregister from Claude Code
    mrbigr install-all              register every MCP whose server.py exists
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

from .mcp.manager import MCPManager
from .mcp.mcp import Status


def _format_row(name: str, status: Status, tool_count: int, description: str) -> str:
    status_col = f"{status.value:<14}"
    tool_col = f"({tool_count:>2} tools)" if tool_count else "         "
    return f"  {name:<12}  {status_col}  {tool_col}  {description}"


def cmd_list(args: argparse.Namespace) -> int:
    mgr = MCPManager()
    mcps = mgr.list()
    if not mcps:
        print(f"No MCPs found in {mgr.registry_path}")
        return 0
    print(f"MRBIGR2 registry: {mgr.registry_path}  ({len(mcps)} MCPs)")
    print()
    print(_format_row("NAME", Status.NOT_INSTALLED, 0, "DESCRIPTION").replace(
        Status.NOT_INSTALLED.value.ljust(14), "STATUS".ljust(14)
    ).replace("(0 tools)", "         "))
    print("  " + "-" * 76)
    for m in mcps:
        print(_format_row(m.name, m.status(mgr.root), m.tool_count, m.description))
    return 0


def cmd_status(args: argparse.Namespace) -> int:
    mgr = MCPManager()
    try:
        mcp = mgr.get(args.name)
    except KeyError as e:
        print(e, file=sys.stderr)
        return 1
    status = mcp.status(mgr.root)
    print(f"{mcp.name}: {status.value}")
    print(f"  path:        {mgr.root / mcp.path}")
    print(f"  entry:       {mcp.entry}")
    print(f"  tool_count:  {mcp.tool_count}")
    print(f"  category:    {mcp.category}")
    print(f"  description: {mcp.description}")
    print(f"  server.py:   {'present' if mcp.is_installed(mgr.root) else 'MISSING'}")
    print(f"  registered:  {mcp.is_registered()}")
    return 0


def cmd_install(args: argparse.Namespace) -> int:
    mgr = MCPManager()
    try:
        new_status = mgr.install(args.name, python_bin=args.python)
    except KeyError as e:
        print(e, file=sys.stderr)
        return 1
    except FileNotFoundError as e:
        print(f"Install failed: {e}", file=sys.stderr)
        return 2
    except RuntimeError as e:
        print(f"Install failed: {e}", file=sys.stderr)
        return 3
    print(f"{args.name}: {new_status.value}")
    return 0


def cmd_uninstall(args: argparse.Namespace) -> int:
    mgr = MCPManager()
    try:
        new_status = mgr.uninstall(args.name)
    except KeyError as e:
        print(e, file=sys.stderr)
        return 1
    print(f"{args.name}: {new_status.value}")
    return 0


def cmd_install_all(args: argparse.Namespace) -> int:
    mgr = MCPManager()
    rc = 0
    for mcp in mgr.list():
        if not mcp.is_installed(mgr.root):
            print(f"skip {mcp.name}: server.py not yet created")
            continue
        try:
            new_status = mgr.install(mcp.name, python_bin=args.python)
            print(f"{mcp.name}: {new_status.value}")
        except Exception as e:  # noqa: BLE001 — surface any failure, continue others
            print(f"{mcp.name}: FAILED ({e})", file=sys.stderr)
            rc = 4
    return rc


def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(prog="mrbigr", description="MRBIGR2 MCP orchestrator")
    sub = p.add_subparsers(dest="command", required=True)

    sub.add_parser("list", help="list all MCPs and their status").set_defaults(func=cmd_list)

    sp = sub.add_parser("status", help="show status of one MCP")
    sp.add_argument("name")
    sp.set_defaults(func=cmd_status)

    sp = sub.add_parser("install", help="register an MCP with Claude Code")
    sp.add_argument("name")
    sp.add_argument("--python", default=None, help="python binary to launch the server with")
    sp.set_defaults(func=cmd_install)

    sp = sub.add_parser("uninstall", help="unregister an MCP")
    sp.add_argument("name")
    sp.set_defaults(func=cmd_uninstall)

    sp = sub.add_parser("install-all", help="install every MCP whose server.py exists")
    sp.add_argument("--python", default=None)
    sp.set_defaults(func=cmd_install_all)

    return p


def main(argv: list[str] | None = None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)
    return args.func(args)


if __name__ == "__main__":
    raise SystemExit(main())
