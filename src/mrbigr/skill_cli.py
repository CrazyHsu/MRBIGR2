"""`mrbigr-skill` CLI — deploy workflow skills for Claude Code.

    mrbigr-skill list                    show skills available in skills/
    mrbigr-skill install <name>          copy skills/<name>.md to ~/.claude/skills/
    mrbigr-skill uninstall <name>        remove it from ~/.claude/skills/
"""
from __future__ import annotations

import argparse
import sys

from .skill.deployer import (
    CLAUDE_SKILLS_DIR,
    available_skills,
    deploy_skill,
    skills_source_dir,
    uninstall_skill,
)


def cmd_list(args: argparse.Namespace) -> int:
    skills = available_skills()
    src = skills_source_dir()
    if not skills:
        print(f"No skills found in {src}")
        return 0
    print(f"Skills available in {src}:")
    installed = {p.stem for p in CLAUDE_SKILLS_DIR.glob("*.md")} if CLAUDE_SKILLS_DIR.is_dir() else set()
    for path in skills:
        marker = "[installed]" if path.stem in installed else "[          ]"
        print(f"  {marker}  {path.stem}")
    return 0


def cmd_install(args: argparse.Namespace) -> int:
    try:
        target = deploy_skill(args.name)
    except FileNotFoundError as e:
        print(e, file=sys.stderr)
        return 1
    print(f"deployed: {target}")
    return 0


def cmd_uninstall(args: argparse.Namespace) -> int:
    removed = uninstall_skill(args.name)
    print(f"{args.name}: {'removed' if removed else 'not installed'}")
    return 0


def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(prog="mrbigr-skill", description="deploy MRBIGR2 workflow skills")
    sub = p.add_subparsers(dest="command", required=True)
    sub.add_parser("list").set_defaults(func=cmd_list)
    sp = sub.add_parser("install"); sp.add_argument("name"); sp.set_defaults(func=cmd_install)
    sp = sub.add_parser("uninstall"); sp.add_argument("name"); sp.set_defaults(func=cmd_uninstall)
    return p


def main(argv: list[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    return args.func(args)


if __name__ == "__main__":
    raise SystemExit(main())
