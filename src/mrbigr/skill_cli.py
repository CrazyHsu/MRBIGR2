"""`mrbigr-skill` CLI — deploy workflow skills for Agent Skills tools.

    mrbigr-skill list                    show skills available in skills/
    mrbigr-skill install <name>          copy one skill to the selected targets
    mrbigr-skill install-all             copy all active workflow skills
    mrbigr-skill uninstall <name>        remove it from the selected targets
"""
from __future__ import annotations

import argparse
import sys

from .skill.deployer import (
    ALL_TARGETS,
    active_skill_infos,
    available_skill_infos,
    deploy_all_skills,
    deploy_skill,
    has_legacy_flat_install,
    is_installed,
    resolve_targets,
    skills_source_dir,
    uninstall_skill,
)


def cmd_list(args: argparse.Namespace) -> int:
    try:
        targets = _targets(args)
    except ValueError as e:
        print(e, file=sys.stderr)
        return 1
    skills = available_skill_infos(include_deprecated=True)
    src = skills_source_dir()
    if not skills:
        print(f"No skills found in {src}")
        return 0
    print(f"Skills available in {src}:")
    print("Targets:")
    for target in targets:
        print(f"  {target.name}: {target.path}")
    active = {info.file_stem for info in active_skill_infos(include_deprecated=False)}
    for info in skills:
        installed = [target.name for target in targets if is_installed(info, target.path)]
        legacy = [target.name for target in targets if has_legacy_flat_install(info, target.path)]
        marker = "[installed]" if installed else "[          ]"
        flags = []
        if info.file_stem in active:
            flags.append("active")
        if info.deprecated:
            flags.append("deprecated")
        if installed:
            flags.append(f"installed={','.join(installed)}")
        if legacy:
            flags.append(f"legacy-flat={','.join(legacy)}")
        suffix = f" ({', '.join(flags)})" if flags else ""
        print(f"  {marker}  {info.skill_name:<24} source={info.source.name}{suffix}")
    return 0


def cmd_install(args: argparse.Namespace) -> int:
    try:
        targets = _targets(args)
        deployed = [
            (target.name, deploy_skill(
                args.name,
                target_dir=target.path,
                include_deprecated=args.include_deprecated,
            ))
            for target in targets
        ]
    except (FileNotFoundError, ValueError) as e:
        print(e, file=sys.stderr)
        return 1
    for target_name, path in deployed:
        print(f"deployed [{target_name}]: {path}")
    return 0


def cmd_install_all(args: argparse.Namespace) -> int:
    try:
        deployed = [
            (target.name, path)
            for target in _targets(args)
            for path in deploy_all_skills(
                target_dir=target.path,
                include_deprecated=args.include_deprecated,
            )
        ]
    except ValueError as e:
        print(e, file=sys.stderr)
        return 1
    if not deployed:
        print("No skills deployed")
        return 0
    for target_name, path in deployed:
        print(f"deployed [{target_name}]: {path}")
    return 0


def cmd_uninstall(args: argparse.Namespace) -> int:
    try:
        removed = [
            (target.name, uninstall_skill(args.name, target_dir=target.path))
            for target in _targets(args)
        ]
    except (FileNotFoundError, ValueError) as e:
        print(e, file=sys.stderr)
        return 1
    for target_name, was_removed in removed:
        print(f"{args.name} [{target_name}]: {'removed' if was_removed else 'not installed'}")
    return 0


def _targets(args: argparse.Namespace):
    return resolve_targets(targets=args.target, target_dirs=args.target_dir)


def _add_targets(parser: argparse.ArgumentParser) -> None:
    known = ", ".join((*ALL_TARGETS, "all"))
    parser.add_argument(
        "--target",
        action="append",
        default=None,
        help=f"built-in target alias ({known}); repeatable",
    )
    parser.add_argument(
        "--target-dir",
        action="append",
        default=None,
        help="custom Agent Skills directory; repeatable",
    )


def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(prog="mrbigr-skill", description="deploy MRBIGR2 workflow skills")
    sub = p.add_subparsers(dest="command", required=True)

    sp = sub.add_parser("list", help="list source skills and install status")
    _add_targets(sp)
    sp.set_defaults(func=cmd_list)

    sp = sub.add_parser("install", help="install one skill by source stem or Claude skill name")
    sp.add_argument("name")
    sp.add_argument("--include-deprecated", action="store_true", help="allow installing deprecated redirect stubs")
    _add_targets(sp)
    sp.set_defaults(func=cmd_install)

    sp = sub.add_parser("install-all", help="install all active MRBIGR2 workflow skills")
    sp.add_argument("--include-deprecated", action="store_true", help="also install deprecated redirect stubs")
    _add_targets(sp)
    sp.set_defaults(func=cmd_install_all)

    sp = sub.add_parser("uninstall", help="remove one installed skill directory")
    sp.add_argument("name")
    _add_targets(sp)
    sp.set_defaults(func=cmd_uninstall)
    return p


def main(argv: list[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    return args.func(args)


if __name__ == "__main__":
    raise SystemExit(main())
