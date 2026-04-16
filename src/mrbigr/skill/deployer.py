"""Deploy workflow skills from skills/*.md into Claude Code.

ProteinMCP's `pskill` does the same: copies curated markdown files into
~/.claude/skills/ where Claude Code picks them up as slash commands.
"""
from __future__ import annotations

import shutil
from pathlib import Path

from ..core.paths import repo_root

CLAUDE_SKILLS_DIR = Path.home() / ".claude" / "skills"


def skills_source_dir() -> Path:
    return repo_root() / "skills"


def available_skills() -> list[Path]:
    src = skills_source_dir()
    if not src.is_dir():
        return []
    return sorted(src.glob("*.md"))


def deploy_skill(name: str, target_dir: Path = CLAUDE_SKILLS_DIR) -> Path:
    """Copy skills/<name>.md into the Claude skills directory.

    Returns the path to the deployed skill. Raises FileNotFoundError when
    the source skill does not exist.
    """
    source = skills_source_dir() / f"{name}.md"
    if not source.is_file():
        raise FileNotFoundError(f"skill not found: {source}")
    target_dir.mkdir(parents=True, exist_ok=True)
    target = target_dir / source.name
    shutil.copy2(source, target)
    return target


def uninstall_skill(name: str, target_dir: Path = CLAUDE_SKILLS_DIR) -> bool:
    target = target_dir / f"{name}.md"
    if target.is_file():
        target.unlink()
        return True
    return False
