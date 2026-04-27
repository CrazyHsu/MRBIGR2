"""Deploy workflow skills from ``skills/*.md`` into agent skill directories.

Claude Code, Codex, Gemini CLI, OpenCode, and other Agent Skills-compatible tools discover
personal skills from directories shaped like
``<skills-root>/<skill-name>/SKILL.md``.  The source files in this repository
are intentionally kept as flat markdown files for easy review, so this module
translates that source layout into each agent's runtime layout.
"""
from __future__ import annotations

import os
import shutil
from dataclasses import dataclass
from pathlib import Path

CLAUDE_SKILLS_DIR = Path.home() / ".claude" / "skills"
GEMINI_SKILLS_DIR = Path.home() / ".gemini" / "skills"
OPENCODE_SKILLS_DIR = Path.home() / ".config" / "opencode" / "skills"
AGENTS_SKILLS_DIR = Path.home() / ".agents" / "skills"
ACTIVE_SKILL_STEMS = (
    "gwas_pipeline",
    "qtl_to_target",
    "causal_network",
    "functional_enrichment",
    "reproduce_v1_case",
)
DEPRECATED_SKILL_STEMS = {"mr_analysis", "qtl_mapping"}
DEFAULT_TARGET = "claude"
ALL_TARGETS = ("claude", "codex", "gemini", "opencode", "agents")
DEFAULT_ALL_TARGETS = ("claude", "codex", "gemini", "opencode")


@dataclass(frozen=True)
class SkillInstallTarget:
    name: str
    path: Path


@dataclass(frozen=True)
class SkillInfo:
    source: Path
    file_stem: str
    skill_name: str
    deprecated: bool = False

    @property
    def target_dir_name(self) -> str:
        return self.skill_name

    def target_file(self, target_dir: Path | str = CLAUDE_SKILLS_DIR) -> Path:
        return _as_path(target_dir) / self.target_dir_name / "SKILL.md"

    def legacy_file(self, target_dir: Path | str = CLAUDE_SKILLS_DIR) -> Path:
        return _as_path(target_dir) / self.source.name

    def matches(self, name: str) -> bool:
        normalized = name.removesuffix(".md")
        return normalized in {
            self.file_stem,
            self.skill_name,
            self.file_stem.replace("_", "-"),
            self.skill_name.replace("-", "_"),
        }


def _as_path(path: Path | str) -> Path:
    return Path(path).expanduser()


def codex_skills_dir() -> Path:
    codex_home = os.environ.get("CODEX_HOME")
    root = Path(codex_home).expanduser() if codex_home else Path.home() / ".codex"
    return root / "skills"


def built_in_target_paths() -> dict[str, Path]:
    return {
        "claude": CLAUDE_SKILLS_DIR,
        "codex": codex_skills_dir(),
        "gemini": GEMINI_SKILLS_DIR,
        "opencode": OPENCODE_SKILLS_DIR,
        "agents": AGENTS_SKILLS_DIR,
    }


def resolve_targets(
    targets: list[str] | tuple[str, ...] | None = None,
    target_dirs: list[Path | str] | tuple[Path | str, ...] | None = None,
) -> list[SkillInstallTarget]:
    """Resolve built-in target aliases and custom directories.

    No explicit target preserves the historical behavior: install into the
    Claude target. The ``all`` alias installs into concrete client targets only;
    the shared ``agents`` directory is explicit opt-in to avoid duplicate skill
    discovery in clients that scan both their own directory and ``~/.agents``.
    """
    target_names = list(targets or [])
    custom_dirs = list(target_dirs or [])
    if not target_names and not custom_dirs:
        target_names = [DEFAULT_TARGET]

    aliases = built_in_target_paths()
    resolved: list[SkillInstallTarget] = []
    for target in target_names:
        name = target.lower()
        if name == "all":
            resolved.extend(SkillInstallTarget(alias, aliases[alias]) for alias in DEFAULT_ALL_TARGETS)
            continue
        if name not in aliases:
            known = ", ".join((*ALL_TARGETS, "all"))
            raise ValueError(f"unknown skill target: {target}. Known targets: {known}; use --target-dir for custom paths")
        resolved.append(SkillInstallTarget(name, aliases[name]))

    for path in custom_dirs:
        expanded = _as_path(path)
        resolved.append(SkillInstallTarget(str(expanded), expanded))

    return _dedupe_targets(resolved)


def _dedupe_targets(targets: list[SkillInstallTarget]) -> list[SkillInstallTarget]:
    seen: set[str] = set()
    deduped: list[SkillInstallTarget] = []
    for target in targets:
        key = str(target.path)
        if key in seen:
            continue
        seen.add(key)
        deduped.append(target)
    return deduped


def repo_root() -> Path:
    """Return the MRBIGR2 repository root without importing the heavy core."""
    env_root = os.environ.get("MRBIGR_ROOT")
    if env_root:
        return Path(env_root).expanduser().resolve()

    here = Path(__file__).resolve()
    for candidate in (here, *here.parents):
        if (candidate / "skills").is_dir() and (candidate / "pyproject.toml").is_file():
            return candidate
    return here.parents[3]


def skills_source_dir() -> Path:
    return repo_root() / "skills"


def available_skills() -> list[Path]:
    """Return source skill markdown files for backward-compatible callers."""
    return [info.source for info in available_skill_infos(include_deprecated=True)]


def available_skill_infos(include_deprecated: bool = True) -> list[SkillInfo]:
    src = skills_source_dir()
    if not src.is_dir():
        return []
    infos = [_skill_info(path) for path in sorted(src.glob("*.md"))]
    if include_deprecated:
        return infos
    return [info for info in infos if not info.deprecated]


def active_skill_infos(include_deprecated: bool = False) -> list[SkillInfo]:
    infos = available_skill_infos(include_deprecated=True)
    if include_deprecated:
        return infos
    active = set(ACTIVE_SKILL_STEMS)
    return [info for info in infos if info.file_stem in active and not info.deprecated]


def resolve_skill(name: str, include_deprecated: bool = False) -> SkillInfo:
    """Resolve a source file stem or runtime skill name to a skill manifest."""
    matches = [info for info in available_skill_infos(include_deprecated=True) if info.matches(name)]
    if not matches:
        raise FileNotFoundError(f"skill not found: {name}")
    info = matches[0]
    if info.deprecated and not include_deprecated:
        raise ValueError(f"skill is deprecated: {info.skill_name}; pass --include-deprecated to install it")
    return info


def deploy_skill(
    name: str,
    target_dir: Path | str = CLAUDE_SKILLS_DIR,
    include_deprecated: bool = False,
) -> Path:
    """Copy one source skill into an Agent Skills-compatible layout.

    Returns the path to the deployed skill. Raises FileNotFoundError when
    the source skill does not exist.
    """
    info = resolve_skill(name, include_deprecated=include_deprecated)
    target = info.target_file(target_dir)
    target.parent.mkdir(parents=True, exist_ok=True)
    shutil.copy2(info.source, target)
    return target


def deploy_all_skills(
    target_dir: Path | str = CLAUDE_SKILLS_DIR,
    include_deprecated: bool = False,
) -> list[Path]:
    """Deploy the default active skills, optionally including deprecated stubs."""
    return [
        deploy_skill(info.skill_name, target_dir=target_dir, include_deprecated=True)
        for info in active_skill_infos(include_deprecated=include_deprecated)
    ]


def is_installed(info: SkillInfo, target_dir: Path | str = CLAUDE_SKILLS_DIR) -> bool:
    return info.target_file(target_dir).is_file()


def has_legacy_flat_install(info: SkillInfo, target_dir: Path | str = CLAUDE_SKILLS_DIR) -> bool:
    return info.legacy_file(target_dir).is_file()


def uninstall_skill(name: str, target_dir: Path | str = CLAUDE_SKILLS_DIR) -> bool:
    info = resolve_skill(name, include_deprecated=True)
    skill_dir = _as_path(target_dir) / info.target_dir_name
    if skill_dir.is_dir():
        shutil.rmtree(skill_dir)
        return True
    return False


def uninstall_all_skills(
    target_dir: Path | str = CLAUDE_SKILLS_DIR,
    include_deprecated: bool = False,
) -> list[tuple[SkillInfo, bool]]:
    """Remove all MRBIGR2-managed active skills from one target directory."""
    removed: list[tuple[SkillInfo, bool]] = []
    for info in active_skill_infos(include_deprecated=include_deprecated):
        skill_dir = _as_path(target_dir) / info.target_dir_name
        if skill_dir.is_dir():
            shutil.rmtree(skill_dir)
            removed.append((info, True))
        else:
            removed.append((info, False))
    return removed


def _skill_info(path: Path) -> SkillInfo:
    frontmatter = _read_frontmatter(path)
    skill_name = frontmatter.get("name") or path.stem.replace("_", "-")
    deprecated = path.stem in DEPRECATED_SKILL_STEMS or skill_name.endswith("-deprecated")
    return SkillInfo(
        source=path,
        file_stem=path.stem,
        skill_name=skill_name,
        deprecated=deprecated,
    )


def _read_frontmatter(path: Path) -> dict[str, str]:
    lines = path.read_text(encoding="utf-8").splitlines()
    if not lines or lines[0].strip() != "---":
        return {}

    data: dict[str, str] = {}
    for line in lines[1:]:
        stripped = line.strip()
        if stripped == "---":
            break
        if ":" not in line:
            continue
        key, value = line.split(":", 1)
        key = key.strip()
        if not key:
            continue
        data[key] = value.strip().strip("'\"")
    return data
