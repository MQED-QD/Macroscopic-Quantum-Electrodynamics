# mqed/utils/file_utils.py

import os
from pathlib import Path
from omegaconf import DictConfig
from loguru import logger

def _resolve_path(p: str) -> Path:
    """Expand env vars and ~, return absolute Path."""
    return Path(os.path.expandvars(os.path.expanduser(p))).resolve()


def _find_newest(pattern: str) -> Path | None:
    """Return newest file matching a glob pattern, or None."""
    logger.debug(f"Finding newest file matching pattern: {pattern}")
    from glob import glob
    hits = sorted(glob(pattern), key=lambda s: Path(s).stat().st_mtime, reverse=True)
    return Path(hits[0]).resolve() if hits else None

def _resolve_input_path(curve_cfg) -> Path:
    """
    Either use an absolute 'path', or if 'use_latest_glob' is set,
    choose newest file matching the glob (relative to MQED_ROOT/PWD if present).
    """
    if getattr(curve_cfg, "path", None):
        return _resolve_path(curve_cfg.path)
    if getattr(curve_cfg, "use_latest_glob", None):
        base = os.environ.get("MQED_ROOT", os.environ.get("PWD", "."))
        newest = _find_newest(os.path.join(base, curve_cfg.use_latest_glob))
        if newest is None:
            raise FileNotFoundError(f"No files for pattern: {curve_cfg.use_latest_glob}")
        return newest
    raise ValueError("Each curve needs either 'path' or 'use_latest_glob'.")
