# mqed/utils/file_utils.py

from pathlib import Path
from omegaconf import DictConfig
from loguru import logger

def find_latest_file(directory: str, pattern: str = "result_*.hdf5"):
    """
    Finds the latest file in a directory that matches a given pattern.
    
    Args:
        directory (str): The directory to search.
        pattern (str): The filename pattern (e.g., 'result_*.hdf5').

    Returns:
        Path: The most recently modified matching file.

    Raises:
        FileNotFoundError: If no matching files are found.
    """
    if isinstance(directory, DictConfig):
        directory = directory.get("dir", str(directory))
    dir_path = Path(str(directory))
    dir_path = Path(directory)
    files = list(dir_path.glob(pattern))
    if not files:
        raise FileNotFoundError(f"No files matching {pattern} found in {dir_path}")

    latest_file = max(files, key=lambda f: f.stat().st_mtime)
    logger.info(f"Found latest file: {latest_file}")
    return latest_file
