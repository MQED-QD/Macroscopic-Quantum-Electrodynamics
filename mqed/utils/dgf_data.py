# mqed/utils/dgf_data.py
from __future__ import annotations
import h5py
import numpy as np
from typing import Dict, Optional
from loguru import logger

def load_gf_h5(h5_path: str) -> Dict[str, np.ndarray]:
    """Load dyadic Green's function arrays from your HDF5 file.
    Args:
        h5_path (str): Path to the HDF5 file to read.
    Returns:
        Dict[str, np.ndarray]: A dictionary containing:
            - "G_total": Total Green's function array of shape (M, N, 3, 3).
            - "G_vac": Vacuum Green's function array of shape (M, N, 3, 3).
            - "energy_eV": Energy array of shape (M,).
            - "Rx_nm": Donor-acceptor distance array in nm of shape (N,).
            - "zD": Fixed z-position of the donor in meters (float).
            - "zA": Fixed z-position of the acceptor in meters (float).
    """
    try:
        with h5py.File(h5_path, "r") as f:
            Gtot = f["green_function_total"][:]      # (M, N, 3, 3)
            Gvac = f["green_function_vacuum"][:]     # (M, N, 3, 3)
            E    = f["energy_eV"][:].astype(float)   # (M,)
            Rxnm = f["Rx_nm"][:].astype(float)       # (N,)
            pos  = f["position_fixed"]
            zD   = float(pos.attrs["zD_meters"])
            zA   = float(pos.attrs["zA_meters"])
            logger.success(f"Loaded Green's function data from {h5_path}")
    except FileNotFoundError as e:
        logger.exception(f"HDF5 file not found: {h5_path}")
        raise
    except KeyError as e:
        logger.exception(f"Missing dataset in HDF5 file: {e}")
        raise
    except Exception as e:
        logger.exception(f"Error loading Green's function data: {e}")
        raise
    return {"G_total": Gtot, "G_vac": Gvac, "energy_eV": E, "Rx_nm": Rxnm, "zD": zD, "zA": zA}

def save_gf_h5(h5_path: str, Gtot: np.ndarray, Gvac: np.ndarray, E: np.ndarray, Rxnm: np.ndarray, zD: float, zA: float) -> None:
    """Save dyadic Green's function arrays to an HDF5 file.
    Args:
        h5_path (str): Path to the HDF5 file to create.
        Gtot (np.ndarray): Total Green's function array of shape (M, N, 3, 3).
        Gvac (np.ndarray): Vacuum Green's function array of shape (M, N, 3, 3).
        E (np.ndarray): Energy array of shape (M,).
        Rxnm (np.ndarray): Donor-acceptor distance array in nm of shape (N,).
        zD (float): Fixed z-position of the donor in meters.
        zA (float): Fixed z-position of the acceptor in meters.
    """
    with h5py.File(h5_path, "w") as f:
        f.create_dataset("green_function_total", data=Gtot)
        f.create_dataset("green_function_vacuum", data=Gvac)
        f.create_dataset("energy_eV", data=E)
        f.create_dataset("Rx_nm", data=Rxnm)
        pos = f.create_group("position_fixed")
        pos.attrs["zD_meters"] = zD
        pos.attrs["zA_meters"] = zA
