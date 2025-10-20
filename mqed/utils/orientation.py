# mqed/utils/orientation.py

import numpy as np
from loguru import logger

def resolve_angle_deg(v):
    MAGIC_DEG = float(np.degrees(np.arccos(1/np.sqrt(3))))  # 54.735610317245346
    if isinstance(v, (int, float)):     # already numeric
        return float(v)
    s = str(v).strip().lower()
    if s in {"magic", "ma", "magic_angle"}:
        return MAGIC_DEG
    if s in {"acos(1/sqrt(3))", "arccos(1/sqrt(3))"}:
        return MAGIC_DEG
    raise ValueError(f"Unrecognized angle spec: {v!r}")


def spherical_to_cartesian_dipole(theta_deg, phi_deg):
    """
    Converts spherical coordinates (theta, azimuthal_angle) to a Cartesian unit vector.

    Args:
        theta_deg (float or array-like): Polar angle in degrees (0 to 180, from +z axis).
        phi_deg (float or array-like): Azimuthal angle in degrees (0 to 360, from +x axis towards +y).

    Returns:
        np.ndarray: A 3-element numpy array representing the Cartesian unit vector [x, y, z].
    """
    theta = np.asanyarray(theta_deg, dtype=float)
    phi = np.asanyarray(phi_deg, dtype= float)

    # Broadcast to a common shape (e.g., (N,))
    theta_b, phi_b = np.broadcast_arrays(theta, phi)

    theta_rad = np.deg2rad(theta_b)
    phi_rad = np.deg2rad(phi_b)

    x = np.sin(theta_rad) * np.cos(phi_rad)
    y = np.sin(theta_rad) * np.sin(phi_rad)
    z = np.cos(theta_rad)

    vec = np.column_stack((x.ravel(), y.ravel(), z.ravel()))  # (N,3)
    # breakpoint()
    if np.isscalar(theta_deg) and np.isscalar(phi_deg):
        logger.debug(f"Theta={theta_deg} deg, Phi={phi_deg} deg -> Vector={vec[0]}")
        return vec[0]
    else:
        logger.debug(f"Generated {vec.shape[0]} vectors; first row={vec[0]}")
        return vec
