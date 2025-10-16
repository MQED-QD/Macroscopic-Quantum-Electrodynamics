from qutip import Qobj, qeye
import numpy as np


def msd_operator(dim: int, d_nm: float, Nmol: int, init_site_index: int) -> Qobj:
    """Mean-square displacement operator in the single-excitation manifold.


    Basis: |0> (ground), |1>,...,|N> (site excitations).
    Positions are 0 for ground and j*d for site j. MSD = (X - x0 I)^2.
    """
    positions = np.zeros(dim)
    positions[1:] = d_nm * np.arange(1, Nmol + 1, dtype=float)
    X = Qobj(np.diag(positions), dims=[[dim], [dim]])
    x0 = positions[init_site_index]
    return (X - x0 * qeye(dim)) ** 2