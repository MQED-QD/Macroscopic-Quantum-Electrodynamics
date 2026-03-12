"""Nearest-neighbour 1D chain dynamics with diagonal + off-diagonal disorder.

Physics
-------
Tight-binding Hamiltonian on N sites with nearest-neighbour coupling:
.. math::
    H_{nn'} = \epsilon_n \delta_{nn'} + J_n \delta_{n,n'+1} + J_n^* \delta_{n,n'-1}

Disorder is Gaussian:
.. math::
    \epsilon_n = \epsilon_0 + \delta\epsilon_n,    \delta\epsilon_n ~ N(0, \sigma_{\epsilon}^2)
    J_n   = J_0   + \delta J_n,    \delta J_n ~ N(0, \sigma_J^2)

MSD analytical solution for off-diagonal disorder for local excitation:
.. math::
    \langle x^2(t) \\rangle = 2 a^2 \frac{J_0^2 + \delta J_n^2}{\hbar^2} t^2

For plane wave excitation with wavevector k_parallel:
.. math::
    c_n(0) = e^{i k_{parallel} n a} \exp(-\frac{(n - n_0)^2}{2 \omega_{0}^2})
MSD analytical solution for off-diagonal disorder for plane wave excitation:
.. math::
    \langle x^2(t) \\rangle = a^2 (\frac{\omega_{0}^2}{2} +  \frac{2 \sigma_J^2}{\hbar^2} t^2)

This program is used to verify if our analytical formula for MSD with off-diagonal disorder matches the numerical solution.
"""

from __future__ import annotations
import os
from pathlib import Path

import numpy as np
import h5py
from matplotlib import pyplot as plt

current_dir = Path(__file__).parent.resolve()

print(f"Current directory: {current_dir}")

project_root = current_dir.parent.parent
print(f"Project root: {project_root}")
breakpoint()

root_path = os.path.dirname(os.path.abspath(__file__))
plane_wave_numerical_data_path = os.path.join(root_path, "NN_cache/nn_chain_sigma_eps0.0_sigma_J0.05_avg.hdf5")

USE_PLANE_WAVE_PHASE = True  # set to False to test the local excitation formula instead
if USE_PLANE_WAVE_PHASE:
    print("Comparing with plane wave excitation analytical formula.")
    numerical_data_path = os.path.join(root_path, "NN_cache/nn_chain_sigma_eps0.0_sigma_J0.05_avg.hdf5")
else:    
    print("Comparing with local excitation analytical formula.")
    numerical_data_path = os.path.join(root_path, "NN_cache/nn_chain_sigma_eps0.0_sigma_J0.05_local_avg.hdf5")
# ---- load numerical data for comparison ----
with h5py.File(numerical_data_path, "r") as f:
    # scalar parameters are stored as file-level attributes
    J_0_eV = float(f.attrs["J_0_eV"])
    sigma_J_eV = float(f.attrs["sigma_J_eV"])
    eps_0_eV = float(f.attrs["eps_0_eV"])
    sigma_eps_eV = float(f.attrs["sigma_eps_eV"])
    t_total_fs = float(f.attrs["t_total_fs"])
    n_steps = int(f.attrs["n_steps"])
    k_parallel = float(f.attrs["k_parallel"])
    sigma_sites = float(f.attrs["sigma_sites"])
    # time axis is a top-level dataset in picoseconds
    t_ps = np.array(f["t_ps"][:])
    # breakpoint()
    # MSD lives inside the "expectations" group
    msd_mean = np.array(f["expectations/msd_mean"][:])


def msd_analytical_formula_offdiag_local_excitation(a, hbar, J_0_eV, sigma_J_eV, t_fs):
    """Analytical formula for MSD with off-diagonal disorder."""
    J_eff_squared = J_0_eV**2 + sigma_J_eV**2
    prefactor = 2 * a**2 * J_eff_squared / hbar**2
    return prefactor * t_fs**2

def msd_analytical_formula_offdiag_plane_wave_excitation(a, hbar, J_0_eV, sigma_J_eV, t_fs, k_parallel, omega_0):
    """Analytical formula for MSD with off-diagonal disorder for plane wave excitation."""
    term1 = omega_0**2 / 2
    term2 = (4 * J_0_eV**2 / hbar**2) * np.sin(k_parallel * a)**2
    term3 = (2 * sigma_J_eV**2) / hbar**2
    prefactor = a**2 * (term1 + ( term3) * t_fs**2)
    return prefactor

def compare_msd_with_analytical_plane_wave_excitation(t_fs, msd_numerical, a, hbar, J_0_eV, sigma_J_eV, k_parallel, omega_0):
    """Compare numerical MSD with analytical formula."""
    if USE_PLANE_WAVE_PHASE:
        msd_analytical = msd_analytical_formula_offdiag_plane_wave_excitation(
            a=a,
            hbar=hbar,
            J_0_eV=J_0_eV,
            sigma_J_eV=sigma_J_eV,
            t_fs=t_fs,
            k_parallel=k_parallel,
            omega_0=omega_0
        )
        
    else:
        msd_analytical = msd_analytical_formula_offdiag_local_excitation(
            a=a,
            hbar=hbar,
            J_0_eV=J_0_eV,
            sigma_J_eV=sigma_J_eV,
            t_fs=t_fs
        )

    fig, ax = plt.subplots(figsize=(8, 6))

    ax.plot(t_fs, msd_numerical, label="Numerical MSD", color="blue", marker="o", linestyle="--")
    ax.plot(t_fs, msd_analytical, label="Analytical MSD", color="red", linestyle="-")
    ax.set_xlabel("Time (fs)")
    ax.set_ylabel("MSD")
    ax.set_xlim(0, 20)
    ax.set_ylim(0, 100)
    ax.legend()
    if USE_PLANE_WAVE_PHASE:
        ax.set_title("Plane wave excitation: Numerical vs Analytical MSD")
    else:
        ax.set_title("Local excitation: Numerical vs Analytical MSD")
    plt.show()

hbar_eVfs = 0.6582119514  # ℏ in eV·fs
a = 1.0                    # lattice constant (lattice units)
omega_0 = sigma_sites      # Gaussian wavepacket width σ_sites (lattice units)

t_fs = np.arange(t_total_fs) # time axis in fs (from the loaded data)

compare_msd_with_analytical_plane_wave_excitation(
    t_fs=t_fs,
    msd_numerical=msd_mean,
    a=a,
    hbar=hbar_eVfs,
    J_0_eV=J_0_eV,
    sigma_J_eV=sigma_J_eV,
    k_parallel=k_parallel,
    omega_0=omega_0,
)
