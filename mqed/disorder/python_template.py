# -*- coding: utf-8 -*-
import numpy as np
import h5py
import scipy.io
import time
from scipy.sparse import diags, csc_matrix
from scipy.sparse.linalg import expm_multiply

start_time = time.time()

# ---- constants ----
hbar = 1.0545718e-34
eV = 1.602176634e-19

# ---- simulation controls ----
N_Emit = 1000
Tot_t_start = 0.0
Seg_t_end = 1e-14          # seconds per segment
Re_t_step = 2e-19          # seconds (physical dt grid, only for reporting)
Num_segment = 150
Record_step = 10000        # record every this many "grid points"

# ---- Gaussian wavepacket initial condition ----
USE_PLANE_WAVE_PHASE = True
k_parallel = 1.0           # radians per site (set 0.0 for purely real Gaussian)
sigma_sites = 10.0         # Gaussian width in site units (std dev)
# If you prefer your paper's exp(-(x^2)/w0^2), then sigma_sites = w0 / sqrt(2)

# Derived counts
N_step = int((Seg_t_end - Tot_t_start) / Re_t_step) + 1
rec_idx = np.arange(0, N_step, Record_step, dtype=int)
n_rec = len(rec_idx)

# Precompute time grids for a segment (same for every segment)
t_span_segment = np.linspace(Tot_t_start, Tot_t_start + Seg_t_end, num=N_step)
t_rec_segment = t_span_segment[rec_idx]

# Dimensionless tau = t * eV/hbar (matches your old code convention)
tau_span_segment = t_span_segment * eV / hbar
tau_rec_segment = tau_span_segment[rec_idx]  # used only to get dtau between records

# ---- load eps and J from MATLAB v7.3 file ----
matfile = '_MAT_FILE_PATH_'   # <-- change to your filename
with h5py.File(matfile, "r") as f:
    eps = np.array(f["eps"]).squeeze()
    J = np.array(f["J"]).squeeze()

eps = eps.astype(np.complex128)  # allow complex
J = J.astype(np.complex128)

N = eps.size
assert N == N_Emit, f"eps length {N} != N_Emit {N_Emit}"
assert J.size == N_Emit - 1, f"J length {J.size} != N_Emit-1 {N_Emit-1}"

# ---- Build Hamiltonian (in eV) ----
H_phys = diags(
    diagonals=[np.conjugate(J), eps, J],
    offsets=[-1, 0, 1],
    format="csc"
)

# ---- Build generator for dimensionless tau ----
# dC/dtau = (-i * H[eV]) C
H_eff = csc_matrix((-1j) * H_phys)

# ---- initial condition: Gaussian wavepacket centered in the middle ----
Index_Ini_Exc = (N_Emit - 1) // 2  # middle site
n = np.arange(N_Emit, dtype=float)
x0 = n - float(Index_Ini_Exc)      # coordinate relative to center

# Gaussian envelope (std dev = sigma_sites)
envelope = np.exp(-(x0 ** 2) / (2.0 * sigma_sites ** 2))

if USE_PLANE_WAVE_PHASE:
    phase = np.exp(1j * k_parallel * x0)
else:
    phase = 1.0

C0 = (envelope * phase).astype(np.complex128)

# Normalize
norm = np.sqrt(np.sum(np.abs(C0) ** 2))
if norm == 0:
    raise ValueError("Initial Gaussian wavepacket has zero norm. Check sigma_sites.")
C0 /= norm

# ---- MSD operator ----
# Define x relative to the *initial center* (middle) to avoid boundary bias.
x = np.arange(N_Emit) - Index_Ini_Exc
Inter_Mol_Distance = 1.0
x_squared = (Inter_Mol_Distance * x) ** 2

# Optional: also track <x> if you want
# x_op = Inter_Mol_Distance * x

# ---- Propagate segment by segment (record only) ----
for seg in range(Num_segment):
    print("segment", seg)

    # Physical time grid for this segment
    Curr_t_start = Tot_t_start + seg * Seg_t_end
    Curr_t_end = Tot_t_start + (seg + 1) * Seg_t_end

    # Reuse segment-relative grid but shift by Curr_t_start
    t_span = np.linspace(Curr_t_start, Curr_t_end, num=N_step)
    t_rec = t_span[rec_idx]

    # tau grid (dimensionless) for this segment
    tau_span = t_span * eV / hbar
    tau_rec = tau_span[rec_idx]

    # Storage only for recorded points
    C_rec = np.zeros((N_Emit, n_rec), dtype=np.complex128)
    MSD_rec = np.zeros(n_rec, dtype=float)
    # Xmean_rec = np.zeros(n_rec, dtype=float)  # if needed

    # start of this segment
    C = C0.copy()
    C_rec[:, 0] = C
    P = np.abs(C) ** 2
    MSD_rec[0] = np.sum(P * x_squared)
    # Xmean_rec[0] = np.sum(P * x_op)

    # propagate between record points using expm_multiply
    for k in range(1, n_rec):
        dtau = tau_rec[k] - tau_rec[k - 1]
        C = expm_multiply(H_eff * dtau, C)
        C_rec[:, k] = C

        P = np.abs(C) ** 2
        MSD_rec[k] = np.sum(P * x_squared)
        # Xmean_rec[k] = np.sum(P * x_op)

    # Save like your previous workflow
    scipy.io.savemat(f'./QD_Time_{seg}.mat', {'Py_time': t_rec})
    scipy.io.savemat(f'./QD_C_coet_{seg}.mat', {'Py_C_coet': C_rec})
    scipy.io.savemat(f'./QD_Init_Gaussian.mat', {
        'k_parallel': float(k_parallel),
        'sigma_sites': float(sigma_sites),
        'Index_Ini_Exc': int(Index_Ini_Exc),
        'C0': C0
    })
    # Optional MSD output
    # scipy.io.savemat(f'./QD_MSD_{seg}.mat', {'Py_MSD': MSD_rec})

    # pass final state to next segment
    C0 = C.copy()

    print("--- %.2f seconds ---" % (time.time() - start_time))