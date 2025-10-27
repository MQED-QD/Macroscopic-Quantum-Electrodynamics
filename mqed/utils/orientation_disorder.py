import numpy as np


def phi_wrapped_normal_deg(N, mu_deg, sigma_deg, seed=None):
    rng = np.random.default_rng(seed)
    return np.mod(rng.normal(mu_deg, sigma_deg, size=N), 360.0)