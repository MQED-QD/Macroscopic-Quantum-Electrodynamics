import pandas as pd
import numpy as np
from pathlib import Path
from scipy.interpolate import interp1d
from loguru import logger
from importlib import resources

# Physical constants for converting wavelength to angular frequency
C_METERS_PER_SEC = 2.99792458e8  # Speed of light in m/s

class DataProvider:
    """
    Handles loading and providing material permittivity. It is optimized to create
    interpolation functions only once and performs interpolation in frequency
    space for better numerical accuracy.
    """

    def __init__(self, material_config):
        self.config = material_config
        self.source_type = self.config.source_type
        logger.info(f"Initializing DataProvider with source type: '{self.source_type}'")

        if self.source_type == 'excel':
            self._setup_interpolator_from_excel()
        elif self.source_type == 'constant':
            self._setup_constant_epsilon()
        else:
            raise ValueError(f"Unknown material source_type: '{self.source_type}'")

    def _setup_interpolator_from_excel(self):
        """Loads data, converts to frequency space, and creates interpolators."""
        try:
            excel_cfg = self.config.excel_config
            filepath = resources.files('mqed') / excel_cfg.filepath
            logger.info(f"Loading dispersive material data from: {filepath}")
            logger.info(f"Excel sheet name: {excel_cfg.sheet_name}")

            df = pd.read_excel(filepath, sheet_name=excel_cfg.sheet_name)
            
            lambda_nm = df.iloc[:, 0].values
            epsilon_real = df.iloc[:, 1].values
            epsilon_imag = df.iloc[:, 2].values
            epsilon_complex = epsilon_real + 1j * epsilon_imag

            # Convert wavelength to angular frequency (omega)
            omega_data = 2 * np.pi * C_METERS_PER_SEC / (lambda_nm * 1e-9)

            # Sort data by omega, as interp1d requires monotonically increasing x-values
            sort_indices = np.argsort(omega_data)
            omega_sorted = omega_data[sort_indices]
            epsilon_sorted = epsilon_complex[sort_indices]

            # Create and store interpolation functions once
            self.interp_real = interp1d(
                omega_sorted, epsilon_sorted.real, kind='cubic',
                bounds_error=False, fill_value="extrapolate"
            )
            self.interp_imag = interp1d(
                omega_sorted, epsilon_sorted.imag, kind='cubic',
                bounds_error=False, fill_value="extrapolate"
            )
            logger.success("Successfully initialized interpolator from Excel data.")

        except Exception as e:
            logger.error(f"Failed to load or process Excel data: {e}")
            raise

    def _setup_constant_epsilon(self):
        """Sets up the provider for a non-dispersive material."""
        self.constant_epsilon = complex(self.config.constant_value)
        logger.info(f"Using constant, non-dispersive epsilon = {self.constant_epsilon}")

    def get_epsilon(self, omega):
        """
        Returns the complex relative permittivity (ε_r) for a given angular frequency.

        Args:
            omega (float): The target angular frequency in rad/s.

        Returns:
            complex: The complex relative permittivity ε_r(ω).
        """
        if self.source_type == 'excel':
            real_part = self.interp_real(omega)
            imag_part = self.interp_imag(omega)
            return complex(real_part, imag_part)
        
        elif self.source_type == 'constant':
            return self.constant_epsilon