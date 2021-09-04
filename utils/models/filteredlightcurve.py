"""

    This module defines a FilteredLightCurve

"""

# Imports 
from utils.visualization.data_viz import view_results, compare_fourier_original_filtered

from dataclasses import dataclass, field
import numpy as np

@dataclass
class FilteredLightCurve():
    """
    Explicação...

    Parameters
    ----------
    time : numpy ndarray
        Time values

    flux : numpy ndarray
        Flux values, related to every time point.

    filtered_flux : numpy ndarray
        Filtered flux values, related to every time point.

    """
    
    # Constructor
    time: np.array = field(default=None)
    flux: np.array = field(default=None)
    filtered_flux: np.array = field(default=None)
    filter_technique: str = field(default=None)
    cutoff_freq: float = field(default=None)
    order: int = field(default=None)
    numNei: int = field(default=None)

    # Class Methods
    def view_filter_results(self):
        view_results(self.time, self.flux, self.filtered_flux, self.filter_technique, self.cutoff_freq, self.order, self.numNei)

    def view_fourier_results(self):
        freq_original = np.fft.fftfreq(self.flux.shape[-1])
        sp_original = np.fft.fft(self.flux)

        freq_filtered = np.fft.fftfreq(self.filtered_flux.shape[-1])
        sp_filtered = np.fft.fft(self.filtered_flux)

        compare_fourier_original_filtered(freq_original, sp_original, freq_filtered, sp_filtered)

    def getFilteredFlux(self):
        return self.filtered_flux







