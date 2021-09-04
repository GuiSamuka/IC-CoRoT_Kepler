"""

    This module defines a PhaseFoldedLightCurve

"""

# Imports 
from utils.visualization.data_viz import view_phase_folded_simulation 

from dataclasses import dataclass, field
import numpy as np
from tabulate import tabulate # Fancy compare results


@dataclass
class PhaseFoldedLightCurve():
    """
    Explicação...

    Parameters
    ----------
    original_time : numpy ndarray
        Original Time values

    original_flux : numpy ndarray
        Original Flux values, related to every time point.

    simulated_time : numpy ndarray
        Simulated/Resampled Time values

    simulated_flux : numpy ndarray
        Simulated/Resampled Flux values, related to every time point.

    chi2 : float
        Error: chi squared

    """

    # Constructor
    original_time: np.array = field(default=None)
    original_flux: np.array = field(default=None)

    simulated_time: np.array = field(default=None)
    simulated_flux: np.array = field(default=None)
    chi2: float = field(default=None)

    # Class Methods 
    def view_simulation_results(self):
        print('Plotting simulation results')
        view_phase_folded_simulation(self.original_time, self.original_flux, self.simulated_time, self.simulated_flux)

    def compare_results(self, see_values=True):
        if see_values:
            print(tabulate(np.c_[self.original_flux, self.simulated_flux], headers=['Original flux', 'Simulated flux'], tablefmt='fancy_grid'))
        print('Chi squared =', round(self.chi2, 4))

