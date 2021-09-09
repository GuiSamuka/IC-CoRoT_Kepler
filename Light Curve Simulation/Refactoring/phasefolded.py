"""

    This module defines a PhaseFoldedLightCurve

"""

# Imports
from Refactoring.utils import data_viz

from dataclasses import dataclass, field
import numpy as np
from tabulate import tabulate # Fancy compare results


@dataclass
class PhaseFoldedLightCurve():
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

    original_time: np.ndarray = field(default=None)
    original_flux: np.ndarray = field(default=None)

    simulated_time: np.ndarray = field(default=None)
    simulated_flux: np.ndarray = field(default=None)
    chi2: float = field(default=None)

    def view_simulation_results(self):
        print('Plotting simulation results')
        data_viz.view_phase_folded_simulation(self.original_time, self.original_flux, self.simulated_time, self.simulated_flux)

    def compare_results(self, see_values=True):
        if see_values:
            print(tabulate(np.c_[self.original_flux, self.simulated_flux], headers=['Original flux', 'Simulated flux'], tablefmt='fancy_grid'))
        print('Chi squared =', round(self.chi2, 4))


    