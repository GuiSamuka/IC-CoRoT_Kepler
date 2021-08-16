"""

    This module defines a LightCurve

"""
from . import data_viz
from . import filter_helper
import os
import pandas as pd

__all__ = ["LightCurve"]


class LightCurve():
    """
    Explicação...

    Parameters
    ----------
    time : numpy ndarray
        Time values

    flux : numpy ndarray
        Flux values, related to every time point.

    """

    # Constructor of the `LightCurve` class
    _required_columns = ["time", "flux"]

    def __init__(self, time=None, flux=None):
        self.time = time
        self.flux = flux

    def plot(self):
        """
        Plot the LightCurve using utils.Data_Viz's `~utils.data_viz.view_lightcurve` method. 
        """
        data_viz.view_lightcurve(self.time, self.flux)

    def ideal_lowpass_filter(self, cutoff_freq, numExpansion=70):
        return filter_helper.apply_filter(self.time, self.flux, filter_technique='ideal', cutoff_freq=cutoff_freq, numExpansion=numExpansion)
        
    def gaussian_lowpass_filter(self, cutoff_freq, numExpansion=70):
        return filter_helper.apply_filter(self.time, self.flux, filter_technique='gaussian', cutoff_freq=cutoff_freq, numExpansion=numExpansion)

    def butterworth_lowpass_filter(self, order, cutoff_freq, numExpansion=70):
        return filter_helper.apply_filter(self.time, self.flux, filter_technique='butterworth', order=order, cutoff_freq=cutoff_freq, numExpansion=numExpansion)

    def bessel_lowpass_filter(self, order, cutoff_freq, numExpansion=70):
        return filter_helper.apply_filter(self.time, self.flux, filter_technique='bessel', order=order, cutoff_freq=cutoff_freq, numExpansion=numExpansion)

    def median_filter(self, numNei, numExpansion=70):
        return filter_helper.apply_filter(self.time, self.flux, filter_technique='median', numNei=numNei, numExpansion=numExpansion)

    def how_to_filter(self, order, cutoff_freq, numNei, numExpansion=70):
        """
        This function describes how to filtering

        Parameters
        ----------
        order:
        cutoff_freq:
        numNei:
        numExpansion:
        """
        pass


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

    # Constructor of the `LightCurve` class
    _required_columns = ["time", "flux", "filtered_flux", "filter_technique", "cutoff_freq", "order", "numNei"]

    def __init__(self, time=None, flux=None, filtered_flux=None, filter_technique="", cutoff_freq=None, order=None, numNei=None):
        self.time = time
        self.flux = flux
        self.filtered_flux = filtered_flux
        self.filter_technique = filter_technique
        self.cutoff_freq = cutoff_freq
        self.order = order
        self.numNei = numNei
        

    def view_filter_results(self):
        """
        Plot the LightCurve using utils.Data_Viz's `~utils.data_viz.view_filter_results` method. 
        """
        data_viz.view_results(self.time, self.flux, self.filtered_flux, self.filter_technique, self.cutoff_freq, self.order, self.numNei)

    def getFilteredFlux(self):
        return self.filtered_flux



def export_results_csv(WHERE_TO_SAVE_PATH, filter_technique, cutoff_freq, order, numNei):
    # Path to resampled csv files
    DATASET_PATH = 'C:/Users/guisa/Google Drive/01 - Iniciação Científica/IC-CoRoT_Kepler/resampled_files'
    
    count = 0
    for root_dir_path, sub_dirs, files in os.walk(DATASET_PATH):
        for j in range(0, len(files)):
            if files[j].endswith('.csv'):
                # print(files[j] + " => Save it!")
                data = pd.read_csv(root_dir_path + "/" + files[j])

                time = data.DATE.to_numpy()
                flux = data.WHITEFLUX.to_numpy()

                curve = LightCurve(time, flux)
                

                if filter_technique.upper() == 'IDEAL':
                    filtered = curve.ideal_lowpass_filter(cutoff_freq)
                    flux_filtered = filtered.getFilteredFlux()

                elif filter_technique.upper() == 'GAUSSIAN':
                    filtered = curve.gaussian_lowpass_filter(cutoff_freq)
                    flux_filtered = filtered.getFilteredFlux()

                elif filter_technique.upper() == 'BUTTERWORTH':
                    filtered = curve.butterworth_lowpass_filter(order, cutoff_freq)
                    flux_filtered = filtered.getFilteredFlux()

                elif filter_technique.upper() == 'BESSEL':
                    filtered = curve.bessel_lowpass_filter(order, cutoff_freq, numExpansion=100)
                    flux_filtered = filtered.getFilteredFlux()  
                
                elif filter_technique.upper() == 'MEDIAN':
                    filtered = curve.median_filter(numNei)
                    flux_filtered = filtered.getFilteredFlux()

                
                # Creating a new `pd.DataFrame`
                concat_dict = {
                  "DATE": pd.Series(time), 
                  "WHITEFLUX": pd.Series(flux_filtered)
                }
                curve_filtered = pd.concat(concat_dict, axis=1)

                # Salving data
                curve_filtered.to_csv(WHERE_TO_SAVE_PATH + "/" + files[j], index=False)
                count += 1
    # print("All files have been saved sucessefuly\n") if count == 33 else print("Something went wrong! Please uncomment the line just under the if statement to see details of what file have been not saved\n")
