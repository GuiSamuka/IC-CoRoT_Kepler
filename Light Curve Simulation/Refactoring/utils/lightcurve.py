"""

    This module defines a LightCurve

"""

from . import data_viz
from . import filter_helper
import os
import pandas as pd
import numpy as np

__all__ = ["LightCurve"]


class LightCurve():
    """
    Explicação...

    Parameters
    ----------
    time : numpy ndarray
        Time values.

    flux : numpy ndarray
        Flux values, related to every time point.

    """

    # Constructor of the `LightCurve` class
    _required_columns = ["time", "flux"]

    def __init__(self, time=None, flux=None):
        """
        Class constructor
        """
        self.time = time
        self.flux = flux

    def plot(self, title='Lightcurve', x_axis='Julian Date', y_axis='Flux', label='Lightcurve'):
        """
        Plot the LightCurve using utils.Data_Viz's `utils.data_viz.view_lightcurve` method. 
        """
        data_viz.view_lightcurve(self.time, self.flux, title=title, x_axis=x_axis, y_axis=y_axis, label=label)
    
    def view_fourier_spectrum(self):
        """
        Computes the Fourier Spectrum, and plot it using utils.Data_Viz's 
        `utils.data_viz.view_fourier` method.
        """
        freq = np.fft.fftfreq(self.flux.shape[-1])
        sp = np.fft.fft(self.flux)

        data_viz.view_fourier(freq_data=freq, sp_data=sp)

    def ideal_lowpass_filter(self, cutoff_freq, numExpansion=70):
        """
        Realize the ideal lowpass filter process using `filter_helper.apply_filter` method.
        """
        return filter_helper.apply_filter(self.time, self.flux, filter_technique='ideal', cutoff_freq=cutoff_freq, numExpansion=numExpansion)
        
    def gaussian_lowpass_filter(self, cutoff_freq, numExpansion=70):
        """
        Realize the gaussian lowpass filter process using `filter_helper.apply_filter` method.
        """
        return filter_helper.apply_filter(self.time, self.flux, filter_technique='gaussian', cutoff_freq=cutoff_freq, numExpansion=numExpansion)

    def butterworth_lowpass_filter(self, order, cutoff_freq, numExpansion=70):
        """
        Realize the butterworth lowpass filter process using `filter_helper.apply_filter` method.
        """
        return filter_helper.apply_filter(self.time, self.flux, filter_technique='butterworth', order=order, cutoff_freq=cutoff_freq, numExpansion=numExpansion)

    def bessel_lowpass_filter(self, order, cutoff_freq, numExpansion=70):
        """
        Realize the bessel lowpass filter process using `filter_helper.apply_filter` method.
        """
        return filter_helper.apply_filter(self.time, self.flux, filter_technique='bessel', order=order, cutoff_freq=cutoff_freq, numExpansion=numExpansion)

    def median_filter(self, numNei, numExpansion=70):
        """
        Realize the median filter process using `filter_helper.apply_filter` method.
        """
        return filter_helper.apply_filter(self.time, self.flux, filter_technique='median', numNei=numNei, numExpansion=numExpansion)

    def how_to_filter(self, order: int, cutoff_freq: float, numNei: int, numExpansion: int):
        """

        This function describes how to filtering using this library

        Parameters
        ----------
        order : int
            Used in Butterworth and Bessel filtering. Matches the filter order. 

        cutoff_freq : float
            Used in Ideal, Gaussian, Butterworth and Bessel filtering. Matches the cutoff frequency.

        numNei : int
            Used in Median. Matches the number of neighbors to consider

        numExpansion : int
            Used in all processes. Corresponds to how much you want to expanded the curve's edges 
            (to avoid some problems caused by the Fast Fourier Transform algorithm). 
            Preliminary tests show that all processes works fine with numExpansion=70, 
            except for Bessel filtering, which required numExpansion=100 
          

        Methods
        -------
        curve.ideal_lowpass_filter(cutoff_freq)

        curve.gaussian_lowpass_filter(cutoff_freq)

        curve.butterworth_lowpass_filter(order, cutoff_freq)

        curve.bessel_lowpass_filter(order, cutoff_freq)

        curve.median_filter(numNei)


        Examples
        --------
        >>> from utils import *
        >>> curve = lightcurve.Lightcurve(time=[1, 2, 3, 4], flux=[10, 15, 12, 20])
        >>> print(curve)
        <utils.lightcurve.LightCurve object>
        >>> curve.plot()

        >>> filtered = curve.bessel_lowpass_filter(order=3, cutoff_freq=0.4, numExpansion=100)
        >>> print(filtered)
        <utils.lightcurve.FilteredLightCurve object>

        >>> filtered.view_filter_results()

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
        Plot the LightCurve using utils.Data_Viz's `utils.data_viz.view_filter_results` method. 
        """
        data_viz.view_results(self.time, self.flux, self.filtered_flux, self.filter_technique, self.cutoff_freq, self.order, self.numNei)

    def view_fourier_results(self):
        """
        Computes the Fourier Spectrum of the oringal and the filtered curve and plot it, using utils.Data_Viz's `utils.data_viz.compare_fourier_original_filtered` method.
        """
        freq_original = np.fft.fftfreq(self.flux.shape[-1])
        sp_original = np.fft.fft(self.flux)

        freq_filtered = np.fft.fftfreq(self.filtered_flux.shape[-1])
        sp_filtered = np.fft.fft(self.filtered_flux)

        data_viz.compare_fourier_original_filtered(freq_original, sp_original, freq_filtered, sp_filtered)
        


    
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
