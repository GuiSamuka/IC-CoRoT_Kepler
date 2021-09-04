"""

    This module defines a LightCurve

"""

# Imports
from utils.models.filteredlightcurve import FilteredLightCurve
from utils.visualization.data_viz import view_lightcurve, view_fourier
from utils.filtering.filter_helper import *

from dataclasses import dataclass, field
import numpy as np
from scipy.signal import medfilt
from math import exp, factorial

@dataclass
class LightCurve():
    """
    Explicação...

    Parameters
    ----------
    time : numpy array
        Time values.

    flux : numpy array
        Flux values, related to every time point.

    """

    # Constructor
    time: np.array = field(default=None)
    flux: np.array = field(default=None)

    # Class Methods 
    def plot(self, title='Lightcurve', x_axis='Julian Data', y_axis='Flux', label='Lightcurve'):
        view_lightcurve(self.time, self.flux, title=title, x_axis=x_axis, y_axis=y_axis, label=label)

    def view_fourier_spectrum(self):
        freq = np.fft.fftfreq(self.flux.shape[-1])
        sp   = np.fft.fft(self.flux)

        view_fourier(freq_data=freq, sp_data=sp)

    def apply_filter(
        self,
        time, 
        array, 
        filter_technique,
        cutoff_freq=0.1,
        order=1,
        numNei=3,
        numExpansion=70):
        pass

        if filter_technique.upper() == 'IDEAL':
            # Extracting info from curve
            n_time = len(array)
            D0 = cutoff_freq*n_time

            # Procedures to apply the ideal lowpass filter
            expanded = expand_edges(array, numExpansion=numExpansion)
            fourier = fourier_transform(expanded)

            i = 0
            for i in range(len(fourier)):
                if fourier[i] > D0:
                    fourier[i] = 0
            
            ifft = inverse_fourier_transform(fourier)
            array_filtered = remove_expand_edges(ifft, numExpansion=numExpansion)
            array_filtered += (array.mean() - array_filtered.mean())


        elif filter_technique.upper() == 'MEDIAN':
            array_filtered = medfilt(array, numNei)
            
        else:
            # Extracting info from curve
            n_time = len(array)
            D0 = cutoff_freq*n_time
            xc = n_time

            # Procedures to filtering on frequency domain
            expanded = expand_edges(array, numExpansion=numExpansion)
            padded = padding(expanded)
            centralized = centralize_fourier(padded)
            fourier = fourier_transform(centralized)

            # Creating the low-pass transfer function array
            len_filter = len(fourier)
            filter_array = np.zeros(len_filter)

            if filter_technique.upper() == 'GAUSSIAN':
                i = 0
                for i in range(len_filter):
                    filter_array[i] = exp((-(i-(xc-1.0))**2)/(2*((D0)**2)))

            
            elif filter_technique.upper() == 'BUTTERWORTH':
                i = 0
                for i in range(len_filter):
                    filter_array[i] = 1.0 / (1.0+(abs(i-(xc-1.0))/D0)**(2.0*order))


            elif filter_technique.upper() == 'BESSEL':
                # Coef ak
                coef = []
                i = 0
                while i <= order:
                    ak = (factorial(2*order - i)) / ( 2**(order - i)*factorial(i)*factorial(order - i) )
                    coef.append(ak)
                    i += 1

                # Computing θn(s)
                s = TransferFunction.s
                theta_array = []
                k=0

                for k in range(order+1):
                    theta_n = coef[k] * (s**k)
                    theta_array.append(theta_n)

                # Computing H(s)
                coef_numerator = theta_array[0]
                list_denominator = theta_array[:]
                denominator = 0

                for item in list_denominator:
                    denominator += item

                # Computing Transfer Function
                G = coef_numerator / denominator

                i=0
                for i in range(len_filter):
                    filter_array[i] = np.real(evalfr(G, ( np.abs(i-(xc-1.0))/D0 )))

                
            raw_filtered = filter_array * fourier
            ifft_raw_filtered = inverse_fourier_transform(raw_filtered)
            no_padded_ifft_raw_filtered = remove_padding(ifft_raw_filtered)
            no_expanded_no_padded_ifft_raw_filtered = remove_expand_edges(
                no_padded_ifft_raw_filtered, numExpansion)

            array_filtered = centralize_fourier(no_expanded_no_padded_ifft_raw_filtered)

            # array_filtered += (array.mean() - array_filtered.mean())
            array_filtered += (np.mean(array) - np.mean(array_filtered))

        return FilteredLightCurve(time, array, array_filtered, filter_technique, cutoff_freq, order, numNei)

    def ideal_lowpass_filter(self, cutoff_freq, numExpansion=70):
        return self.apply_filter(self.time, self.flux, filter_technique='ideal', cutoff_freq=cutoff_freq, numExpansion=numExpansion)

    def gaussian_lowpass_filter(self, cutoff_freq, numExpansion=70):
        return self.apply_filter(self.time, self.flux, filter_technique='gaussian', cutoff_freq=cutoff_freq, numExpansion=numExpansion)

    def butterworth_lowpass_filter(self, order, cutoff_freq, numExpansion=70):
        return self.apply_filter(self.time, self.flux, filter_technique='butterworth', order=order, cutoff_freq=cutoff_freq, numExpansion=numExpansion)

    def bessel_lowpass_filter(self, order, cutoff_freq, numExpansion=70):
        return self.apply_filter(self.time, self.flux, filter_technique='bessel', order=order, cutoff_freq=cutoff_freq, numExpansion=numExpansion)

    def median_filter(self, numNei, numExpansion=70):
        return self.apply_filter(self.time, self.flux, filter_technique='median', numNei=numNei, numExpansion=numExpansion)
    
    def how_to_filter(self, order: int, cutoff_freq: float, numNei: int, numExpansion: int): #TODO
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
        

        """
        pass

    


