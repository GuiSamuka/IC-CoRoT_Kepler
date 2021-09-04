"""

    This module contains the methods 
    related to the `02 - Filters` notebook.
    It has the implementation of each filtering
    processes used on this research.

"""

# Imports
# from utils.models.filteredlightcurve import FilteredLightCurve

import numpy as np
# from scipy.signal import medfilt
# from math import exp, factorial
from control import *

def expand_edges(array, numExpansion) -> np.ndarray:
    """
    Receive an array and expanded their edges.

    Given an array of length N, it returns 
    an array with length (N + 2*numExpansion).
    The procedure add a certain `numExpansion`
    points before the first point of the array, 
    and the same `numExpansion` after the last 
    point of the array. This values added are
    respectively equals to the first point and 
    the last point of the array.

    Parameters
    ----------
    array : numpy ndarray
        Array to be expanded.

    numExpansion : int
        Number of points to be added at the beginning 
        and at the end of the array.

    Returns
    -------
    `np.array`

    """
    aux_pre = np.zeros(numExpansion)
    aux_pos = np.zeros(numExpansion)
    i = 0

    for i in range(numExpansion):
        aux_pre[i] = array[0]
        aux_pos[i] = array[-1]

    return np.concatenate((aux_pre, array, aux_pos)).ravel()


def padding(array) -> np.ndarray:
    """
    Receive an array and apply the 
    zero padding algorithm.

    Given an array of length N, it returns 
    an array with length (2*N), filled 
    with zeros.

    Parameters
    ----------
    array : numpy ndarray
        Array to be zero padded.

    Returns
    -------
    `np.array`

    """
    return np.append(array, np.zeros(len(array)))


def centralize_fourier(array) -> np.ndarray:
    """
    Receive an array and multiply
    each value of it by the factor:
    (-1)^i, being i the array index.

    Parameters
    ----------
    array : numpy ndarray
        Array to be multiplied

    Returns
    -------
    `np.array`
    """
    multiplied = np.ones(len(array))
    i = 0

    for i in range(len(array)):
        multiplied[i] = array[i] * ((-1)**i)

    return multiplied


def fourier_transform(array) -> np.ndarray:
    """
    Receive an array and computes 
    the fourier transform of it.

    Parameters
    ----------
    array : numpy ndarray
        Array to be transformed

    Returns
    -------
    `np.array`
    """
    return np.fft.fft(array)


def inverse_fourier_transform(array) -> np.ndarray:
    """
    Receive an array and computes 
    the inverse fourier transform of it.

    Parameters
    ----------
    array : numpy ndarray
        Array to be [inversed] transformed

    Returns
    -------
    `np.array`
    """
    return np.real(np.fft.ifft(array))


def remove_padding(array) -> np.ndarray:
    """
    Receive an array and remove
    the padding transformation 
    previously applied.

    Parameters
    ----------
    array : numpy ndarray
        Array to lose their padded

    Returns
    -------
    `np.array`
    
    """
    return array[:int(len(array)/2)]


def remove_expand_edges(array, numExpansion) -> np.ndarray:
    """
    Receive an array and remove the 
    expantion edges.

    The description of what the expanded
    borded do is at `filter_helper.expand_edges`.
    This method undo this procedure.

    Parameters
    ----------
    array : numpy ndarray
        Array to lose the expansion.

    numExpansion : int
        Number of points previously choose.

    Returns
    -------
    `np.array`

    """
    aux = np.delete(array, np.s_[:numExpansion])

    return np.delete(aux, np.s_[-numExpansion:])


# def apply_filter(
#     time, 
#     array, 
#     filter_technique,
#     cutoff_freq=0.1,
#     order=1,
#     numNei=3,
#     numExpansion=70):
    
#     if filter_technique.upper() == 'IDEAL':
#         # Extracting info from curve
#         n_time = len(array)
#         D0 = cutoff_freq*n_time

#         # Procedures to apply the ideal lowpass filter
#         expanded = expand_edges(array, numExpansion=numExpansion)
#         fourier = fourier_transform(expanded)

#         i = 0
#         for i in range(len(fourier)):
#             if fourier[i] > D0:
#                 fourier[i] = 0
        
#         ifft = inverse_fourier_transform(fourier)
#         array_filtered = remove_expand_edges(ifft, numExpansion=numExpansion)
#         array_filtered += (array.mean() - array_filtered.mean())


#     elif filter_technique.upper() == 'MEDIAN':
#         array_filtered = medfilt(array, numNei)
        
#     else:
#         # Extracting info from curve
#         n_time = len(array)
#         D0 = cutoff_freq*n_time
#         xc = n_time

#         # Procedures to filtering on frequency domain
#         expanded = expand_edges(array, numExpansion=numExpansion)
#         padded = padding(expanded)
#         centralized = centralize_fourier(padded)
#         fourier = fourier_transform(centralized)

#         # Creating the low-pass transfer function array
#         len_filter = len(fourier)
#         filter_array = np.zeros(len_filter)

#         if filter_technique.upper() == 'GAUSSIAN':
#             i = 0
#             for i in range(len_filter):
#                 filter_array[i] = exp((-(i-(xc-1.0))**2)/(2*((D0)**2)))

        
#         elif filter_technique.upper() == 'BUTTERWORTH':
#             i = 0
#             for i in range(len_filter):
#                 filter_array[i] = 1.0 / (1.0+(abs(i-(xc-1.0))/D0)**(2.0*order))


#         elif filter_technique.upper() == 'BESSEL':
#             # Coef ak
#             coef = []
#             i = 0
#             while i <= order:
#                 ak = (factorial(2*order - i)) / ( 2**(order - i)*factorial(i)*factorial(order - i) )
#                 coef.append(ak)
#                 i += 1

#             # Computing θn(s)
#             s = TransferFunction.s
#             theta_array = []
#             k=0

#             for k in range(order+1):
#                 theta_n = coef[k] * (s**k)
#                 theta_array.append(theta_n)

#             # Computing H(s)
#             coef_numerator = theta_array[0]
#             list_denominator = theta_array[:]
#             denominator = 0

#             for item in list_denominator:
#                 denominator += item

#             # Computing Transfer Function
#             G = coef_numerator / denominator

#             i=0
#             for i in range(len_filter):
#                 filter_array[i] = np.real(evalfr(G, ( np.abs(i-(xc-1.0))/D0 )))

            
#         raw_filtered = filter_array * fourier
#         ifft_raw_filtered = inverse_fourier_transform(raw_filtered)
#         no_padded_ifft_raw_filtered = remove_padding(ifft_raw_filtered)
#         no_expanded_no_padded_ifft_raw_filtered = remove_expand_edges(
#             no_padded_ifft_raw_filtered, numExpansion)

#         array_filtered = centralize_fourier(no_expanded_no_padded_ifft_raw_filtered)

#         # array_filtered += (array.mean() - array_filtered.mean())
#         array_filtered += (np.mean(array) - np.mean(array_filtered))

#     return FilteredLightCurve(time, array, array_filtered, filter_technique, cutoff_freq, order, numNei)
