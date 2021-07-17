"""

    This module ...

"""

# Imports
import numpy as np
from numpy import savetxt
from math import exp
import sys
from math import factorial
from control import *
import time
import os
import pandas as pd

class FrequencyDomainFiltering:
    """
        This object ...
    """

    def __init__(self) -> None:
        pass

    def expand_borders(self, array, numExpansion):
        """
            This method ...
        """
        aux_pre = np.zeros(numExpansion)
        aux_pos = np.zeros(numExpansion) 
        i=0        

        for i in range(numExpansion):
            aux_pre[i] = array[0]
            aux_pos[i] = array[-1]

        self.array_expanded = np.concatenate((aux_pre, array, aux_pos)).ravel()

    def padding(self, array):
        """
            This method ...
        """
        self.padded = np.append(array, np.zeros(len(array)))

    def multiplying_by_minus_one_to_index(self, array):
        """
            This method ...
        """
        multiplied = np.ones(len(array))
        i = 0
        
        for i in range(len(array)):
            multiplied[i] = array[i] * ((-1)**i)

        self.multiplied = multiplied

    def fourier_transform(self, array):
        """
            This method ...
        """
        self.fft = np.fft.fft(array)

    def filter_array(self, array, fourier_transform, algorithm, cutoff_freq, order):
        """
            This method ...
        """
        if algorithm.upper() == 'BUTTERWORTH':
            """
            Descrição do filtro
            """
            
            # Extracting features from signal
            n_time = len(array)
            D0 = cutoff_freq * n_time
            xc = n_time

            # Creating the butterworth low-pass transfer function array
            len_filter = len(fourier_transform)
            filter = np.zeros(len_filter)
            i=0

            for i in range(len_filter):
                filter[i] = 1.0 / (1.0+(abs(i-(xc-1.0))/D0)**(2.0*order))
            
            self.array_filter = filter
       

        if algorithm.upper() == 'BESSEL':
            """
            Descrição do filtro
            """

            # Computing ak
            coef = []
            i=0
            
            while i <= order:
                ak = (factorial(2*order - i)) / ( 2**(order - i)*factorial(i)*factorial(order - i) )
                coef.append(ak)
                i+=1

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

            # Extracting features from signal
            n_time = len(array)
            D0 = cutoff_freq * n_time
            xc = n_time

            # Creating the bessel transfer function array
            len_filter = len(fourier_transform)
            filter = np.zeros(len_filter)
            i=0

            for i in range(len_filter):
                filter[i] = np.real(evalfr(G, ( np.abs(i-(xc-1.0))/D0 )))

            self.array_filter = filter


        if algorithm.upper() == 'GAUSSIAN':
            """
            Descrição do filtro
            """
            
            order = None
            # Extracting features from signal
            n_time = len(array)
            D0 = cutoff_freq * n_time
            xc = n_time

            # Creating the gaussian low-pass transfer function array            
            len_filter = len(fourier_transform)
            filter = np.zeros(len_filter)
            i=0
        
            for i in range(len_filter):
                filter[i] = exp( (-(i-(xc-1.0))**2)/(2*((D0)**2)) )

            self.array_filter = filter


        if algorithm.upper() == 'IDEAL':
            """
            No caso do filtro ideal, as componentes dentro de 
            um certo intervalo de frequências são totalmente 
            removidas da Transformada de Fourier dos dados (o que, 
            como é sabido, pode ter como consequência a introdução 
            de artefatos nos dados).
            """
            print("Ideal filtering")
            sys.exit()

    def apply_filter(self, array_filter, fourier_transform):
        """
            This method ...
        """
        self.applied_filter = array_filter * fourier_transform

    def inverse_fourier_transform(self, array):
        """
            This method ...
        """
        self.ifft = np.real(np.fft.ifft(array))

    def remove_padding(self, array):
        """
            This method ...
        """
        self.no_padded = array[:int(len(array)/2)]
    
    def remove_expanded_borders(self, array, numExpansion):
        """
            This method ...
        """ 
        aux = np.delete(array, np.s_[:numExpansion])

        self.no_expanded = np.delete(aux, np.s_[-numExpansion:])

    def filter(self, array, filter_technique, numExpansion, cutoff_freq, order):
        """
            This method ...
        """

        if filter_technique.upper() == 'IDEAL':
            n_time = len(array)
            D0 = cutoff_freq * n_time

            # Expand borders           

            # Padding

            # Fourier transform

            # Filtering
            for i in range(len(y_fourier)):
                if y_fourier[i] > D0:
                    y_fourier[i] = 0
            
            # Inverse Fourier Transform

            # Remove Padding

            # Remove expanded borders

            # self.result = self.no_borders
            



        else: 
            self.expand_borders(array, numExpansion)
            self.padding(self.array_expanded)
            self.multiplying_by_minus_one_to_index(self.padded)
            self.fourier_transform(self.multiplied)
            self.filter_array(array, self.fft, filter_technique, cutoff_freq, order)
            self.apply_filter(self.array_filter, self.fft)
            self.inverse_fourier_transform(self.applied_filter)
            self.remove_padding(self.ifft)
            self.remove_expanded_borders(self.no_padded, numExpansion)
            self.multiplying_by_minus_one_to_index(self.no_expanded)
            self.result = self.multiplied

    
    # Getters
    @property
    def getExpandedBorders(self):
      return self.array_expanded

    @property
    def getZeroPadded(self):
      return self.padded

    @property
    def getFourier(self):
      return self.fft
    
    @property
    def getFiltered(self):
      return self.result
    
    @property
    def getFilterArray(self):
      return self.array_filter

    @property
    def getFilterApplied(self):
      return self.applied_filter

    @property
    def getInverseFourier(self):
      return self.ifft

    @property
    def getNoPadded(self):
      return self.no_padded

    @property 
    def getNoExpanded(self):
      return self.no_expanded


class NonLinearFilter:
    """
        This object ...
    """

    def __init__(self) -> None:
        pass


    def MedianFilter(self):
        pass
    
def export_results_csv(PATH_DIR, filter_technique, cutoff_freq, order):
    print(f"Saving filtered data for {filter_technique}, order = {order} and cutoff frequency = {cutoff_freq}")
    
    # Path to raw csv files
    DATA_DIR = 'C:/Users/guisa/Google Drive/01 - Iniciação Científica/02 - Datasets/csv_files'
    t_o = time.time()
    count = 0
    for root_dir_path, sub_dirs, files in os.walk(DATA_DIR):
        for j in range(0, len(files)):
            if files[j] != 'desktop.ini':
                # print(files[j] + " => Save it!")
                data = pd.read_csv(root_dir_path + "/" + files[j])
                y = data.WHITEFLUX.to_numpy()
                Filter = FrequencyDomainFiltering()
                Filter.filter(array=y, filter_technique=filter_technique, numExpansion=70, cutoff_freq=cutoff_freq, order=order)
                y_filtered = Filter.getFiltered
                y_filtered += (y.mean() - y_filtered.mean())
                savetxt(PATH_DIR + "/" + files[j], y_filtered, delimiter=',', header="WHITEFLUX")
                count += 1
    t_f = time.time()
    print("It takes:", round(t_f-t_o, 2), "seconds to save all")
    print("All files have been save sucessefuly\n") if count == 37 else print("Something went wrong! Please uncomment the line just under the if statement to see details of what file have been not saved\n")

