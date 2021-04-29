"""

    This module ...

"""

# Imports
import numpy as np
from math import exp
import sys


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
            #print("Butterworth filtering")

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
            #print("Bessel filtering")
            sys.exit()


        if algorithm.upper() == 'GAUSSIAN':
            """
            O filtro gaussiano possui apenas um parâmetro livre, 
            que é a frequência de corte da filtragem (que basicamente 
            dá a largura do filtro)
            """
            #print("Gaussian filtering")
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
                filter[i] = exp( (-(i-(xc-1.0))**2)/(2*((cutoff_freq * n_time)**2)) )

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
    def getFiltered(self):
        return self.result
    
    @property
    def getFourier(self):
        return self.fft

    @property
    def getFilterArray(self):
        return self.array_filter

    @property
    def getInverseFourier(self):
        return self.ifft


class NonLinearFilter:
    """
        This object ...
    """

    def __init__(self) -> None:
        pass


    def MedianFilter(self):
        pass
    