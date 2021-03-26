import numpy as np
from os import system
system("cls")

class FrequencyDomainFiltering(object):
  def __init__(self):
    pass

  def get_algorithm(self):
    return self._algorithm

  def set_algorithm(self, x):
    self._algorithm = x

  def expand_borders(self, array, numExpansion):
    i = 0
    aux_pre = np.zeros(numExpansion)
    aux_pos = np.zeros(numExpansion)
    
    for i in range(numExpansion):
      aux_pre[i] = array[0]
      aux_pos[i] = array[-1]
    
    self.array_expanded = np.concatenate((aux_pre, array, aux_pos)).ravel()
    
  def padding(self, array):
    self.padded = np.append(array, np.zeros(len(array)))

  def multiplying_by_minus_one_to_index(self, array):
    i = 0
    aux = np.ones(len(array))

    for i in range(len(array)):
      aux[i] = array[i] * ((-1)**i)
    self.multiplied = aux

  def fourier_transform(self, array):
    self.fft = np.fft.fft(array)

  def filter_technique(self, array, algorithm):
    if algorithm.upper() == 'BUTTERWORTH':
      print("Butterworth filtering")

  def inverse_fourier_transform(self, array):
    self.ifft = np.real(np.fft.ifft(array))

  def no_padding(self, array):
    self.no_padded = array[:int(len(array)/2)]

  def remove_expanded_borders(self, array, numExpansion):
    aux = np.delete(array, np.s_[:numExpansion])
    self.no_expanded = np.delete(aux, np.s_[-numExpansion:])
    self.filtered = self.no_expanded

   
    

  def filter(self, array, filter, numExpansion):
    self.expand_borders(array, numExpansion)
    self.padding(self.array_expanded)
    self.fourier_transform(self.padded)

    self.filter_technique(self.fft, filter)

    self.inverse_fourier_transform(self.fft)
    self.no_padding(self.ifft)
    self.remove_expanded_borders(self.no_padded, numExpansion)

  

y = [1, 2, 3, 4]
filter = FrequencyDomainFiltering()
filter.filter(y, 'butterworth', 3)
print(filter)