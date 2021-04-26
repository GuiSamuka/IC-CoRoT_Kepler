import numpy as np
class FrequencyDomainFiltering(object):
  def __init__(self):
    pass
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
    
  def filter_array(self, array, fourier_transform, algorithm, cutoff_freq, order):
    if algorithm.upper() == 'BUTTERWORTH':
      print("Butterworth filtering")
      # Extracting information of the signal
      n_time = len(array)
      D0 = cutoff_freq * n_time
      xc = n_time
      # Creating the filter array
      len_filter = len(fourier_transform)
      filter = np.zeros(len_filter)
      for i in range(len_filter):
        filter[i] = 1.0 / (1.0+(abs(i-(xc-1.0))/D0)**(2.0*order))
      # self.butter = filter * fourier_transform
      self.filter_array = filter
      
  def apply_filter(self, array_filter, fourier_transform):
    self.applied_filter = array_filter * fourier_transform
  def inverse_fourier_transform(self, array):
    self.ifft = np.real(np.fft.ifft(array))
  def no_padding(self, array):
    self.no_padded = array[:int(len(array)/2)]
  def remove_expanded_borders(self, array, numExpansion):
    aux = np.delete(array, np.s_[:numExpansion])
    self.no_expanded = np.delete(aux, np.s_[-numExpansion:])
    
  def filter(self, array, filter, numExpansion, cutoff_freq, order):
    self.expand_borders(array, numExpansion)
    self.padding(self.array_expanded)
    self.multiplying_by_minus_one_to_index(self.padded)
    self.fourier_transform(self.multiplied)
    #def filter_technique(self, array, fourier_transform, algorithm, cutoff_freq, order):
    self.filter_array(array, self.fft, filter, cutoff_freq, order)
    #def apply_filter(self, array_filter, fourier_transform):
    self.apply_filter(self.filter_array, self.fft)
    self.inverse_fourier_transform(self.applied_filter)
    self.no_padding(self.ifft)
    self.remove_expanded_borders(self.no_padded, numExpansion)
    self.multiplying_by_minus_one_to_index(self.no_expanded)
  
  @property #TODO Create getter for all attributes
  def filtered(self):
    return self.multiplied
  @property
  def filter_response(self):
    return self.filter_array