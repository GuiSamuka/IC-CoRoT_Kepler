from os import system
system("cls")

#%%
from utils import *
import pandas as pd

data = pd.read_csv('resampled_files\RESAMPLED_0100725706_20070516T060226.csv')
time = data.DATE.to_numpy()
flux = data.WHITEFLUX.to_numpy()

curve = lightcurve.LightCurve(time, flux)
# curve.plot()
# curve.view_fourier_spectrum()
curve.median_filter(3).view_fourier_results()



# import numpy as np

# freq = np.fft.fftfreq(curve.flux.shape[-1])
# freq_shifted = np.fft.fftshift(freq)
# sp = np.fft.fft(curve.flux)


# filtered = curve.butterworth_lowpass_filter(9, 0.3)
# freq_filtered = np.fft.fftfreq(filtered.filtered_flux.shape[-1])
# freq_filtered_shifted = np.fft.fftshift(freq_filtered)

# sp_filtered = np.fft.fft(filtered.filtered_flux)


# Bokeh plotting
# from bokeh.io import output_notebook
# from bokeh.plotting import figure, show

# p_real = figure(title='Original and Filtered Fourier Spectrum - Real Part',
#             y_axis_type='log',
#             plot_width=650, plot_height=400,
#             background_fill_color='#fafafa')

# p_real.line(freq, np.real(sp), 
#         line_color='blue', 
#         legend_label='Original spectrum',
#         line_width=2)

# p_real.line(freq_filtered, np.real(sp_filtered), 
#         line_color='red',
#         legend_label='Filtered spectrum',
#         line_width=2)

# p_real.xaxis[0].axis_label = 'Frequency'
# p_real.yaxis[0].axis_label = 'Magnitude'

# output_notebook()
# show(p_real)


# p_imag = figure(title='Original and Filtered Fourier Spectrum - Imag Part',
#             y_axis_type='log',
#             plot_width=650, plot_height=400,
#             background_fill_color='#fafafa')

# p_imag.line(freq, np.imag(sp), 
#         line_color='blue', 
#         legend_label='Original spectrum',
#         line_width=2)

# p_imag.line(freq_filtered, np.imag(sp_filtered), 
#         line_color='red',
#         legend_label='Filtered spectrum',
#         line_width=2)

# p_imag.xaxis[0].axis_label = 'Frequency'
# p_imag.yaxis[0].axis_label = 'Magnitude'

# output_notebook()
# show(p_imag)




# %%
