"""

    This module contains the methods 
    related to visualization

"""


# Imports
from bokeh.plotting import figure, output_file, show
from bokeh.io import output_notebook
# import pandas as pd
import numpy as np

# output_file('image.html')
output_notebook()


def view_lightcurve(
  x_data=None, 
  y_data=None,
  title='LightCurve',
  x_axis='Julian Date',
  y_axis='Flux',
  label='Lightcurve'):

  p = figure(title=title,
            plot_width=650, plot_height=400)

  p.xaxis[0].axis_label = x_axis
  p.yaxis[0].axis_label = y_axis

  p.line(x_data, y_data, line_width=2, legend_label=label)

  show(p)

def view_results(
  time=None, 
  flux=None,
  filtered_flux=None,
  filter_technique="",
  cutoff_freq=None,
  order=None,
  numNei=None,
  title='LightCurve',
  x_axis='Julian Date',
  y_axis='Flux'):

  if filter_technique.upper() == 'IDEAL':
    title = f"{filter_technique.capitalize()} filter with Cutoff frequency = {cutoff_freq}"

  elif filter_technique.upper() == 'MEDIAN':
    title = f"{filter_technique.capitalize()} filter with {numNei} neighbors"

  elif filter_technique.upper() == 'GAUSSIAN':
    title = f"{filter_technique.capitalize()} filter with Cutoff frequency = {cutoff_freq}"

  elif filter_technique.upper() == 'BUTTERWORTH':
    title = f"{filter_technique.capitalize()} filter with Order = {order} and Cutoff frequency = {cutoff_freq}"

  elif filter_technique.upper() == 'BESSEL':
    title = f"{filter_technique.capitalize()} filter with Order = {order} and Cutoff frequency = {cutoff_freq}"


  p = figure(title=title,
            plot_width=650, plot_height=400)

  p.xaxis[0].axis_label = x_axis
  p.yaxis[0].axis_label = y_axis

  xs = [time, time]
  ys = [flux, filtered_flux]

  p.multi_line(xs, ys, color=["blue", "red"], line_width=2)

  show(p)

def view_fourier(
  freq_data=None,
  sp_data=None,
  title='Original Fourier Spectrum'):

  p = figure(title=title,
        y_axis_type='log',
        plot_width=650, plot_height=400,
        background_fill_color='#fafafa')

  p.line(freq_data, np.real(sp_data),
        line_color='blue',
        legend_label='Original Spectrum',
        line_width=2)

  p.xaxis[0].axis_label = 'Frequency'
  p.yaxis[0].axis_label = 'Magnitude'

  output_notebook()
  show(p)

def compare_fourier_original_filtered(
  freq_original_data=None,
  sp_original_data=None,
  freq_filtered_data=None,
  sp_filtered_data=None,
  title='Original and Filtered Fourier Spectrum'):
  
  # print('Original frequency data: ', freq_original_data)
  # print('Original sp data: ', sp_original_data)
  # print('Filtered frequency data: ', freq_filtered_data)
  # print('Filtered sp data:', sp_filtered_data)

  p = figure(title=title, 
        y_axis_type='log',
        plot_width=650, plot_height=400,
        background_fill_color='#fafafa')
  
  p.line(freq_original_data, np.real(sp_original_data),
        line_color='blue',
        legend_label='Original Spectrum',
        line_width=2)
  
  p.line(freq_filtered_data, np.real(sp_filtered_data),
        line_color='red',
        legend_label='Filtered Spectrum',
        line_width=2)

  p.xaxis[0].axis_label = 'Frequency'
  p.yaxis[0].axis_label = 'Magnitude'

  output_notebook()
  show(p)

