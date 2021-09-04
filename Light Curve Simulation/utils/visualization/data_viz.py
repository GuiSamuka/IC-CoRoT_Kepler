"""

    This module contains the methods 
    related to visualization

"""

# Imports
from bokeh.plotting import figure, show#, output_file
from bokeh.io import output_notebook
from bokeh.models import Legend, LegendItem

import numpy as np
from scipy.stats import norm
import scipy.special


output_notebook()

def view_lightcurve(
    x_data=None,
    y_data=None,
    title='LightCurve',
    x_axis='Julian Data',
    y_axis='Flux',
    label='Lightcurve'):
    
    p = figure(title=title,
            plot_width=650, plot_height=400)

    p.xaxis[0].axis_label = x_axis
    p.yaxis[0].axis_label = y_axis

    p.line(x_data, y_data, line_width=2, legend_label=label)

    show(p)

def view_fourier(
    freq_data, 
    sp_data,
    title='Original Fourier Spectrum'):

    y_min = 1
    y_max = round(np.real(max(sp_data)), 0)
    y_max = 10**(len(str(y_max))-2)

    p = figure(title=title,
            y_axis_type='log',
            plot_width=650, plot_height=400,
            background_fill_color='#fafafa',
            #x_range=(x_min, x_max),
            y_range=(y_min, y_max))

    p.line(freq_data, np.real(sp_data),
            line_color='blue',
            legend_label='Original Spectrum',
            line_width=2)

    p.xaxis[0].axis_label = 'Frequency'
    p.yaxis[0].axis_label = 'Magnitude'
    
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

  r = p.multi_line(xs, ys, color=["blue", "red"], line_width=2)

  legend = Legend(items=[
    LegendItem(label="Original", renderers=[r], index=0),
    LegendItem(label="Filtered", renderers=[r], index=1),
  ])
  p.add_layout(legend)

  show(p)

def compare_fourier_original_filtered(
  freq_original_data=None,
  sp_original_data=None,
  freq_filtered_data=None,
  sp_filtered_data=None,
  title='Original and Filtered Fourier Spectrum'):

  y_min = 1
  y_max = round(np.real(max(sp_original_data)), 0)
  y_max = 10**(len(str(y_max))-2)

  p = figure(title=title, 
        y_axis_type='log',
        plot_width=650, plot_height=400,
        background_fill_color='#fafafa',
        #x_range=(x_min, x_max),
        y_range=(y_min, y_max)
        )
  
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

def view_phase_folded_simulation(
  original_time=None,
  original_flux=None,
  simulated_time=None,
  simulated_flux=None,
  x_axis='Julian Date',
  y_axis='Flux',
  title='Phase-Folded Compared'):

  p = figure(title=title,
        plot_width=650, plot_height=400)

  p.xaxis[0].axis_label = x_axis
  p.yaxis[0].axis_label = y_axis

  time = [original_time, simulated_time]
  flux = [original_flux, simulated_flux]

  r = p.multi_line(time, flux, color=["blue", "red"], line_width=2)

  legend = Legend(items=[
    LegendItem(label="Original", renderers=[r], index=0),
    LegendItem(label="Simulated", renderers=[r], index=1),
  ])
  p.add_layout(legend)

  show(p)


def make_histogram(
  data=None,
  bins=50,
  gaussian_approx=False,
  plot_pdf=True,
  plot_cdf=True):

  hist, edges = np.histogram(data, density=True, bins=bins)

  p = figure(title='Histogram plot',
          plot_width=650, plot_height=400,
          background_fill_color='#fafafa')

  p.quad(top=hist, bottom=0, left=edges[:-1], right=edges[1:],
          fill_color='navy', line_color='white', alpha=0.5)

  p.y_range.start = 0

  if gaussian_approx:
    x = np.linspace(min(data)-1, max(data)+1, len(data))

    mu_computed, sigma_computed = norm.fit(data)

    pdf = 1/(sigma_computed * np.sqrt(2*np.pi)) * np.exp(-(x-mu_computed)**2 / (2*sigma_computed**2))

    cdf = (1+scipy.special.erf((x-mu_computed)/np.sqrt(2*sigma_computed**2)))/2

    p.title.text = f'Normal Distribution Approximation (μ≈{round(mu_computed,4)}, σ≈{round(sigma_computed,4)})'

    if plot_pdf:
      p.line(x, pdf, line_color='#ff8888', line_width=4, alpha=0.7, legend_label='PDF')
  
    if plot_cdf:
      p.line(x, cdf, line_color="orange", line_width=2, alpha=0.7, legend_label='CDF')

    p.y_range.start = 0
    p.legend.location = "center_right"
    p.xaxis.axis_label = 'x'
    p.yaxis.axis_label = 'Pr(x)'
    # p.grid.grid_line_color="white"

    show(p)
    return mu_computed, sigma_computed

  show(p)

