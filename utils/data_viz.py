"""

    This module ...

"""


# Imports
from bokeh.plotting import figure, output_file, show
from bokeh.io import output_notebook
import pandas as pd

output_file('image.html')
# output_notebook()


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