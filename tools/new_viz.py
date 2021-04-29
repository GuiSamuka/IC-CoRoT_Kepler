"""

    This module ...

"""

# Imports
from bokeh.plotting import figure, output_file, show
from bokeh.io import output_notebook
import pandas as pd

output_notebook()

def view_lightcurve(
  x_data=None, 
  y_data=None,
  title='Sample Light Curve',
  x_axis='Date',
  y_axis='Whiteflux',
  label='Lightcurve'):
  
  x_data = pd.to_datetime(x_data)

  #output_file('line.html')

  p = figure(x_axis_type="datetime", 
            title=title,
            plot_width=650, plot_height=400)

  p.xaxis[0].axis_label = x_axis
  p.yaxis[0].axis_label = y_axis

  p.line(x_data, y_data, line_width=2, legend_label=label)

  show(p)

def view_filter_results(
  x_original=None,
  y_original=None,
  x_filtered=None,
  y_filtered=None,
  title='Filtering result',
  x_axis='Date',
  y_axis='Whiteflux'
  ):
  
  x_original = pd.to_datetime(x_original)
  x_filtered = pd.to_datetime(x_filtered)

  #output_file('line.html')

  p = figure(x_axis_type="datetime", 
            title=title,
            plot_width=650, plot_height=400)

  p.xaxis[0].axis_label = x_axis
  p.yaxis[0].axis_label = y_axis

  xs = [x_original, x_filtered]
  ys = [y_original, y_filtered]

  p.multi_line(xs, ys, color=["blue", "red"], line_width=2)

  show(p)
   
def line_plot():
  pass

