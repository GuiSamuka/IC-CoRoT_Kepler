from bokeh.plotting import figure, output_file, show
from bokeh.models import MultiLine, Plot
import pandas as pd

def view_lightcurve(
  x_data=None, 
  y_data=None,
  title='Sample Light Curve',
  x_axis='Date',
  y_axis='Whiteflux',
  label='Lightcurve'):
  
  x_data = pd.to_datetime(x_data)

  output_file('line.html')

  p = figure(x_axis_type="datetime", 
            title=title,
            plot_width=1200, plot_height=600)

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

  output_file('line.html')

  pass

   


