import plotly.graph_objs as go

def view_lightcurve(
  x_data=None, 
  y_data=None,
  title='Sample Light Curve',
  x_axis='Date',
  y_axis='Whiteflux',
  label='label'):

  fig = go.Figure()

  fig.update_layout(title=title,
                    xaxis_title=x_axis,
                    yaxis_title=y_axis)

  fig.add_trace(go.Scatter(x=x_data, y=y_data,
                           mode='lines',
                           name=label))
  fig.show()

def view_filter_results(
  x_original=None,
  y_original=None,
  x_filtered=None,
  y_filtered=None,
  title='Filtering result',
  x_axis='Date',
  y_axis='Whiteflux'
  ):

  fig = go.Figure()

  fig.update_layout(title=title,
                    xaxis_title=x_axis,
                    yaxis_title=y_axis)

  fig.add_trace(go.Scatter(x=x_original, y=y_original,
                          mode='lines',
                          name='Original Light Curve'))

  fig.add_trace(go.Scatter(x=x_filtered, y=y_filtered,
                          mode='lines',
                          name='Light Curve Filtered'))
  fig.show()