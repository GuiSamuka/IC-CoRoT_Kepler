import plotly.graph_objs as go

def view_lightcurve(
  x_data=None, 
  y_data=None,
  title='Title example',
  x_axis='x axis',
  y_axis='y axis',
  label='label'):


  fig = go.Figure()

  fig.update_layout(title=title,
                    xaxis_title=x_axis,
                    yaxis_title=y_axis)

  fig.add_trace(go.Scatter(x=x_data, y=y_data,
                           mode='lines',
                           name=label))
  fig.show()