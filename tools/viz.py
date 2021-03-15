import plotly.graph_objs as go

def view_lightcurve(x, y, name='title'):
    fig = go.Figure()

    fig.add_trace(go.Scatter(x=x, y=y,
                            mode='lines',
                            name=name))

    fig.update_layout(title='Lightcurve',
                    xaxis_title='Date',
                    yaxis_title='Whiteflux')

    fig.show()
