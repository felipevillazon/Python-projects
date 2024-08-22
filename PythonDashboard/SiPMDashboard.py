## Dashboard for SiPM Control and QSL database connection and reading

import dash
from dash import dcc, html
import plotly.graph_objs as go
import numpy as np
import random
from dash.dependencies import Input, Output
import plotly.express as px
import dash_bootstrap_components as dbc

# Initialize the Dash app
external_stylesheets = [dbc.themes.CERULEAN]
app = dash.Dash(__name__, external_stylesheets=external_stylesheets)

# Define color for the text
colors = {
    'background': '#DADADA',
    'text': '#ffc600'
}

ColorAlarmLL = '#F9FBFF' 
ColorAlarmL = '#F9FBFF' 
ColorAlarmH = '#F9FBFF' 
ColorAlarmHH = '#F9FBFF' 



# Define the layout of the dashboard
app.layout = html.Div(style={
    'height': '100vh',        # Full viewport height
    'margin': '0',           # Remove default margin
    'display': 'flex',       # Enable Flexbox for layout
    'flex-direction': 'column',  # Arrange children in a column (top to bottom)
    'backgroundColor': colors['background'],
}, children=[
    
    # Title Section
    html.Div(style={
        'textAlign': 'center',    # Center the text horizontally
        'padding': '10px',        # Add padding around the title
        'backgroundColor': '#0091D5', # Optional background color for the title area
    }, children=[
        html.H1(
            children='SiPM Temperature Control Data',
            style={
                'textAlign': 'center',    # Center the text horizontally within the div
                'color': '#F9FBFF'   # Set the text color
            }
        )
    ]),

    # Graphs Section
    html.Div(style={
        'height': 'calc(100vh - 60px)',  # Height adjusted to account for the title section
        'display': 'flex',              # Enable Flexbox for layout within this section
        'flex-direction': 'row',        # Arrange children in a row (left to right)
        'flex-wrap': 'wrap',             # Allow wrapping to the next line if needed
    }, children=[
        
        # Top-Left Graph
        html.Div(style={
            'width': '50%',    # Each graph takes up half of the container's width
            'height': '50%',   # Each graph takes up half of the container's height
            'padding': '5px', # Add padding around the graph
            'box-sizing': 'border-box'  # Include padding in the width and height calculation
        }, children=[
            dcc.Graph(
                id='top-left-graph',
                style={'height': '100%'}  # Graph stretches to fill the container's height
            )
        ]),

        # Top-Right Section with 9 small boxes
        html.Div(style={
            'width': '50%',         # Take up half of the width of the parent container
            'height': '50%',        # Take up half of the height of the parent container
            'padding': '5px',      # Add padding around the entire section
            'box-sizing': 'border-box',  # Include padding in the width and height calculation
            'display': 'flex',      # Enable Flexbox layout
            'flex-direction': 'column',  # Arrange child elements in a column
            'flex-wrap': 'wrap'     # Allow wrapping to the next line if needed
        }, children=[
            
            # Top row with three larger boxes
            html.Div(style={
                'display': 'flex',        # Arrange child elements in a row
                'flex-direction': 'row',  # Arrange child elements in a row
                'flex-wrap': 'wrap',      # Allow wrapping
                'width': '100%',          # Take up full width of the parent container
                'height': '50%',          # Take up half of the height of the section
                'padding': '5px',         # Add padding around the boxes
                'box-sizing': 'border-box' # Include padding in the width and height calculation
               
            }, children=[
                html.Div(id='top-left-box', style={
                    'flex': '1',            # Each box takes up equal space
                    'margin': '5px',        # Add margin between boxes
                    'background-color': '#F9FBFF',  # Background color for better visibility
                    'height': '100%',        # Full height of the row
                    'display': 'flex',
                    'align-items': 'center',
                    'justify-content': 'center'
                }),

                html.Div(id='top-center-box', style={
                    'flex': '1',            
                    'margin': '5px',        
                    'background-color': '#F9FBFF',
                    'height': '100%',
                    'display': 'flex',
                    'align-items': 'center',
                    'justify-content': 'center'
                }),

                html.Div(id='top-right-box', style={
                    'flex': '1',            
                    'margin': '5px',        
                    'background-color': '#F9FBFF',
                    'height': '100%',
                    'display': 'flex',
                    'align-items': 'center',
                    'justify-content': 'center'
                })
            ]),

            # Middle row with two medium boxes
            html.Div(style={
                'display': 'flex',        # Arrange child elements in a row
                'flex-direction': 'row',  # Arrange child elements in a row
                'flex-wrap': 'wrap',      # Allow wrapping
                'width': '100%',          # Take up full width of the parent container
                'height': '25%',          # Take up a quarter of the height of the section
                'padding': '5px',         # Add padding around the boxes
                'box-sizing': 'border-box'  # Include padding in the width and height calculation
            }, children=[
                html.Div(id='middle-left-box', style={
                    'flex': '1',            # Each box takes up equal space
                    'margin': '5px',        # Add margin between boxes
                    'background-color': '#F9FBFF',  # Background color for better visibility
                    'height': '100%',        # Full height of the row
                    'display': 'flex',
                    'align-items': 'center',
                    'justify-content': 'center'
                }),

                html.Div(id='middle-right-box', style={
                    'flex': '1',            
                    'margin': '5px',        
                    'background-color': '#F9FBFF',
                    'height': '100%',
                    'display': 'flex',
                    'align-items': 'center',
                    'justify-content': 'center'
                })
            ]),

            # Bottom row with four smaller boxes
            html.Div(style={
                'display': 'flex',        # Arrange child elements in a row
                'flex-direction': 'row',  # Arrange child elements in a row
                'flex-wrap': 'wrap',      # Allow wrapping
                'width': '100%',          # Take up full width of the parent container
                'height': '25%',          # Take up a quarter of the height of the section
                'padding': '5px',         # Add padding around the boxes
                'box-sizing': 'border-box'  # Include padding in the width and height calculation
            }, children=[
                html.Div(id='bottom-left-box-1', style={
                    'flex': '1',            # Each box takes up equal space
                    'margin': '5px',        # Add margin between boxes
                    'background-color': ColorAlarmLL,  # Background color for better visibility
                    'height': '100%',        # Full height of the row
                    'display': 'flex',
                    'align-items': 'center',
                    'justify-content': 'center'
                }),

                html.Div(id='bottom-left-box-2', style={
                    'flex': '1',            
                    'margin': '5px',        
                    'background-color': ColorAlarmL,
                    'height': '100%',
                    'display': 'flex',
                    'align-items': 'center',
                    'justify-content': 'center'
                }),

                html.Div(id='bottom-right-box-1', style={
                    'flex': '1',            
                    'margin': '5px',        
                    'background-color':ColorAlarmH,
                    'height': '100%',
                    'display': 'flex',
                    'align-items': 'center',
                    'justify-content': 'center'
                }),

                html.Div(id='bottom-right-box-2', style={
                    'flex': '1',            
                    'margin': '5px',        
                    'background-color':ColorAlarmHH,
                    'height': '100%',
                    'display': 'flex',
                    'align-items': 'center',
                    'justify-content': 'center'
                })
            ])
        ]),

        # Bottom-Left Graph
        html.Div(style={
            'width': '50%',
            'height': '50%',
            'padding': '5px',
            'box-sizing': 'border-box'
        }, children=[
            dcc.Graph(
                id='bottom-left-graph',
                style={'height': '100%', 'width': '100%'}
            )
        ]),

        # Bottom-Right Graph
        html.Div(style={
            'width': '50%',
            'height': '50%',
            'padding': '5px',
            'box-sizing': 'border-box'
        }, children=[
            dcc.Graph(
                id='bottom-right-graph',
                style={'height': '100%'}
            )
        ])
    ])
])

# Callback to update random data for the plots and boxes
@app.callback(
    [
        Output('top-left-graph', 'figure'),
        Output('bottom-left-graph', 'figure'),
        Output('bottom-right-graph', 'figure'),
        Output('top-left-box', 'children'),
        Output('top-center-box', 'children'),
        Output('top-right-box', 'children'),
        Output('middle-left-box', 'children'),
        Output('middle-right-box', 'children'),
        Output('bottom-left-box-1', 'children'),
        Output('bottom-left-box-2', 'children'),
        Output('bottom-right-box-1', 'children'),
        Output('bottom-right-box-2', 'children')
    ],
    Input('interval-component', 'n_intervals')
)
def update_data(n_intervals):
    # Generate random data
    time = np.linspace(0, 50, 10)
    
     # Generate temperature data (ln1, ln2, ln3) within the range 230 to 260
    ln1 = np.random.uniform(230, 260, 10)
    ln2 = np.random.uniform(230, 260, 10)
    ln3 = np.random.uniform(230, 260, 10)

    # Generate flux data (fi, fo) within the range 0 to 24
    fi = np.random.uniform(0, 24, 10)
    fo = np.random.uniform(0, 24, 10)
    Alarm_LL = random.choice([True, False])
    Alarm_L = random.choice([True, False])
    Alarm_H = random.choice([True, False])
    Alarm_HH = random.choice([True, False])
    
    # Create figures for the graphs using plotly.express

    # Create figures for the graphs
    top_left_fig = {
    'data': [
        go.Scatter(
            x=time,
            y=ln1,
            mode='lines',
            name='LN1',
            line=dict(color='blue', width=2, dash='solid'),
            yaxis='y1'
        ),
        go.Scatter(
            x=time,
            y=ln2,
            mode='lines',
            name='LN2',
            line=dict(color='green', width=2, dash='dash'),
            yaxis='y1'
        ),
        go.Scatter(
            x=time,
            y=ln3,
            mode='lines',
            name='LN3',
            line=dict(color='orange', width=2, dash='dot'),
            yaxis='y1'
        ),
        go.Scatter(
            x=time,
            y=fi,
            mode='lines',
            name='Flux In',
            line=dict(color='red', width=2, dash='dashdot'),
            yaxis='y2'
        ),
        go.Scatter(
            x=time,
            y=fo,
            mode='lines',
            name='Flux Out',
            line=dict(color='purple', width=2, dash='solid'),
            yaxis='y2'
        )
    ],
    'layout': go.Layout(
        title='Temperature and Flux Over Time',
        title_font=dict(size=20, color='black', family='Arial'),
        xaxis=dict(
            title='Time (s)',
            titlefont=dict(size=16, family='Arial'),
            tickfont=dict(size=14, family='Arial'),
            showline=True,
            linecolor='black',
            linewidth=2,
            fixedrange=False,  # Prevent the x-axis from moving
            rangemode='normal',  # Ensure the range doesn't auto-adjust unnecessarily
            range=[0, 50],  # Set the x-axis range
            type='linear'  # Set the scale type (linear, log, date, category)
        ),
        yaxis=dict(
            title='Temperature (°K)',
            titlefont=dict(size=16, family='Arial'),
            tickfont=dict(size=14, family='Arial'),
            showline=True,
            linecolor='black',
            linewidth=2,
            fixedrange=True,  # Prevent the y-axis from moving
            zeroline=True,
            zerolinecolor='gray',
            zerolinewidth=1.5,
            range=[220, 270],  # Set the y-axis range for temperature
            type='linear'  # Set the scale type (linear, log, date, category)
        ),
        yaxis2=dict(
            title='Flux (V)',
            titlefont=dict(size=16, family='Arial'),
            tickfont=dict(size=14, family='Arial'),
            overlaying='y',
            side='right',
            showline=True,
            linecolor='black',
            linewidth=2,
            fixedrange=True,  # Prevent the second y-axis from moving
            zeroline=True,
            zerolinecolor='gray',
            zerolinewidth=1.5,
            range=[0, 24],  # Set the y-axis range for flux
            type='linear'  # Set the scale type (linear, log, date, category)
        ),
        uirevision='constant',  # Preserve the state between updates
        legend=dict(
            x=0.5,
            y=1.1,
            xanchor='center',
            yanchor='top',
            orientation='h',
            bgcolor='rgba(255, 255, 255, 0.7)',
            bordercolor='black',
            borderwidth=1
        ),
        margin=dict(l=70, r=70, t=70, b=70),
    )
}



    bottom_left_fig = {
    'data': [
        go.Scatter(
            x=time,
            y=  abs((fi - fo))/fi*100,
            mode='lines',
            name='LN1',
            line=dict(color='black', width=2, dash='solid'),
            yaxis='y1'
        )],
    'layout': go.Layout(
        title='Flux Input/Output Difference',
        title_font=dict(size=20, color='black', family='Arial'),
        xaxis=dict(
            title='Time (s)',
            titlefont=dict(size=16, family='Arial'),
            tickfont=dict(size=14, family='Arial'),
            showline=True,
            linecolor='black',
            linewidth=2,
            fixedrange=False,  # Prevent the x-axis from moving
            rangemode='normal',  # Ensure the range doesn't auto-adjust unnecessarily
            range=[0, 50],  # Set the x-axis range
            type='linear'  # Set the scale type (linear, log, date, category)
        ),
        yaxis=dict(
            title='Flux Error (%)',
            titlefont=dict(size=16, family='Arial'),
            tickfont=dict(size=14, family='Arial'),
            showline=True,
            linecolor='black',
            linewidth=2,
            fixedrange=True,  # Prevent the y-axis from moving
            zeroline=True,
            zerolinecolor='gray',
            zerolinewidth=1.5,
            range=[0, 100],  # Set the y-axis range for temperature
            type='linear'  # Set the scale type (linear, log, date, category)
        ),
        yaxis2=dict(
            title='Flux (V)',
            titlefont=dict(size=16, family='Arial'),
            tickfont=dict(size=14, family='Arial'),
            overlaying='y',
            side='right',
            showline=True,
            linecolor='black',
            linewidth=2,
            fixedrange=True,  # Prevent the second y-axis from moving
            zeroline=True,
            zerolinecolor='gray',
            zerolinewidth=1.5,
            range=[0, 24],  # Set the y-axis range for flux
            type='linear'  # Set the scale type (linear, log, date, category)
        ),
        uirevision='constant',  # Preserve the state between updates
        legend=dict(
            x=0.5,
            y=1.1,
            xanchor='center',
            yanchor='top',
            orientation='h',
            bgcolor='rgba(255, 255, 255, 0.7)',
            bordercolor='black',
            borderwidth=1
        ),
        margin=dict(l=70, r=70, t=70, b=70),
    )
}

    bottom_right_fig = {
       'data': [
            go.Scatter(
                x=fi,
                y=ln1,
                mode= 'markers',
                name='LN1',
                 marker=dict( size=12, opacity=0.3), # Transparency level,,
                yaxis='y1'
            ),
        go.Scatter(
            x=fi,
            y=ln2,
            mode= 'markers',
            name='LN2',
             marker=dict( size=12, opacity=0.3), # Transparency level,,
            yaxis='y1'
        ),
        go.Scatter(
            x=fi,
            y=ln3,
            mode= 'markers',
            name='LN3',
             marker=dict( size=12, opacity=0.3), # Transparency level,,
            yaxis='y1'
        )
    ],
    'layout': go.Layout(
        title='Temperature vs Input Flow',
        title_font=dict(size=20, color='black', family='Arial'),
        xaxis=dict(
            title='Flow (V)',
            titlefont=dict(size=16, family='Arial'),
            tickfont=dict(size=14, family='Arial'),
            showline=True,
            linecolor='black',
            linewidth=2,
            fixedrange=False,  # Prevent the x-axis from moving
            rangemode='normal',  # Ensure the range doesn't auto-adjust unnecessarily
            range=[0, 24],  # Set the x-axis range
            type='linear'  # Set the scale type (linear, log, date, category)
        ),
        yaxis=dict(
            title='Temperature (°K)',
            titlefont=dict(size=16, family='Arial'),
            tickfont=dict(size=14, family='Arial'),
            showline=True,
            linecolor='black',
            linewidth=2,
            fixedrange=True,  # Prevent the y-axis from moving
            zeroline=True,
            zerolinecolor='gray',
            zerolinewidth=1.5,
            range=[220, 270],  # Set the y-axis range for temperature
            type='linear'  # Set the scale type (linear, log, date, category)
        ),
        yaxis2=dict(
            title='Flux (V)',
            titlefont=dict(size=16, family='Arial'),
            tickfont=dict(size=14, family='Arial'),
            overlaying='y',
            side='right',
            showline=True,
            linecolor='black',
            linewidth=2,
            fixedrange=True,  # Prevent the second y-axis from moving
            zeroline=True,
            zerolinecolor='gray',
            zerolinewidth=1.5,
            range=[0, 24],  # Set the y-axis range for flux
            type='linear'  # Set the scale type (linear, log, date, category)
        ),
        uirevision='constant',  # Preserve the state between updates
        legend=dict(
            x=0.5,
            y=1.1,
            xanchor='center',
            yanchor='top',
            orientation='h',
            bgcolor='rgba(255, 255, 255, 0.7)',
            bordercolor='black',
            borderwidth=1
        ),
        margin=dict(l=70, r=70, t=70, b=70),
    )
}

    #----- SMALL BOX DECOR ------#

    LN_fontsize = '40px'
    LN_numbersize = '40px'

    F_fontsize = '30px'
    F_numbersize = '30px'

    A_fontsize = '20px'
    A_numbersize = '20px'

    # Random values for the boxes
    box_values = {
        'top-left-box': html.Div([
            html.Span('LN1', style={'display': 'block', 'text-align': 'center', 'font-size': LN_fontsize}),
            html.Span(f'{random.random():.2f} K', style={'font-size': LN_numbersize, 'text-align': 'center'})
        ], style={'display': 'flex', 'flex-direction': 'column', 'align-items': 'center', 'justify-content': 'center', 'height': '100%'}),
        
        'top-center-box': html.Div([
            html.Span('LN2', style={'display': 'block', 'text-align': 'center', 'font-size': LN_fontsize}),
            html.Span(f'{random.random():.2f} K', style={'font-size': LN_numbersize, 'text-align': 'center'})
        ], style={'display': 'flex', 'flex-direction': 'column', 'align-items': 'center', 'justify-content': 'center', 'height': '100%'}),
        
        'top-right-box': html.Div([
            html.Span('LN3', style={'display': 'block', 'text-align': 'center', 'font-size': LN_fontsize}),
            html.Span(f'{random.random():.2f} K', style={'font-size': LN_numbersize, 'text-align': 'center'})
        ], style={'display': 'flex', 'flex-direction': 'column', 'align-items': 'center', 'justify-content': 'center', 'height': '100%'}),
        
        'middle-left-box': html.Div([
            html.Span('Flow In', style={'display': 'block', 'text-align': 'center', 'font-size': F_fontsize}),
            html.Span(f'{random.random():.2f} V', style={'font-size': F_numbersize , 'text-align': 'center'})
        ], style={'display': 'flex', 'flex-direction': 'column', 'align-items': 'center', 'justify-content': 'center', 'height': '100%'}),
        
        'middle-right-box': html.Div([
            html.Span('Flow Out', style={'display': 'block', 'text-align': 'center', 'font-size': F_fontsize}),
            html.Span(f'{random.random():.2f} V', style={'font-size': F_numbersize , 'text-align': 'center'})
        ], style={'display': 'flex', 'flex-direction': 'column', 'align-items': 'center', 'justify-content': 'center', 'height': '100%'}),
        
        'bottom-left-box-1': html.Div([
            html.Span('Alarm LL', style={'display': 'block', 'text-align': 'center', 'font-size': A_fontsize}),
        ], style={'display': 'flex', 'flex-direction': 'column', 'align-items': 'center', 'justify-content': 'center', 'height': '100%',
                 'width': '100%', 'background-color': 'red' if Alarm_LL else 'green', 'box-sizing': 'border-box'}),
        
        'bottom-left-box-2': html.Div([
            html.Span('Alarm L', style={'display': 'block', 'text-align': 'center', 'font-size': A_fontsize}),
        ], style={'display': 'flex', 'flex-direction': 'column', 'align-items': 'center', 'justify-content': 'center', 'height': '100%',
                 'width': '100%', 'background-color': 'red' if Alarm_L else 'green', 'box-sizing': 'border-box'}),
        
        'bottom-right-box-1': html.Div([
            html.Span('Alarm H', style={'display': 'block', 'text-align': 'center', 'font-size': A_fontsize}),
        ], style={'display': 'flex', 'flex-direction': 'column', 'align-items': 'center', 'justify-content': 'center', 'height': '100%',
                 'width': '100%', 'background-color': 'red' if Alarm_H else 'green', 'box-sizing': 'border-box'}),
        
        'bottom-right-box-2': html.Div([
            html.Span('Alarm HH', style={'display': 'block', 'text-align': 'center', 'font-size': A_fontsize}),
        ], style={'display': 'flex', 'flex-direction': 'column', 'align-items': 'center', 'justify-content': 'center', 'height': '100%',
                  'width': '100%', 'background-color': 'red' if Alarm_HH else 'green', 'box-sizing': 'border-box'})
    }
    
    return (
        top_left_fig,
        bottom_left_fig,
        bottom_right_fig,
        box_values['top-left-box'],
        box_values['top-center-box'],
        box_values['top-right-box'],
        box_values['middle-left-box'],
        box_values['middle-right-box'],
        box_values['bottom-left-box-1'],
        box_values['bottom-left-box-2'],
        box_values['bottom-right-box-1'],
        box_values['bottom-right-box-2']
    )

# Interval component to trigger updates
app.layout.children.append(
    dcc.Interval(
        id='interval-component',
        interval=1*1000,  # Update every 1 second
        n_intervals=0
    )
)

if __name__ == '__main__':
    app.run_server(debug=True)
