#%%
# Imports
from re import S
from models.lightcurve import LightCurve
import pandas as pd

LIGHTCURVE = 'RESAMPLED_0102912369_20070203T130553'

data = pd.read_csv('https://raw.githubusercontent.com/Guilherme-SSB/IC-CoRoT_Kepler/main/resampled_files/' + LIGHTCURVE + '.csv')

time = data.DATE.to_numpy()
flux = data.WHITEFLUX.to_numpy()

curve = LightCurve(time, flux)
# curve.plot() # Works
# curve.view_fourier_spectrum() # Works

###### Filtering tests
# curve.ideal_lowpass_filter(0.2).view_filter_results() # Works
# curve.ideal_lowpass_filter(0.2).view_fourier_results() # Works

# curve.gaussian_lowpass_filter(0.2).view_filter_results() # Works
# curve.gaussian_lowpass_filter(0.2).view_fourier_results() # Works

# curve.butterworth_lowpass_filter(6, 0.4).view_filter_results() # Works
# curve.butterworth_lowpass_filter(6, 0.4).view_fourier_results() # Works

# curve.bessel_lowpass_filter(2, 0.4).view_filter_results() # Works
# curve.bessel_lowpass_filter(2, 0.4).view_fourier_results() # Works
 
# curve.median_filter(5).view_filter_results() # Works
# curve.median_filter(5).view_fourier_results() # Works

# %%
# Imports
from models.simulation import Simulate

import numpy as np

###### Loading parameters to simulate
path = 'C:/Users/guisa/Google Drive/01 - Iniciação Científica/IC-CoRoT_Kepler/Light Curve Simulation/'

# Planet coordinate, along the x-axis, as a function of the start's radius
x_values_path = path+'files\Valores_x_simulacao.txt'
x_values = np.loadtxt(x_values_path, dtype='float', delimiter='\n')

# Lightcurve phase-folded
observed_curve_path = path+'files\curva_luz_eclipse_medio_ID100725706_Butterworth_n2_f02_autocalibrada.txt'
vetor = [x.split('\n')[0] for x in open(observed_curve_path).readlines()]

Nobs = 0
for line in vetor:
    Nobs += 1

observed_curve = np.zeros((Nobs, 3))

time = []
flux = []
flux_error = []

for line in vetor:
    splitted_line = line.split(' ')
    splitted_line_iterator = filter(None, splitted_line)
    splitted_line_filtered = list(splitted_line_iterator)
    time.append(splitted_line_filtered[0])
    flux.append(splitted_line_filtered[1])
    flux_error.append(splitted_line_filtered[2])

observed_curve[:, 0] = time
observed_curve[:, 1] = flux
observed_curve[:, 2] = flux_error

# Creating a Simulate object
SimulateObject = Simulate()

## Simulate a lightcurve with given parameters
simulated_curve = SimulateObject.simulate_lightcurve(observed_curve=observed_curve, b_impact=0.92, p=0.1483, period=13.240160, adivR=22.36, x_values=x_values)
simulated_curve.view_simulation_results()
simulated_curve.compare_results(see_values=True)


# %%
# Imports
from models.simulation import Simulate

import numpy as np

###### Loading parameters to simulate
path = 'C:/Users/guisa/Google Drive/01 - Iniciação Científica/IC-CoRoT_Kepler/Light Curve Simulation/'

# Planet coordinate, along the x-axis, as a function of the start's radius
x_values_path = path+'files\Valores_x_simulacao.txt'
x_values = np.loadtxt(x_values_path, dtype='float', delimiter='\n')

# Transit impact parameter
b_values_path = path+'files\Valores_b_simulacao_rodada2.txt'
b_values = np.loadtxt(b_values_path, dtype='float', delimiter='\n')

# Radius values of the planet compared to the star
p_values_path = path+'files\Valores_p_simulacao_rodada2.txt'
p_values = np.loadtxt(p_values_path, dtype='float', delimiter='\n')

# Orbital period values to be considered
period_values_path = path+'files\Valores_periodo_simulacao_rodada2.txt'
period_values = np.loadtxt(period_values_path, dtype='float', delimiter='\n')

# Orbital radius values compared to star radius
adivR_values_path = path+'files\Valores_adivR_simulacao_rodada2.txt'
adivR_values = np.loadtxt(adivR_values_path, dtype='float', delimiter='\n')

# Lightcurve phase-folded
observed_curve_path = path+'files\curva_luz_eclipse_medio_ID100725706_Butterworth_n2_f02_autocalibrada.txt'
vetor = [x.split('\n')[0] for x in open(observed_curve_path).readlines()]

Nobs = 0
for line in vetor:
    Nobs += 1

observed_curve = np.zeros((Nobs, 3))

time = []
flux = []
flux_error = []

for line in vetor:
    splitted_line = line.split(' ')
    splitted_line_iterator = filter(None, splitted_line)
    splitted_line_filtered = list(splitted_line_iterator)
    time.append(splitted_line_filtered[0])
    flux.append(splitted_line_filtered[1])
    flux_error.append(splitted_line_filtered[2])

observed_curve[:, 0] = time
observed_curve[:, 1] = flux
observed_curve[:, 2] = flux_error


# Creating a Simulate object
SimulateObject = Simulate()

## Simulate values
final_table_sorted_by_chi2 = SimulateObject.simulate_values(observed_curve, [0.91, 0.92], [0.13, 0.14], [13, 14], [18.3, 18.4], x_values, results_to_csv=False)
print(final_table_sorted_by_chi2.head())

## Build the lightcurve with the best parameters computeds
simulated_curve = SimulateObject.simulate_lightcurve(observed_curve=observed_curve, x_values=x_values)
simulated_curve.view_simulation_results()
simulated_curve.compare_results(see_values=True)

# %%
