from os import system
from types import FrameType
system("cls")

# observed_curve_path = 'files\curva_luz_eclipse_medio_ID100725706_Butterworth_n2_f02_autocalibrada.txt'
# flux = []

# vetor = [x.split('\n')[0] for x in open(observed_curve_path).readlines()]
# for line in vetor:
#     splitted_line = line.split(' ')
#     splitted_line_iterator = filter(None, splitted_line)
#     splitted_line_filtered = list(splitted_line_iterator)
#     flux.append(float(splitted_line_filtered[1]))

# print(flux)
# print(primeira_linha_splitada_filtrada)
# flux.append(primeira_linha_splitada_filtrada[1])
# print(flux)



# import numpy as np
# Nmax = 10
# new_final_table = np.zeros((5, Nmax))
# new_sorted_final_table = np.zeros((5, Nmax))

# print(new_final_table)
# print()
# print(new_final_table[4,:])


# # %%
# import numpy as np
# from scipy.stats import norm
# import matplotlib.pyplot as plt
# import pandas as pd

# dataset = pd.read_csv('final_table.csv')
# print(dataset.columns)
# data = dataset['qui2']

# mean,std=norm.fit(data)


# plt.hist(data, bins=30, density=True, histtype='stepfilled', alpha=0.2)
# xmin, xmax = plt.xlim()
# x = np.linspace(xmin, xmax, 100)
# y = norm.pdf(x, mean, std)
# plt.plot(x, y)
# plt.title(f'{data.name} histogram')
# plt.show()



# # %%
# import numpy as np
# import pandas as pd

# dataset = pd.read_csv('final_table.csv')

# sorted = dataset.sort_values(by='qui2')

# b = sorted.loc[sorted.index[0]][0]
# p = sorted.loc[sorted.index[0]][1]
# period = sorted.loc[sorted.index[0]][2]
# adivR = sorted.loc[sorted.index[0]][3]
# qui2 = sorted.loc[sorted.index[0]][4]

# print(sorted.loc[sorted.index[0]])
# print(qui2)




################################
## Testando classe
#%%

import utils
import numpy as np


# Planet coordinate, along the x-axis, as a function of the start's radius
x_values_path = 'files\Valores_x_simulacao.txt'
x_values = np.loadtxt(x_values_path, dtype='float', delimiter='\n')

# Transit impact parameter
b_values_path = 'files\Valores_b_simulacao_rodada2.txt'
b_values = np.loadtxt(b_values_path, dtype='float', delimiter='\n')

# Radius values of the planet compared to the star
p_values_path = 'files\Valores_p_simulacao_rodada2.txt'
p_values = np.loadtxt(p_values_path, dtype='float', delimiter='\n')

# Orbital period values to be considered
period_values_path = 'files\Valores_periodo_simulacao_rodada2.txt'
period_values = np.loadtxt(period_values_path, dtype='float', delimiter='\n')

# Orbital radius values compared to star radius
adivR_values_path = 'files\Valores_adivR_simulacao_rodada2.txt'
adivR_values = np.loadtxt(adivR_values_path, dtype='float', delimiter='\n')


observed_curve_path = 'files\curva_luz_eclipse_medio_ID100725706_Butterworth_n2_f02_autocalibrada.txt'
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


simulation_object = utils.Simulate()
final_table_sorted_by_chi2 = simulation_object.simulate_values(observed_curve, [0.91, 0.92, 0.12, 0.13, 0.123, 0.123], [0.13, 0.14], [13, 14], [18.3, 18.4], x_values, results_to_csv=False)
print(final_table_sorted_by_chi2.head())
# final_table_sorted_by_chi2 = simulation_object.simulate_values(observed_curve, b_values, p_values, period_values, adivR_values, x_values, results_to_csv=True)
# print(final_table_sorted_by_chi2.head())

simulated_curve = simulation_object.simulate_lightcurve(observed_curve=observed_curve, x_values=x_values)
simulated_curve.view_simulation_results()
simulated_curve.compare_results(see_values=True)


# %%
## Calculo do erro por parametro
## Capítulo 14.1 do livro: Numerical recipe : the art of scientific computing

# https://docs.bokeh.org/en/latest/docs/gallery/histogram.html

import utils
import numpy as np
from math import sqrt

adivR_values_path = 'files\Valores_adivR_simulacao_rodada2.txt'
data = np.loadtxt(adivR_values_path, dtype='float', delimiter='\n')

# mu, sigma = 0, 0.5
# data = np.random.normal(mu, sigma, 1000)

# utils.make_histogram(data, bins=int(round(sqrt(len(data)), 0)))
# mu, sigma = utils.make_histogram(data, bins=int(round(sqrt(len(data)), 0)), gaussian_approx=True, plot_cdf=False)
# print(f'Computed values: μ≈{round(mu,4)} and σ≈{round(sigma,4)}')



import numpy as np
import statsmodels.stats.api as sms

sms.DescrStatsW(data).tconfint_mean()


#%%
from uncertainties import ufloat
import numpy as np

array = [1, 2]
print('Original array =', array)

std = np.std(array)

# uncert_arr = np.array([ufloat(array[0], std), ufloat(array[1], std)])
uncert_arr = []
for i in range(len(array)):
    uncert_arr.append(ufloat(array[i], std))
uncert_arr = np.array(uncert_arr)


print('Array with uncertainties =', uncert_arr)


# %%
import utils
import pandas as pd
import numpy as np
np.set_printoptions(suppress=True)

# path = 'files\Valores_b_simulacao_rodada2.txt'
# data = np.loadtxt(path, dtype='float', delimiter='\n')

table = pd.read_csv('final_table.csv')
# Index(['b_impact', 'p', 'period', 'adivR', 'chi2'], dtype='object')
data = table['p'].to_numpy()


# (unique, counts) = np.unique(data, return_counts=True)
# frequencies = np.asarray((unique, counts)).T
# print(frequencies)


mu, sigma = utils.make_histogram(data, bins=10, gaussian_approx=True, plot_cdf=False)
print(f'σ = {round(sigma, 4)} => Confidence = 68.67%')
print(f'2σ = {round(2*sigma, 4)} => Confidence = 95.45%')
print(f'3σ = {round(3*sigma, 4)} => Confidence = 99.73%')
print(f'1.645σ = {round(1.645*sigma, 4)} => Confidence = 90.00%')
print(f'2.576σ = {round(2.576*sigma, 4)} => Confidence = 99.00%')
print(f'∆=0.6745σ = {round(0.6745*sigma, 4)} => Confidence = 50.00%')


# mu, sigma = utils.make_histogram(data, bins=int(round(sqrt(len(data)))), gaussian_approx=True, plot_cdf=False)




# %%

