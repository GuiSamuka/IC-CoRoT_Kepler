import filters
import viz

import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
from os import system
system("cls")


FILE_PATH = 'C:/Users/guisa/Google Drive/01 - Iniciação Científica/02 - Datasets/csv_files/EN2_STAR_CHR_0101086161_20070516T060226_20071005T074409.csv'
data_sample = pd.read_csv(FILE_PATH)
x = data_sample.DATE.to_numpy()
y = data_sample.WHITEFLUX.to_numpy()


import new_filters
import new_viz

Filter = new_filters.FrequencyDomainFiltering()
Filter.filter(array=y, filter='butterworth', numExpansion=70, cutoff_freq=0.05, order=2)
filtered = Filter.getFiltered
#new_viz.view_lightcurve(x, y)
new_viz.view_filter_results(x, y, x, filtered)



# Filter.expand_borders(y, 70)
# expanded = Filter.getExpanded
# viz.view_lightcurve(x, expanded)

#Filter.filter(array=y, filter='butterworth', numExpansion=70, cutoff_freq=0.05, order=2)
#Filter.filter(array=y, filter='gaussian', numExpansion=70, cutoff_freq=0.25, order=None)

#y_filtered = Filter.getFiltered
#viz.view_filter_results(x, y, x, y_filtered)




"""
# Median Filter
vetor = np.array([2, 3, 80, 6, 2])

# Algoritmo
filtrado = np.zeros(len(vetor))

primeiro_pixel = vetor[0]
ultimo_pixel = vetor[-1]

num_neibo = 4   # num_neibo no maximo tem que ser a divisão inteira de 

output = []
for i in range(len(vetor) // num_neibo + len(vetor) % num_neibo):
    a = np.median(vetor[i:num_neibo+i])
    i += 1
    # print(a)  
    output.append(a)  
    
output = np.array(output)

filtered = np.concatenate(([primeiro_pixel], output, [ultimo_pixel]))

# print(filtered)

# 5/5
# 5/4
# 5/3
# 5/2
# 5/1

"""