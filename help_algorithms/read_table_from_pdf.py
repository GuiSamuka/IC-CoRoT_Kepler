from os import system
system("cls")

# pip install tabula-py
#%%
import pandas as pd
import tabula
from tabula.io import read_pdf
import numpy as np

path = r"C:\Users\guisa\Google Drive\01 - Iniciação Científica\01 - Referências\Deleuil et. al., 2018.pdf"

lista_table_5 = tabula.read_pdf(path, pages='20', multiple_tables=True)
# lista_table_7 = tabula.read_pdf(path, pages='22', multiple_tables=True)
# lista_table_8 = tabula.read_pdf(path, pages='23', multiple_tables=True)

# tabela_7 = pd.DataFrame(lista_table_7[0])  # 46 valores
# eclipsing_binaries_ids = tabela_7['CoRoT-ID'].values

# tabela_8 = pd.DataFrame(lista_table_8[0])  # 48 valores
# eclipsing_binaries_ids = np.append(eclipsing_binaries_ids, tabela_8['CoRoT-ID'].values)
# eclipsing_binaries_ids = eclipsing_binaries_ids[~np.isnan(eclipsing_binaries_ids)]

tabela_5 = pd.DataFrame(lista_table_5[0])
tabela_5 = tabela_5.drop(0, axis=0)

#%%

tabela_5['Period'].min() # 0.520945 ± 0.00000   -> period
tabela_5['Period'].min() # 79.967743 ± 0.00297  -> period

tabela_5['Rp/R?'].min()  # 0.0271 ± 0.002 -> p
tabela_5['Rp/R?'].max()  # 0.5690 ± 0.000 -> p

tabela_5['a/R?'].min()   # 1.36 ± 119.91  -> adivR
tabela_5['a/R?'].max()   # 9.67 ± 74.42   -> adivR

tabela_5['b'].min()      # 0.00 ± 2.78    -> b_impact
tabela_5['b'].max()      # 1.50 ± 0.00    -> b_impact


"""Conclusões:

period_values = np.arange(0, 80, 0.01)   -> 8000 valores, 2 casas decimais
p_values = np.arange(0, 0.6, 0.01)       -> 60 valores, 2 casas decimais
adivR_values = np.arange(1, 10, 0.1)     -> 90 valores, 1 casa decimal
b_values = np.arange(0, 1.6, 0.01)       -> 160 valores, 2 casas decimais

500_000 -> 10min
6_912_000_000 -> x

x = 6_912_000_000 * 10 / 500_000

"""

"""CoRoT-ID que eu tenho:

Precisao:
periodo: 2 casas decimais -> 50 valores
p: 3 casas decimais       -> 6 valores
adivR: 2 casas decimais   -> 50 valores
b: 2 casas decimais       -> 5 valores
                             75.000 valores


tabela_5[tabela_5['CoRoT-ID'] == ID]

period, p e b -> 2 casas decimais
adivR -> 1 casa

100725706 
    period = [13.22, 13.23, 13.24, 13.25, 13.26]
    p      = [0.12, 0.13, 0.14, 0.15, 0.16]     
    adivR  = 
    b      = 


101086161  

101206560 

101368192

-----------
102671819
102671819

102708694
102708694

102725122
102725122
-----------

102764809
102890318
102912369
105118236
105209106
105228856
105793995
105819653
105833549
105891283
106017681
110839339
110864907
221686194
300001097
310247220
311519570
315198039
315211361
315239728
630831435
652180928
652180991


"""

#%%
"""Automatizando a macacada"""


def get_true_value(corot_id: int, parameter: str):
    return float(tabela_5[tabela_5['CoRoT-ID'] == corot_id][parameter].values[0].split('±')[0])

def define_interval_period(period: float):
    return np.arange(round(period, 2)-0.02, round(period, 2)+0.03, 0.01)

def define_interval_p(p: float):
    return np.arange(round(p, 2)-0.02, round(p, 2)+0.03, 0.01)

def define_interval_adivR(adivR: float):
    return np.arange(round(adivR, 2)-0.02, round(adivR, 2)+0.03, 0.01)
    
def define_interval_b(b: float):
    return np.arange(round(b, 2)-0.02, round(b, 2)+0.03, 0.01)
    


id = 101086161
period = get_true_value(id, 'Period')
p = get_true_value(id, 'Rp/R?')
adivR = get_true_value(id, 'a/R?')
b = get_true_value(id, 'b')

period_values = define_interval_period(period)
print(round(period, 2))
print(period_values, end='\n\n')

p_values = define_interval_p(p)
print(round(p, 2))
print(p_values, end='\n\n')

adivR_values = define_interval_adivR(adivR)
print(round(adivR, 2))
print(adivR_values, end='\n\n')

b_values = define_interval_b(b)
print(round(b, 2))
print(b_values)


#%%

# import os

# ids = []

# my_dir = r"C:\Users\guisa\Google Drive\01 - Iniciação Científica\IC-CoRoT_Kepler\resampled_files"
# for root_dir_path, sub_dirs, files in os.walk(my_dir):
#     for j in range(0, len(files)):
#         if files[j].endswith('.csv'):
#             id = files[j].split('_')[1]
#             ids.append(id)

            

# # print("Eclipsing Binaries IDs:\n", eclipsing_binaries_ids)
# # print("Total of Eclipsing Binaries:", len(eclipsing_binaries_ids))

# print(tabela_5['FU'].value_counts())

# # np.savetxt('tests/eclipsing_binaries_ids.txt', eclipsing_binaries_ids, fmt='%9.0f', delimiter=',', newline=',')


# %%
