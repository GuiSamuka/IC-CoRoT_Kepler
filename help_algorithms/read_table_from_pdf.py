from os import system
system("cls")

# pip install tabula-py

import pandas as pd
import tabula
from tabula.io import read_pdf
import numpy as np

path = r"C:\Users\guisa\Google Drive\01 - Iniciação Científica\01 - Referências\Deleuil et. all (2018).pdf"

lista_table_7 = tabula.read_pdf(path, pages='22', multiple_tables=True)
lista_table_8 = tabula.read_pdf(path, pages='23', multiple_tables=True)

tabela_7 = pd.DataFrame(lista_table_7[0])  # 46 valores
eclipsing_binaries_ids = tabela_7['CoRoT-ID'].values

tabela_8 = pd.DataFrame(lista_table_8[0])  # 48 valores
eclipsing_binaries_ids = np.append(eclipsing_binaries_ids, tabela_8['CoRoT-ID'].values)
eclipsing_binaries_ids = eclipsing_binaries_ids[~np.isnan(eclipsing_binaries_ids)]

print("Eclipsing Binaries IDs:\n", eclipsing_binaries_ids)
print("Total of Eclipsing Binaries:", len(eclipsing_binaries_ids))


# np.savetxt('tests/eclipsing_binaries_ids.txt', eclipsing_binaries_ids, fmt='%9.0f', delimiter=',', newline=',')

