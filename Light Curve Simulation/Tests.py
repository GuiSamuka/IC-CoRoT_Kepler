observed_curve_path = 'files\curva_luz_eclipse_medio_ID100725706_Butterworth_n2_f02_autocalibrada.txt'

vetor = [x.split('\n')[0] for x in open(observed_curve_path).readlines()]
# vetor contém toda a informação do arquivo

primeira_linha = vetor[0]
primeira_linha_splitada = primeira_linha.split(' ')
primeira_linha_splitada_filtrada_iterator = filter(None, primeira_linha_splitada)
primeira_linha_splitada_filtrada = list(primeira_linha_splitada_filtrada_iterator)

print(primeira_linha_splitada_filtrada)
