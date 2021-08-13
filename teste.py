from os import system
system("cls")

from utils import *

PATH_DIR = 'C:/Users/guisa/Desktop/filters/ideal/cutoffFreq0'

cutoff_freqs = [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9]
for cutoff_freq in cutoff_freqs:
    data_helper.export_results_csv(PATH_DIR+str(int(cutoff_freq*10)), 'ideal', cutoff_freq=cutoff_freq, order=None, numNei=None)
