"""
Started in: August 04, 2021
Finished in: 

@author: Guilherme Samuel
"""
from os import system

from scipy.signal.signaltools import detrend
from scipy.signal.spectral import periodogram
system("cls")

RESAMPLED_ECLIPSING_BINARIES_FILES_FOLDER_PATH = 'C:/Users/guisa/Desktop/resampled_eclipsing_binaries'
RESAMPLED_CONFIRMED_EXOPLANETS_FILES_FOLDER_PATH = 'C:/Users/guisa/Google Drive/01 - Iniciação Científica/02 - Datasets/exoplanets_confirmed/resampled_files'

import os
import pandas as pd
from datetime import datetime
import scipy.signal as ssg
import numpy as np

def get_periodograms_from_dataset(DIR_PATH):
    '''
    Calculates the periodograms for all curves on a given
    dataset, create a `pd.DataFrame` to append all the information
    and return this `pd.DataFrame`

    :param str DIR_PATH: Path to the dataset to extract
    the periodograms

    :return: `pd.DataFrame`
    '''

    DF = pd.DataFrame()

    for root_dir_path, sub_dirs, files in os.walk(DIR_PATH):
        for j in range(0, len(files)):
            if files[j].endswith('.csv'):
                path = root_dir_path + "/" + files[j]
                data = pd.read_csv(path)
                time = data.DATE
                flux = data.WHITEFLUX
                
                try: 
                    time = [datetime.strptime(i, '%Y-%m-%d %H:%M:%S.%f') for i in time]
                except:
                    time = [datetime.strptime(i, '%Y-%m-%d %H:%M:%S') for i in time]

                sample_time = pd.Series(time).diff().min()
                sample_frequency = 1 / sample_time.seconds

                # Detrend data
                detrend_flux = ssg.detrend(flux, type='linear')

                # Create the periodogram
                freq, spec = ssg.periodogram(detrend_flux, fs=sample_frequency)

                # Save the data on a pd.DataFrame
                DF = DF.append(pd.Series(spec), ignore_index=True)

        return DF


"""""""""""""""

Labeling data

"""""""""""""""

# Creating features matrix 
periodograms_confirmed_exoplanets = get_periodograms_from_dataset(RESAMPLED_CONFIRMED_EXOPLANETS_FILES_FOLDER_PATH)
periodograms_eclipsing_binaries = get_periodograms_from_dataset(RESAMPLED_ECLIPSING_BINARIES_FILES_FOLDER_PATH)

# Creating labels:
# 0: confirmed exoplanets
# 1: eclipsing binaries
#
labels_confirmed_exoplanet = np.zeros(periodograms_confirmed_exoplanets.shape[0], dtype='int')
y_confirmed_exoplanet = pd.Series(labels_confirmed_exoplanet)

labels_eclipsing_binaries = np.ones(periodograms_eclipsing_binaries.shape[0], dtype='int')
y_eclipsing_binaries = pd.Series(labels_eclipsing_binaries)


# Creating a `pd.DataFrame` with Confirmed Exoplanets Periodograms + Labels
data_labeled_confirmed_exoplanets = pd.DataFrame(periodograms_confirmed_exoplanets)
data_labeled_confirmed_exoplanets['label'] = labels_confirmed_exoplanet


# Creating a `pd.DataFrame` with Eclipsing Binaries Periodograms + Labels
data_labeled_eclipsing_binaries = pd.DataFrame(periodograms_eclipsing_binaries)
data_labeled_eclipsing_binaries['label'] = labels_eclipsing_binaries


# Creating a `pd.DataFrame` with both data labeled
data = pd.DataFrame(data_labeled_confirmed_exoplanets)
data = data.append(data_labeled_eclipsing_binaries)

print("Data shape: ", data.shape)
print(data.sample(15))

# Saving feature
where_to_save_path = 'C:/Users/guisa/Desktop/feature_periodograms.csv'
data.to_csv(where_to_save_path, index=False)

