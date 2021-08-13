"""
Started in: August 04, 2021
Finished in: 

@author: Guilherme Samuel
"""
from os import system
system("cls")

ECLIPSING_BINARIES_FILES_FOLDER_PATH = 'C:/Users/guisa/Desktop/eclipsing_binaries_csv'
CONFIRMED_EXOPLANETS_FILES_FOLDER_PATH = 'C:/Users/guisa/Google Drive/01 - Iniciação Científica/02 - Datasets/exoplanets_confirmed/csv_files'

import os
import pandas as pd
from statistics import median
from datetime import datetime
import scipy.signal as ssg
import shutil

def get_median_sample_size(DIR_PATH):
    '''
    Computes the sample size for each curve and 
    returns the median one

    :param str DIR_PATH: Path to the dataset to extract
    the median sample size

    :return: Median sample size
    '''
    sample_size = []
    for root_dir_path, sub_dirs, files in os.walk(DIR_PATH):
        for j in range(0, len(files)):
            if files[j].endswith('.csv'): 
                path = root_dir_path + "/" + files[j]
                data = pd.read_csv(path)
                
                sample_size.append(len(data.WHITEFLUX))

    # Uncomment to see details
    # print("Minimum size:", min(sample_size))
    # print("Median size:", median(sample_size))
    # print("Maximum size:", max(sample_size))

    return int(median(sample_size))

def resampling_dataset(DIR_PATH, RESAMPLE_DIR_PATH, sample_size):
    '''
    Resamples an entire given dataset and save the results
    on a given directory

    :param str DIR_PATH: Path to the dataset to resample
    :param str RESAMPLE_DIR_PATH: Path to the new resampled dataset
    :param int sample_size: Sample size to resampled all dataset
    '''
    count = 0
    for root_dir_path, sub_dirs, files in os.walk(DIR_PATH):
        for j in range(0, len(files)):
            if files[j].endswith('.csv'):
                path = root_dir_path + "/" + files[j]
                data = pd.read_csv(path)
                                
                # If sample size is less than 584, the lightcurve will be discarded
                if (data.WHITEFLUX.size > 584):
                    flux = data.WHITEFLUX
                    time = data.DATE

                    try: 
                        time = [datetime.strptime(i, '%Y-%m-%d %H:%M:%S.%f') for i in time]
                    except:
                        time = [datetime.strptime(i, '%Y-%m-%d %H:%M:%S') for i in time]

                    flux_resampled, time_resampled = ssg.resample(flux, sample_size, time)

                    # Creating a new pd.DataFrame
                    concat_dict = {
                    "DATE": pd.Series(time_resampled), 
                    "WHITEFLUX": pd.Series(flux_resampled)
                    }
                    data_resampled = pd.concat(concat_dict, axis=1)

                    # Creating folder with lightcurves resampled
                    if not os.path.isdir(RESAMPLE_DIR_PATH):
                        os.mkdir(RESAMPLE_DIR_PATH)

                    # Renaming lightcurve
                    file_name = 'RESAMPLED_' + files[j]

                    # Saving lightcurves resampled 
                    FILE_DIR = file_name
                    data_resampled.to_csv(file_name, index=False)

                    shutil.move(FILE_DIR, RESAMPLE_DIR_PATH)
                    
                    count += 1
                    # Comment to hide details
                    print('Resampled and saved: ' + files[j])
    print("\nTotal of files resampled:", count)

def main():
    """
    Function that apply the `resampling_dataset`
    to an decided dataset
    """
    # Sample size default will be the median sample size for the confirmed exoplanets dataset
    sample_size = get_median_sample_size(CONFIRMED_EXOPLANETS_FILES_FOLDER_PATH)

    # Calling the function to resample the dataset
    resampling_dataset(ECLIPSING_BINARIES_FILES_FOLDER_PATH, 'C:/Users/guisa/Desktop/resampled_eclipsing_binaries', sample_size)

if __name__ == '__main__':
    main()