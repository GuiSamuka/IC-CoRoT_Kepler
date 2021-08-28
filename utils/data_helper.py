"""

    This module implements some algorithms that help us
    read and manipulate `.fits` files and `.csv` datasets 

"""

from os import system
system("cls")

import os
from astropy.io import fits
import numpy as np
import pandas as pd
import shutil
from statistics import median
import scipy.signal as ssg
# from lightcurve import *

def fits_to_csv(FITS_PATH, CSV_PATH):
    """
    Given a dataset, formed by `.fits` files,
    convert it all to a new dataset, but with 
    equivalent `.csv` files

    Parameters
    ----------
    FITS_PATH : String
        Path to .fits data-set folder

    CSV_PATH : String
        Path to the new folder, with the csv files
    """

    # if csv folder do not exist, create it
    if not os.path.isdir(CSV_PATH):
        os.mkdir(CSV_PATH)

    for root_dir_path, sub_dirs, files in os.walk(FITS_PATH):
        for i in range(0, len(files)):
            path = FITS_PATH + '/' + files[i]
            
            if path.endswith('.fits'):
                # ``image_file`` contains all the information storage into .fits file
                image_file = fits.open(path)
                
                # ``scidata`` contains the information of the third table of the .fits file
                scidata = image_file[3].data

                # ``data_aux`` contains the same info as ``scidata``, but in a np.array format
                data_aux = np.array(scidata).byteswap().newbyteorder()

                # ``data`` contains the same info as ``data_aux``, but in a pd.DataFrame format
                data = pd.DataFrame(data_aux)

                # We don't need STATUSSYS column
                data.drop('STATUSSYS', axis=1, inplace=True)

                # Renaming columns
                data.rename(columns={'DATEBARTT': 'DATE'}, inplace=True)
                data.rename(columns={'WHITEFLUXSYS': 'WHITEFLUX'}, inplace=True)

                # Verify if there's any Not a Number (NaN) value
                if (data.isnull().values.any()):
                    print("There's some NaN value on data. Check the data!")
                    print("Name file:", files[i])
                    break

                # Renaming .csv files
                name = path[path.rfind('/')+1:path.rfind('.')] + '.csv'
                name = name.split('_')[3] + "_" + name.split('_')[4] + ".csv"

                # Saving data
                data.to_csv(name, index=False)
                print("Data saved:", name)

                # Move to .csv folder
                shutil.move(name, CSV_PATH)

    # Verify if all files have been converted                
    path, dirs, files = next(os.walk(FITS_PATH))
    num_fits = 0
    for i in files:
        if i.endswith('.fits'):
            num_fits += 1

    path, dirs, files = next(os.walk(CSV_PATH))
    num_csv = 0
    for i in files:
        if i.endswith('.csv'):
            num_csv += 1

    if num_fits == num_csv:
        print("All files have been converted successfully!")
    else:
        print("Not all files have been converted! Uncomment the `print` statement on line 64 to see details")
    

def get_median_sample_size(CSV_PATH):
    """
    Computes the sample size for each curve and 
    returns the median one

    Parameters
    ----------
    CSV_PATH : String
        Path to .csv data-set folder

    Return
    ------
    Median sample size
    """
    sample_size = []
    for root_dir_path, sub_dirs, files in os.walk(CSV_PATH):
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


def resampling_dataset(CSV_PATH, RESAMPLE_PATH, sample_size):
    """
    Resamples an entire given dataset and save the results
    on a given directory

    Parameters
    ----------
    CSV_PATH : String
        Path to .csv data-set folder
    """

    count = 0
    for root_dir_path, sub_dirs, files in os.walk(CSV_PATH):
        for j in range(0, len(files)):
            if files[j].endswith('.csv'):
                path = root_dir_path + "/" + files[j]
                data = pd.read_csv(path)

                # If sample size is less than 584, the lightcurve will be discarded
                if (data.WHITEFLUX.size > 584):
                    flux = data.WHITEFLUX
                    time = data.DATE

                    flux_resampled, time_resampled = ssg.resample(flux, sample_size, time)

                    # Creating a new pd.DataFrame
                    concat_dict = {
                    "DATE": pd.Series(time_resampled), 
                    "WHITEFLUX": pd.Series(flux_resampled)
                    }
                    data_resampled = pd.concat(concat_dict, axis=1)

                    # Creating folder with lightcurves resampled
                    if not os.path.isdir(RESAMPLE_PATH):
                        os.mkdir(RESAMPLE_PATH)

                    # Renaming lightcurve
                    file_name = 'RESAMPLED_' + files[j]

                    # Saving lightcurves resampled 
                    FILE_DIR = file_name
                    data_resampled.to_csv(file_name, index=False)

                    shutil.move(FILE_DIR, RESAMPLE_PATH)
                    count += 1

                    # Comment to hide details
                    # print('Resampled and saved: ' + files[j])

    print("\nTotal of files resampled:", count)


def main():
    FITS_FILES_PATH = r'C:\Users\guisa\Google Drive\01 - Iniciação Científica\02 - Datasets\exoplanets_confirmed\dataset_exoplanets_confirmed'

    CSV_FILES_PATH = r'C:\Users\guisa\Google Drive\01 - Iniciação Científica\02 - Datasets\exoplanets_confirmed\csv_files'

    RESAMPLED_FILES_PATH = r'resampled_files'

    "First step: convert all files to a csv format"
    # fits_to_csv(FITS_FILES_PATH, 'C:/Users/guisa/Desktop/csv_files')

    "Second step: Resampling data"
    # sample_size = get_median_sample_size(CSV_FILES_PATH)
    # resampling_dataset(CSV_FILES_PATH, RESAMPLED_FILES_PATH, sample_size)

    

if __name__ == '__main__':
    main()
