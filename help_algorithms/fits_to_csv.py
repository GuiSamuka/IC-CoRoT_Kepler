"""
Started in: August 03, 2021
Finished in: 

@author: Guilherme Samuel
"""
import shutil
import os
import time
from datetime import datetime
import julian
import pandas as pd
import numpy as np
from astropy.io import fits
from os import system
system("cls")


FITS_FILES_FOLDER_PATH = 'C:/Users/guisa/Desktop/Eclipsing_Binaries'

FILE_NAME = FITS_FILES_FOLDER_PATH + \
    '/EN2_STAR_CHR_0102582649_20071023T223035_20080303T093534.fits'


def julian_to_stdtime(old_date):
    """
    This is the function used to convert Julian
    date to Gregorian date

    :param numpy.float64 old_date: Represents Julian date of float type
    """
    aux_1 = julian.from_jd(old_date, fmt='mjd')
    try:
        aux_2 = datetime.strptime(str(aux_1), '%Y-%m-%d %H:%M:%S.%f')
    except:
        aux_2 = datetime.strptime(str(aux_1), '%Y-%m-%d %H:%M:%S')

    new_date = str(aux_2)

    return new_date


def fits_to_csv(path, CSV_DIR):
    '''
    Normalize one .fits files and convert it
    into a .csv file

    :param str path: Path to .fits data-set folder

    :param str CSV_DIR: Path to the new folder, with the csv files
    '''

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
        print("Name file:", FILE_NAME)

    data.DATE = data.DATE.apply(lambda x: julian_to_stdtime(x))

    # Creating folder with .csv files
    if not os.path.isdir(CSV_DIR):
        os.mkdir(CSV_DIR)

    # Renaming .csv files
    name = path[path.rfind('/')+1:path.rfind('.')] + '.csv'

    # Saving data
    data.to_csv(name, index=False)
    print("Data saved:", name)

    # Move to .csv folder
    shutil.move(name, CSV_DIR)


def main():
    """
    Function that apply the `fits_to_csv`
    to all .fits dataset (FITS_FILES_FOLDER_PATH)
    """
    i = 0
    my_dir = FITS_FILES_FOLDER_PATH
    t0 = time.time()

    for root_dir_path, sub_dirs, files in os.walk(my_dir):
        for i in range(0, len(files)):
            path = my_dir + '/' + files[i]

            if path.endswith('.fits'):
                fits_to_csv(
                    path, "C:/Users/guisa/Desktop/eclipsing_binaries_csv")
                # break

    tf = time.time()
    print("\nIt takes:", round(tf-t0, 2), "seconds to apply")


if __name__ == '__main__':
    main()
