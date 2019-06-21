import re
import os
from datetime import date
from dateutil import rrule
from Months import Month
import numpy as np
import directories


"""
Script that contains useful functions
"""


def get_ens_num(file):
    """
    Return the ensemble number in file name
    :param file: file name, string
    :return: ensemble number, int
    """
    f = 'ens' + '(\d+)'
    match = re.search(f, file)
    if match:
        return int(match.group(1))


def get_file_two_years(file):
    """
    Returns the years that data is within in file
    :param file: file name, string
    :return: years, ints
    """
    f = '_' + '(\d+)' + '_' + '(\d+)'
    match = re.search(f, file)
    if match:
        return int(match.group(1)), int(match.group(2))


def ens_to_indx(ens_num, max_start=1000000):
    """
    Get the index related to the ensemble number : e.g 101 => 0
    :param ens_num: ensemble number, int
    :param max_start: max number of ensembles, int
    :return: index, int
    """
    start = 100
    while start < max_start:
        ind = ens_num % start
        if ind < start:
            return ind - 1
        # Otherwise, try with bigger number of ensembles
        start *= 10

    print("Error: ens_to_index function: ensemble number cannot be converted to index")


def get_diff_start_end(start_date, end_date, min_yr=None, monthly=False):
    """
    Returns the number of days (or months) between two dates
    :param start_date: start date ([day, month, year])
    :param end_date: end date ([day, month, year])
    :param min_yr: minimum year, int, default = None
    :param monthly: if set, then calculate number of months, otherwise calculate number of days
    :return: the number of days (or month) between beginning of start year to start date
             the number of days (or month) between beginning of start year to end date
    """
    day_s, mon_s, yr_s = start_date[0], start_date[1], start_date[2]
    day_e, mon_e, yr_e = end_date[0], end_date[1], end_date[2]

    if not min_yr:
        min_yr = yr_s

    start, end = date(yr_s, mon_s, day_s), date(yr_e, mon_e, day_e)

    # For daily date
    if not monthly:
        # Calculate the days till the start and end
        till_start_days = (start - date(min_yr, Month.January, 1)).days
        till_end_days = (end - date(min_yr, Month.January, 1)).days
        return till_start_days, till_end_days + 1

    # For monthly data
    start, end = date(yr_s, mon_s, day_s), date(yr_e, mon_e, day_e)
    till_start_mon = len(list(rrule.rrule(rrule.MONTHLY, dtstart=date(min_yr, Month.January, 1), until=start)))
    till_end_mon = len(list(rrule.rrule(rrule.MONTHLY, dtstart=date(min_yr, Month.January, 1), until=end)))
    if mon_s == Month.January and yr_s == min_yr:
        till_start_mon = 0
    return till_start_mon, till_end_mon


def overlaps(x1, x2, y1, y2):
    """
    Returns true if array [x1, x2] overlaps with [y1, y2]
    :param x1: int
    :param x2: int, assume x1 <= x2
    :param y1: int
    :param y2: int, assume y1 <= y2
    :return: boolean
    """

    return x1 <= y2 and y1 <= x2


def find_nearest(array, value):
    """
    Returns the index of the element closest to value in array
    :param array: numpy array
    :param value: float
    :return: index, int
    """
    array = np.asarray(array)
    idx = (np.abs(array - value)).argmin()
    return idx


def generate_example_mask(lat_size, lon_size):
    """
    Generates random mask and saves to a mask_example.txt
    :param lat_size: latitude
    :param lon_size: longitude
    :return: random boolean 2d array, saves output to file
    """
    mask_rand = np.random.randint(2, size=(lat_size, lon_size))
    save_path = directories.ANALYSIS + '/'
    np.savetxt(os.path.join(save_path, 'mask_example.out'), mask_rand)


def get_files_time_period(prefix, yr_s, yr_e):
    """
    Get nc files within the time period, with the prefix
    :param prefix: prefix of file names
    :param yr_s: start year
    :param yr_e: end year
    :return: files in folder, max and min year
    """

    # Get path and folder
    path = directories.CLIMATE_DATA + '/'
    folder = os.listdir(path)

    # Files should be automatically ordered by year assuming that the format of files is what we expect
    files = []

    # List of years to extract
    years = list(range(yr_s, yr_e + 1))

    # Save lowest and highest year in data for later - only used if multiple years are in the same file
    min_yr = yr_s
    max_yr = yr_e

    # Go through the files in the folder and get the relevant files within the time frame
    for file in folder:
        if os.path.isfile(os.path.join(path, file)) and file.startswith(prefix):
            # If file with just one year in it
            if not get_file_two_years(file):
                for year in years:
                    if str(year) in file:
                        files.append(file)
            else:  # file has multiple years in it
                fst_yr, snd_yr = get_file_two_years(file)
                # Get files that have data within the years
                if overlaps(fst_yr, snd_yr, yr_s, yr_e):
                    files.append(file)
                    if fst_yr < min_yr:
                        min_yr = fst_yr
                    if snd_yr > max_yr:
                        max_yr = snd_yr

    return files, min_yr, max_yr
