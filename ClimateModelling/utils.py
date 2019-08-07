import re
import os
from datetime import date, datetime
from dateutil import rrule
from Months import Month
import numpy as np
import directories
import ast
import sys
from collections import defaultdict


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

    # Check if files are empty
    if len(files) == 0:
        print("Error in function get_files_time_period: No NetCDF data files given within selected time period.")
        sys.exit()

    return files, min_yr, max_yr


def get_polygons(mask_file):
    """
    Load mask file which contains masks - the polygons indices (longitude, latitude) to separate the grids into regions
    :param mask_file: mask file name, string
    :return:
        nested list of  polygons indices, each list represents each mask
        level of depth

    """

    # Mask file expected to have only one list of polygons
    # Comments in file can have '#' at the start of the line

    converted, level = None, None

    # open a file using with statement
    with open(mask_file, 'r') as fh:
        for curline in fh:
            # check if the current line
            # starts with "#"
            if "#" not in curline and len(curline) > 1:
                # Get level if level key seen
                if 'level' in curline:
                    try:
                        level = curline.split(':')[1].strip()
                    except Exception:
                        print("Error in function get_polygons: Argument may be missing a semi-colon.")
                        print(
                            "Please see mask_example.out for an example of what kind of string is expected to construct "
                            "polygons.")
                        sys.exit()
                    if '-' in level:
                        try:
                            level = level.split('-')
                            level[0] = int(level[0].strip())
                            level[1] = int(level[1].strip())
                        except Exception:
                            print("Error in function get_polygons: Level number is not recognised as an integer.")
                            sys.exit()
                    else:
                        try:
                            level = int(level)
                        except Exception:
                            print("Error in function get_polygons: Level number is not recognised as an integer.")
                            sys.exit()
                # convert string to nested list of tuples
                else:
                    try:
                        converted = ast.literal_eval(curline)
                    except Exception:
                        print("Error in function get_polygons: List not constructed properly in mask file.")
                        print("Please see mask_example.out for an example of what kind of string is expected to construct "
                              "polygons.")
                        sys.exit()

    return converted, level


def make_into_file_name(str):
    """
    Remove spaces in string and replace with underscore
    :param str: string
    :return: string
    """

    # Replace all runs of whitespace with a single dash
    str = re.sub(r"\s+", '_', str)

    return str


def check_list_date(date_list):
    """
    Check that date is a list with 3 ints (day, month, yr)
    :param date: list of length of 3
    :return: boolean
    """

    return len(date_list) == 3 and all(isinstance(item, int) for item in date_list)


def get_date_from_cftime(full_date_str):
    """
    Get date without time
    :param full_date_str: string
    :return: string
    """

    match = re.search(r'(\d+-\d+-\d+)', full_date_str)
    if match:
        return str(match.group(1))


def convert_cftime_datetime(cftime_date):
    """
    Convert cftime.Datetime to datetime
    :param cftime_date: cftime.Datetime object
    :return: python datetime object
    """

    # Convert cftime date to string
    cftime_date_str = cftime_date.strftime()

    # Get only base date (with no time)
    date_str = get_date_from_cftime(cftime_date_str)

    # Convert date to datetime
    dt = datetime.strptime(date_str, "%Y-%m-%d")

    return dt


def is_nested_list(l):
    """
    Check if list is nested
    :param l: list
    :return: boolean
    """
    return any(isinstance(i, list) for i in l)


def find_middle(arr):
    """
    Gets the middle of an array
    :param arr: array
    :return: middle of array, middle index
    """
    middle = float(len(arr))/2
    if middle % 2 != 0:
        return arr[int(middle - .5)], int(middle)
    return arr[int(middle)], int(middle)


def get_shift_value(old_centre, new_centre):
    """
    Calculates how much to shift old_centre to the new_centre
    :param old_centre: float
    :param new_centre: float
    :return:
    """

    diff = old_centre - new_centre
    if old_centre < new_centre:  # <=
        if diff > 0:
            return -diff
    if new_centre > old_centre:  # =>
        if diff < 0:
            return -diff
    return diff


def shift_by_index(values, new_centre):
    """
    Get number of shifts necessary to shift values to new centre and difference between old and new centre
    :param values: values of floats
    :param new_centre: float
    :return: int
    """

    mid, mid_indx = find_middle(values)

    count = 0

    # Shift to the left
    if mid < new_centre:
        for i in range(mid_indx, len(values)):
            if values[i] <= new_centre:
                count += 1
        count = -count
    # Shift to the right
    elif mid > new_centre:
        for i in range(mid_indx + 1, -1, -1):
            if values[i] >= new_centre:
                count += 1

    return count, mid - new_centre


def get_bins_from_file(hist_file, variables):

    bin_dict = defaultdict(list)

    with open(hist_file, 'r') as hf:
        var_seen = False
        vars_in_file, bins = None, None
        for curline in hf:
            # check if the current line
            # starts with "#"
            if "#" not in curline and len(curline) > 1:
                if not var_seen:
                    # The first element should be the variable name
                    vars_in_file = curline.split(',')
                    vars_in_file = [v.strip() for v in vars_in_file]
                    # Check that variables in file match what has already been given
                    if not set(variables) == set(vars_in_file):
                        print("Error in function get_bins_from_file: Variable given in file is not recognised or "
                              "not all variables are given.")
                        print("  Please make sure that all variables given in input" + str(variables) + " are included in")
                        sys.exit()
                    var_seen = True

                else:  # Get bins
                    bins = curline.split(',')
                    if len(bins) != len(vars_in_file):
                        print("Error in function get_bins_from_file: Number of variables does not match columns of bins.")
                        print(" Number of variables: " + str(len(vars_in_file)) + ", number of bin columns: " + str(len(bins)) + ".")
                        sys.exit()
                    if len(bins) == 1:
                        try:
                            bins = [float(bins[0].strip())]
                        except Exception:
                            print("Error in function get_bins_from_file: Bin number " + bins[0] +
                                  " is not recognised as a float.")
                            sys.exit()
                    else:  # More than one bin given
                        try:
                            bins = [float(b.strip()) for b in bins]
                        except Exception:
                            print("Error in function get_bins_from_file: Bin number " + str(bins) +
                                  " is not recognised as a float.")
                            sys.exit()

                    # Add to dictionary
                    for i in range(len(bins)):
                        bin_dict[vars_in_file[i]].append(bins[i])

    return bin_dict




