"""
Created by Adanna Akwataghibe (Github: AdannaAkwats)
"""
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
from calendar import isleap
"""
Script that contains useful functions
"""


def check_valid_order(start_date, end_date):
    """
    Checks that end date is after start date
    :param start_date: [day, month, year]
    :param end_date: [day, month, year]
    :return: true if end_date is after start date
    """

    day_s, mon_s, yr_s = start_date[0], start_date[1], start_date[2]
    day_e, mon_e, yr_e = end_date[0], end_date[1], end_date[2]

    start = date(yr_s, mon_s, day_s)
    end = date(yr_e, mon_e, day_e)

    return (end - start).days > 0


def check_analysis(a):
    """
    Check analysis a is valid
    :param a: list of strings
    :return: boolesn or throws error
    """

    assert(isinstance(a, list))

    # List of analysis that can be selected
    ans = ['mean', 'rms', 'std', 'median']

    # Check that analysis a is in ans
    all_elements_contained = all(elem.lower() in ans for elem in a)

    if not all_elements_contained:
        print("ERROR in function check_analysis: Analysis given cannot be selected. "
              "Ensure that analysis is one or more of [mean, std, rms, median].")
        sys.exit()
    return True


def check_variables_covary(varbs):
    """
    Check only 2 variables given - this function is called if covariance is set
    :param varbs: list of variables (string)
    :return: boolean or throws error
    """
    assert (isinstance(varbs, list))
    if len(varbs) != 2:  # If not 2 variables, then cannot perform covariance
        print("ERROR in function check_variables_covary: Selecting covariance required two variables.")
        sys.exit()
    return True


def check_valid_indices(index):
    """
    Check that index given is valid
    :param index: string
    :return: boolean or throws error
    """

    # List of indices that can be selected
    indices = ['enso', 'nino12', 'nino4', 'tni', 'iod', 'amo', 'pdo', 'ao', 'aao', 'nao']

    if index not in indices:
        print("ERROR in function check_valid_indices: Index " + str(index) + " given cannot be selected. "
              "Ensure that index is one of " + str(indices) + ".")
        sys.exit()
    return True



def get_ens_num(file):
    """
    Return the ensemble number in file name
    :param file: file name, string
    :return: ensemble number, int
    """
    f = r'ens(\d+)'
    match = re.search(f, file)
    if match:
        return int(match.group(1))
    # If no 'ens' seen, assume only one ensemble given
    return 101


def get_file_two_years(file):
    """
    Returns the years that data is within in file
    :param file: file name, string
    :return: years, ints
    """
    f = r'_(\d+)_(\d+)'

    match = re.search(f, file)
    if match:
        # Check strings are length 4 - years
        if len(match.group(1)) >= 4 and len(match.group(2)) >= 4:
            return int(match.group(1)[:4]), int(match.group(2)[:4])

    f = r'_(\d+)-(\d+)'
    match = re.search(f, file)
    if match:
        # Check strings are length 4 - years
        if len(match.group(1)) >= 4 and len(match.group(2)) >= 4:
            return int(match.group(1)[:4]), int(match.group(2)[:4])
    return False


def ens_to_indx(ens_num, number_of_ensembles, max_start=1000000):
    """
    Get the index related to the ensemble number : e.g 101 => 0
    :param ens_num: ensemble number, int
    :param number_of_ensembles, number of total ensembles used, int
    :param max_start: max number of ensembles, int
    :return: index, int
    """
    start = 100
    res = -1
    while start < max_start:
        ind = ens_num % start
        if ind < start:
            res = ind - 1
            break
        # Otherwise, try with bigger number of ensembles
        start *= 10

    # Check if its within the number of ensembles given
    if res == -1:
        print("ERROR: ens_to_index function: ensemble number cannot be converted to index")
        sys.exit()
    elif res + 1 > number_of_ensembles:
        res = -1
    return res


def get_diff_start_end(start_date, end_date, min_yr=None, monthly=False, num_leap_year_input=None):
    """
    Returns the number of days (or months) between two dates
    :param start_date: start date ([day, month, year])
    :param end_date: end date ([day, month, year])
    :param min_yr: minimum year, int, default = None
    :param monthly: if set, then calculate number of months, otherwise calculate number of days
    :param num_leap_year_input: number of leap years in the input dataset (in extract)
    :return: the number of days (or month) between beginning of start year to start date
             the number of days (or month) between beginning of start year to end date
    """
    day_s, mon_s, yr_s = start_date[0], start_date[1], start_date[2]
    day_e, mon_e, yr_e = end_date[0], end_date[1], end_date[2]

    if not min_yr:
        min_yr = yr_s

    start, end = date(yr_s, mon_s, day_s), date(yr_e, mon_e, day_e)

    # Calculate the number of leap years between min date and start date
    start_num_leap_date, end_num_leap_date = 0, 0
    if not monthly:
        for i in range(min_yr, yr_s):
            if isleap(i):
                start_num_leap_date += 1
        #  Calculate the number of leap years between start date and end date
        for i in range(yr_s, yr_e+ 1):
            if isleap(i):
                end_num_leap_date += 1

    # For daily date
    if not monthly:
        # Calculate the days till the start and end
        till_start_days = (start - date(min_yr, Month.January, 1)).days
        till_end_days = (end - date(min_yr, Month.January, 1)).days
        if num_leap_year_input == 0:  # If calendar is NOLEAP (365-day)
            # remove leap year day from days
            till_start_days -= start_num_leap_date
            till_end_days -= end_num_leap_date
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
    path = directories.DATA + '/'
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
        print("ERROR in function get_files_time_period: No NetCDF data files given within selected time period.")
        print("  - Please ensure that the start and end years given are the same as in the file name.")
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
                        print("ERROR in function get_polygons: Argument may be missing a semi-colon.")
                        print(
                            "Please see mask_example.txt for an example of what kind of string is expected to construct "
                            "polygons.")
                        sys.exit()
                    if '-' in level:
                        try:
                            level = level.split('-')
                            level[0] = int(level[0].strip())
                            level[1] = int(level[1].strip())
                        except Exception:
                            print("ERROR in function get_polygons: Level number is not recognised as an integer.")
                            sys.exit()
                    else:
                        try:
                            level = [int(level)]
                        except Exception:
                            print("ERROR in function get_polygons: Level number is not recognised as an integer.")
                            sys.exit()
                # convert string to nested list of tuples
                else:
                    try:
                        converted = ast.literal_eval(curline)
                    except Exception:
                        print("ERROR in function get_polygons: List not constructed properly in mask file.")
                        print("Please see mask_example.txt for an example of what kind of string is expected to construct "
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
    cftime_date_str = cftime_date.strftime("%Y-%m-%d")

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
        return arr[int(middle - .5)], int(middle - .5)
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
    """
    Get bin edges from file
    :param hist_file: txt file
    :param variables: list of variables
    :return: dictionary key (variable) : value (bin edges)
    """

    bin_dict = defaultdict(list)
    x, y = None, None

    with open(hist_file, 'r') as hf:
        var_seen, cur_var_name = False, None
        vars_seen_in_file = []
        for curline in hf:
            # check if the current line
            # starts with "#"
            if "#" not in curline and len(curline.strip()) > 1:
                key_var = curline.split(':')
                if key_var[0].strip() == 'var':
                    # Get variable name
                    cur_var_name = key_var[1].strip()
                    if len(cur_var_name) < 1:
                        print("ERROR in function get_bins_from_file: No variable name given.")
                        sys.exit()
                    # check variable is in list of variables
                    if cur_var_name not in variables:
                        print("ERROR in function get_bins_from_file: Variable '" + cur_var_name+ "' is not recognised or " +
                              "not all variables are given.")
                        print("  Please make sure that all variables given in input " + str(
                            variables) + " are included in the histogram bins file.")
                        sys.exit()
                    # Check if variable defined more than once in file
                    if cur_var_name in vars_seen_in_file:
                        print("ERROR in function get_bins_from_file: Multiple definitions of " +
                              str(cur_var_name) + " in the histogram bins file.")
                        sys.exit()
                    vars_seen_in_file.append(cur_var_name)
                    continue
                elif key_var[0].strip() == 'x':
                    x = key_var[1].strip()
                    if x not in variables:
                        print("ERROR in function get_bins_from_file: Variable '" + x + "' is not recognised.")
                        print("  Please make sure that x variable is one of " + str(variables))
                        sys.exit()
                elif key_var[0].strip() == 'y':
                    y = key_var[1].strip()
                    if y not in variables:
                        print("ERROR in function get_bins_from_file: Variable '" + y + "' is not recognised.")
                        print("  Please make sure that x variable is one of " + str(variables))
                        sys.exit()
                else:
                    try:
                        bin_ = float(curline.strip())
                        bin_dict[cur_var_name].append(bin_)
                    except ValueError:
                        print("ERROR in function get_bins_from_file: bin number " + str(curline) + " is not recognised "
                                                                                                   "as a number.")
                        sys.exit()


    if x is None and y is None:
        print("ERROR in function get_bins_from_file: x-axis and y-axis is not defined.")
    return bin_dict, x, y


def print_end_statement():
    print(">> PROGRAM FINSISHED.")
    print(">> Progress and potential errors are logged in output.log file.")
    print(">> To open output.log, type in cmd line: less output.log")


def open_txt_points(filename):
    """
    Gets latitudes and longitudes from sample/grid point file
    :param filename: points txt file to open
    :return: lons, lats
    """
    converted = None
    with open(filename, 'r') as f:
        for curline in f:
            # check if the current line
            # starts with "#"
            if "#" not in curline and len(curline) > 1:
                try:
                    converted = ast.literal_eval(curline)
                except Exception:
                    print("ERROR in function open_txt_points: List not constructed properly in sample/grid points file.")
                    print("Please see sample_points.txt for an example of what kind of string is expected to construct "
                          "sample points.")
                    sys.exit()


    # Get list of latitude and longitude
    lons, lats = zip(*converted)
    return lons, lats


def get_sample_grid_points(filename):
    """
    Get sample/grid points from filename
    :param filename: file to open
    :return: longitudes and latitudes and setting if nc file given
    """
    lons, lats, nc_true = None, None, False
    # Check if txt file
    if filename.endswith('.txt'):
        filename = os.path.join(directories.INPUT, filename)
        lons, lats = open_txt_points(filename)
    elif filename.endswith('.nc'):
        nc_true = True
    return lons, lats, nc_true


def check_sample_grid_one_arg(arg, func):
    """
    Check one argument of sample / grid
    :param arg: argument
    :param func: function calling this
    :return: error or boolean
    """
    sg = arg[0]
    if not (sg.endswith(".nc") or sg.endswith(".txt")):
        print("ERROR in function " + str(func) + ": Sample or grid argument invalid.")
        print("  - Argument must be (lat, lon) or .txt file or .nc file")
        sys.exit()
    return True