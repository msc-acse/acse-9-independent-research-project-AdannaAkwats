import re
from datetime import date
from dateutil import rrule
from Months import Month

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


def ens_to_indx(ens_num):
    """
    Get the index related to the ensemble number : e.g 101 => 0
    :param ens_num: ensemble number, int
    :return: index, int
    """
    start = 100
    while(True):
        ind = ens_num % start
        if  ind < start:
            return ind - 1
        # Otherwise, try with bigger number of ensembles
        start *= 10

    print("Error: ens_to_index function: ensemble number cannot be converted to index")
    return None


def get_diff_start_end(start_date, end_date, monthly=False):
    """
    Returns the number of days (or months) between two dates
    :param start_date: start date ([day, month, year])
    :param end_date: end date ([day, month, year])
    :param monthly: if set, then calculate number of months, otherwise calculate number of days
    :return: the number of days (or month) between beginning of start year to start date
             the number of days (or month) between beginning of start year to end date
    """
    day_s, mon_s, yr_s = start_date[0], start_date[1], start_date[2]
    day_e, mon_e, yr_e = end_date[0], end_date[1], end_date[2]

    start, end = date(yr_s, mon_s, day_s), date(yr_e, mon_e, day_e)

    # For daily date
    if not monthly:
        # Calculate the days till the start and end
        till_start_days = (start - date(yr_s, Month.January, 1)).days
        till_end_days = (end - date(yr_s, Month.January, 1)).days
        return till_start_days, till_end_days + 1

    # For monthly data
    start, end = date(yr_s, mon_s, day_s), date(yr_e, mon_e, day_e)
    till_start_mon = len(list(rrule.rrule(rrule.MONTHLY, dtstart=date(yr_s, Month.January, 1), until=start)))
    till_end_mon = len(list(rrule.rrule(rrule.MONTHLY, dtstart=date(yr_s, Month.January, 1), until=end)))
    if mon_s == Month.January:
        till_start = 0
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