import argparse
from Months import Month
from calendar import monthrange
from datetime import datetime, date
from extract import *
import directories


# Booleans set when user gives a just a year (or a year and month)
class start_bools:
    just_start_year = False
    just_start_year_month = False

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

def get_date(date_entry, start=True):
    """
    Separate date into day, month and year
    :param date_entry: string containing date e.g. 2020-04
    :param start: if set, then it is the start date of the analysis, else it is the end date
    :return: day, month, year (all integers)
    """

    date_list = []
    try:
        date = map(int, date_entry.split('-'))
        date_list = list(date)
    except ValueError:
        print("Error in function get_date(): Date written in unrecognisable format. Please try again.")
        return None

    len_d = len(date_list)
    day = 1
    month = Month.January
    year = date_list[0]
    if len_d == 1: # Only the year is given
        start_bools.just_start_year = True
        if not start:
            day = 31
            month = Month.December
    elif len_d == 2: # year and month are given
        start_bools.just_start_year_month = True
        month = date_list[1]
        if not start:
            day = monthrange(year, month)[1]
    elif len_d == 3: # day, year and month are given
        day = date_list[2]
        month = date_list[1]
    else:
        print("Error in function get_date(): too many split arguments")

    # check that these are valid dates
    try:
        datetime(year, month, day)
    except ValueError:
        print("Error in function get_date(): invalid date")
        return None

    return day, month, year


def user_entry():
    """
    Get user input
        - algae type
        - start_date
        - end_date
        - variables
        - plot option
    """
    parser = argparse.ArgumentParser(prog='CLIMATE_ANALYSIS',
                                     formatter_class=argparse.RawTextHelpFormatter,
                                     description="""The functions will give statistical analysis of the climate data presented
    FILENAMES FORMAT
    ----------------
    - The filenames should be in the format "{START OF FILENAME}_ens{NUM}_{YEAR}.nc", where {START OF FILENAME} is the prefix 
    of the file, this can be the algae type etc, {NUM} is the ensemble number and {YEAR} is the year. 
   OR if you have mutiple years stored in one file then:
   - The filenames should be in the format "{START OF FILENAME}_ens{NUM}_{YEAR 1}_{YEAR 2}.nc", where {START OF FILENAME} 
   is the prefix of the file, this can be the algae type etc, {NUM} is the ensemble number and {YEAR 1} and {YEAR 2} are 
   the start and end year of the data in the file. 
   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    - You have a mixture of files with only one year or with multiple years (i.e. the data folder can have files with the 
    different formats at the same time).
    - We assume files do not have overlapped data.
    - Some example files are in the data folder.
    - We assume daily increments of data, except if the monthly tag is set in the arguments.
    """)
    parser._optionals.title = "other arguments"
    parser.add_argument('algae_type', help="Type of algae e.g. dic_deltap, fndet_100, jprod_ndi_100. This can also just be the prefix of the file.")
    parser.add_argument('start_date', help="""Start date of analysis 
    Can be in the following formats:
    ----------------------------------
    YYYY-MM-DD : e.g. 2020-04-12
    YYYY-MM    : e.g. 2020-04
    YYYY       : e.g. 2020 
    - If day is not given, the 1st of the given month will be used i.e 2020-04 => 2020-04-01
    - If day and month is not given, 1st Jan will be used as the start date i.e 2020 => 2020-01-01""")
    parser.add_argument('end_date', nargs='?', help=""" <Not required> End date of analysis - format is the same as start_date
    -----------------------------------end_date not given-------------------------------------
    - If only start year is given, the end_date is automatically set to the 31 Dec of start year
    - If start year and month is given, then end_date is set to the end of the start month
       -----------------------------------end_date given-------------------------------------
    - If day is not given, the end of the given month will be used i.e 2020-04 => 2020-04-30
    - If day and month is not given, 31 Dec will be used as the end date i.e 2020 => 2020-12-31""")
    parser.add_argument('-v', '--vars', nargs='+', metavar=("variables"), help="<Required> Variables of data to analyse", required=True)
    parser.add_argument('-p', '--plot', action="store_true", help="Save plots of analysis in " + directories.ANALYSIS + " as a .png file.")
    parser.add_argument('-m', '--monthly', action="store_true", help="Data in file is stored in monthly increments.")
    group = parser.add_mutually_exclusive_group()
    group.add_argument('-g', '--grid', nargs=2, type=float, metavar=("lat","lon"), help="Uses grid point that latitude and longitude lies in.")
    group.add_argument('-s', '--sample', nargs=2, type=float, metavar=("lat","lon"), help="Uses sample point given by latitude and longitude using interpolation.")
    parser.add_argument('-mk', '--mask', nargs=1, metavar=("filename"), help="Uses masking grid given as a file (contains boolean array to be imposed on the global grid).")
    parser.add_argument('-o', '--output', action="store_true", help="Save data output of histogram and timeseries analysis in " + directories.ANALYSIS + " as a .dat file.")
    parser.add_argument('-cv', '--covary', action="store_true", help="Analysis on how the variables given in -v vary with each other.")
    parser.add_argument('-e', '--ens', nargs='?', type=int, metavar="number_of_ensembles", help="The number of ensemebles of the data. If not set, the default value = 1", const=1, default=1)

    # Arguments
    args = parser.parse_args()
    type = args.algae_type
    vars = args.vars
    start = args.start_date
    end = args.end_date

    # Get split start date
    day_s, mon_s, yr_s = get_date(start)
    if not end: # If end date not given, use the end of start year
        if start_bools.just_start_year:
            end = str(yr_s)
        elif start_bools.just_start_year_month:
            end = str(yr_s) + "-" + str(mon_s)

    # Get split end date
    day_e, mon_e, yr_e = get_date(end, start=False)

    # Print user input
    print("Arguments:")
    print("- algae type: ", type)
    print("- variables: ", vars)
    print("- start date: " + str(yr_s) + "-" + str(mon_s) + "-" + str(day_s))
    print("- end date: " + str(yr_e) + "-" + str(mon_e) + "-" + str(day_e))

    # Check that dates are in valid order
    is_valid = check_valid_order([day_s, mon_s, yr_s], [day_e, mon_e, yr_e])
    if not is_valid:
        print("Error: Invalid start and end date")
        print("  - The end date is earlier than the start date")
        sys.exit()
    print("Number of ensembles:", args.ens)
    if args.plot:
        print("Plotting option selected.")
    if args.monthly:
        print("Monthly date expected.")

    lat, lon = None, None
    if args.grid:
        lat, lon = args.grid[0], args.grid[1]
        print("Grid point option selected.")
    if args.sample:
        lat, lon = args.sample[0], args.sample[1]
        print("Sample point option selected.")

    mask = None
    if args.mask:
        mask = args.mask[0]
        print("Masking grid option selected.")
    if args.output:
        print("Save analysis data output selected.")
    if args.covary:
        print("Co-varying option selected.")

    # Call functions to perform analysis
    from utils import overlaps

    saved, units = extract_data(type, vars, [day_s, mon_s, yr_s], [day_e, mon_e, yr_e], num_ens=args.ens, lat=lat, lon=lon, mask=mask)
    # analysis(saved, units, monthly=args.monthly, save_out=args.output, cov=args.covary)

    # if args.plot:
    #     plot()



    # vars = ['air_temperature']
    # saved, units = extract_data(vars, None)
    # analysis(saved, units)
    # plot("results/means.nc") # This should be automatically off - but put as an option

