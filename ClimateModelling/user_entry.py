import sys
import argparse
from Months import Month
from calendar import monthrange
from datetime import datetime, date
from extract import *
import directories


# Booleans set when user gives a just a year (or a year and month)
class StartBools:
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

    try:
        date_ = map(int, date_entry.split('-'))
        date_list = list(date_)
    except ValueError:
        print("Error in function get_date(): Date written in unrecognisable format. Please try again.")
        return None

    len_d = len(date_list)
    day = 1
    month = Month.January
    year = date_list[0]
    if len_d == 1:  # Only the year is given
        StartBools.just_start_year = True
        if not start:
            day = 31
            month = Month.December
    elif len_d == 2:  # year and month are given
        StartBools.just_start_year_month = True
        month = date_list[1]
        if not start:
            day = monthrange(year, month)[1]
    elif len_d == 3:  # day, year and month are given
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


def file_entry(example=False):
    filename = directories.INPUT_FILE
    if example:
        filename = directories.INPUT_EXAMPLE_FILE

    # Save arguments
    args = []
    # open a file using with statement
    with open(filename, 'r') as fh:
        for curline in fh:
            # check if the current line
            # starts with "#"
            if "#" not in curline:
                arg = curline.split(':')[1].strip()
                args.append(arg)

    # Check if any required arguments are not filled in
    for i in range(4):
        if len(args[i]) == 0:
            print("Error: the input file has missing required arguments.")
            sys.exit()

    # Check if any optional arguments are not filled in
    for i in range(4, 12):
        if args[i].lower() == 'false' or args[i].lower() == 'f':
            args[i] = False
        elif args[i].lower() == 'true' or args[i].lower() == 't':
            args[i] = True

    # histogram option
    if len(args[12]) == 0:
        hist = 'fd'

    algae_type, start, ens, end = args[0], args[1], int(args[3]), args[4]
    plot, monthly, grid, sample, mask, output, covary, hist = args[5], args[6], args[7], args[8], args[9], args[10], args[11], args[12]

    # Check for grid and sample
    if args[7]:
        grid = args[7].split(',')
        grid[0] = float(grid[0].strip())
        grid[1] = float(grid[1].strip())
    elif args[8]:
        sample = args[8].split(',')
        sample[0] = float(sample[0].strip())
        sample[1] = float(sample[1].strip())

    # Get variables and put in list
    varbs = args[2].split(',')
    for i in range(len(varbs)):
        varbs[i] = varbs[i].strip()

    return algae_type, start, varbs, ens, end, plot, monthly, grid, sample, mask, output, covary, hist


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
                                     description="""The functions will give statistical analysis of the climate data 
                                     presented
    FILENAMES FORMAT
    ----------------
    - The filenames should be in the format "{START OF FILENAME}_ens{NUM}_{YEAR}.nc", where {START OF FILENAME} is 
    the prefix of the file, this can be the algae type etc, {NUM} is the ensemble number and {YEAR} is the year. 
   OR if you have multiple years stored in one file then:
   - The filenames should be in the format "{START OF FILENAME}_ens{NUM}_{YEAR 1}_{YEAR 2}.nc", where 
   {START OF FILENAME} is the prefix of the file, this can be the algae type etc, {NUM} is the ensemble number and 
   {YEAR 1} and {YEAR 2} are the start and end year of the data in the file. 
   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    ASSUMPTIONS
    ------------
    - Files do not have overlapped data.
    - Daily increments of data, except if the monthly tag is set in the arguments.
    - Grids have constant latitude and longitude.
    ------------
    - Some example files are in the data folder.
    """)
    parser._optionals.title = "other arguments"
    parser.add_argument('prefix', help="This is the prefix of the file. This can be the type of algae "
                                       "e.g. dic_deltap, fndet_100, jprod_ndi_100. ")
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
    parser.add_argument('-v', '--vars', nargs='+', metavar="variables", help="<Required> Variables of data to analyse",
                        required=True)
    parser.add_argument('-p', '--plot', action="store_true", help="Save plots of analysis in " + directories.ANALYSIS +
                                                                  " as a .png file.")
    parser.add_argument('-m', '--monthly', action="store_true", help="Data in file is stored in monthly increments.")
    group = parser.add_mutually_exclusive_group()
    group.add_argument('-g', '--grid', nargs=2, type=float, metavar=("lat", "lon"), help="Uses grid point that "
                                                                                         "latitude and longitude lies "
                                                                                         "in.")
    group.add_argument('-s', '--sample', nargs=2, type=float, metavar=("lat", "lon"), help="Uses sample point given by"
                                                                                           " latitude and longitude "
                                                                                           "using interpolation.")
    parser.add_argument('-mk', '--mask', nargs=1, metavar="filename", help="Uses masking grid given as a file "
                                                                           "(contains boolean array to be imposed on "
                                                                           "the global grid).")
    parser.add_argument('-o', '--output', action="store_true", help="Save data output of histogram and timeseries "
                                                                    "analysis in "
                                                                    + directories.ANALYSIS + " as a .dat file.")
    parser.add_argument('-cv', '--covary', action="store_true", help="Analysis on how the variables given in -v "
                                                                     "vary with each other.")
    parser.add_argument('-e', '--ens', nargs=1, type=int,
                        metavar="number_of_ensembles", help="<Required> The number of ensembles of the data. "
                                                            "If not set, the default value = 1", required=True)
    parser.add_argument('-ht', '--hist', nargs='?', const='fd', default='fd',
                        metavar="number_of_bins_in_histogram", help=" Options for bin size selection. If not set, the "
                                                                    "default value = fd (Freedman "
                                                                    "Diaconis Estimator). The list of the potential "
    "options are listed in: \n"
    "https://docs.scipy.org/doc/numpy/reference/generated/numpy.histogram_bin_edges.html#numpy.histogram_bin_edges")

    # If no arguments are given, use input file
    if len(sys.argv) == 1:
        algae_type, start, varbs, ens, end, plot, monthly, grid, sample, mask, output, covary, hist = file_entry()
    elif len(sys.argv) == 2 and (sys.argv[1] == '-ex' or sys.argv[1] == '--example'):
        algae_type, start, varbs, ens, end, plot, monthly, grid, sample, mask, output, covary, hist = file_entry(example=True)
    else:
        # Arguments
        args = parser.parse_args()
        algae_type = args.prefix
        start = args.start_date
        varbs = args.vars
        ens = args.ens[0]
        end = args.end_date
        plot = args.plot
        monthly = args.monthly
        grid = args.grid
        sample = args.sample
        mask = args.mask
        output = args.output
        covary = args.covary
        hist = args.hist

    # Get command line arguments
    argv = 'python main.py ' + algae_type + ' ' + str(start)
    if end:
        argv = argv + ' ' + end
    av = ' '.join(varbs)
    argv = argv + ' -v ' + av + ' -e ' + str(ens)

    # Get split start date
    day_s, mon_s, yr_s = get_date(start)
    if not end:  # If end date not given, use the end of start year
        if StartBools.just_start_year:
            end = str(yr_s)
        elif StartBools.just_start_year_month:
            end = str(yr_s) + "-" + str(mon_s)

    # Get split end date
    day_e, mon_e, yr_e = get_date(end, start=False)

    # Print user input
    print("Arguments:")
    print("- algae type: ", algae_type)
    print("- variables: ", varbs)
    print("- start date: " + str(yr_s) + "-" + str(mon_s) + "-" + str(day_s))
    print("- end date: " + str(yr_e) + "-" + str(mon_e) + "-" + str(day_e))

    # Check that dates are in valid order
    is_valid = check_valid_order([day_s, mon_s, yr_s], [day_e, mon_e, yr_e])
    if not is_valid:
        print("Error: Invalid start and end date")
        print("  - The end date is earlier than the start date")
        sys.exit()
    print("Number of ensembles:", ens)
    if plot:
        print("Plotting option selected.")
        argv = argv + ' -p'
    if monthly:
        print("Monthly date expected.")
        argv = argv + ' -m'

    lat, lon = None, None
    if grid:
        lat, lon = grid[0], grid[1]
        print("Grid point option selected.")
        argv = argv + ' -g ' + str(grid[0]) + ' ' + str(grid[1])
    if sample:
        lat, lon = sample[0], sample[1]
        print("Sample point option selected.")
        argv = argv + ' -s ' + str(sample[0]) + ' ' + str(sample[1])

    if mask:
        if isinstance(mask, list):
            mask = mask[0]
        print("Masking grid option selected.")
        argv = argv + ' -mk ' + mask
    if output:
        print("Save analysis data output selected.")
        argv = argv + ' -o'
    if covary:
        print("Co-varying option selected.")
        argv = argv + ' -cv'

    if not hist:
        hist = 'fd'
    elif hist:
        argv = argv + ' -h ' + hist

    print("Histogram bin selection option:", hist)

    # Call functions to perform analysis
    start = [day_s, mon_s, yr_s]
    end = [day_e, mon_e, yr_e]

    # Extract data from files
    saved, units, files, nan_values = extract_data(algae_type, varbs, start, end, ens,
                                                       monthly=monthly, lat=lat, lon=lon, grid=grid,
                                                       mask=mask)

    # Compute averages
    ens_means = compute_average(saved, nan_values)

    # Write to netcdf file, pass in means and variable names
    write_means_to_netcdf_file(files, ens_means, varbs, start, end, argv, test=True)

    # create_histogram(saved, units, start, end, nan_values, sel=hist, plot=plot)
