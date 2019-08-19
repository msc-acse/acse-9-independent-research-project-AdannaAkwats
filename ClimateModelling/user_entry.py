import argparse
from calendar import monthrange
from Extract import *
from Analysis import *
from WriteOutput import *
from plots import *
import directories
from utils import check_valid_order, check_analysis, check_variables_covary, print_end_statement, check_valid_indices
from calculate_indices import *

# Booleans set when user gives a just a year (or a year and month)
class StartBools:
    just_start_year = False
    just_start_year_month = False


def get_date(date_entry, start=True):
    """
    Separate date  d-m-y into day, month and year
    :param date_entry: string containing date e.g. 2020-04
    :param start: if set, then it is the start date of the analysis, else it is the end date
    :return: day, month, year (all integers)
    """

    try:
        date_ = map(int, date_entry.split('-'))
        date_list = list(date_)
    except ValueError:
        print("ERROR in function get_date(): Date written in unrecognisable format. Please try again.")
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
        print("ERROR in function get_date(): too many split arguments")

    # check that these are valid dates
    try:
        datetime(year, month, day)
    except ValueError:
        print("ERROR in function get_date(): invalid date")
        return None

    return day, month, year


def file_entry(example=False):
    """
    Get arguments from input or input_example file
    :param example: set if input_example is used
    :return: arguments
    """
    filename = os.path.join(directories.INPUT, directories.INPUT_FILE)
    if example:
        filename = os.path.join(directories.INPUT, directories.INPUT_EXAMPLE_FILE)

    # Save arguments
    args = []
    # open a file using with statement
    with open(filename, 'r') as fh:
        for curline in fh:
            # check if the current line
            # starts with "#"
            if "#" not in curline and len(curline) > 1:
                arg = curline.split(':')[1].strip()
                args.append(arg)

    # Check if any required arguments are not filled in
    for i in range(4):
        if len(args[i]) == 0:
            print("ERROR in function file_entry: the input file has missing required arguments.")
            sys.exit()

    # Check if any optional arguments are not filled in
    for i in range(4, 20):
        if len(args[i]) == 0:
            args[i] = False
        elif args[i].lower() == 'false' or args[i].lower() == 'f':
            args[i] = False
        elif args[i].lower() == 'true' or args[i].lower() == 't':
            args[i] = True

    # Liat of all arguments
    # Required
    # Preix of file
    algae_type = args[0]
    algae_type = algae_type.split(',')
    if len(algae_type) == 2:
        algae_type[0] = algae_type[0].strip()
        algae_type[1] = algae_type[1].strip()
    elif len(algae_type) > 2:
        print("ERROR in function file_entry: Too many arguments in 'Prefix' in input file.")
        sys.exit()
    elif len(algae_type) < 1:
        print("ERROR in function file_entry: Required argument 'Prefix' is missing in input file.")
        sys.exit()

    # Start date
    start = args[1]
    start = start.split(',')
    if len(start) == 2:
        start[0] = start[0].strip()
        start[1] = start[1].strip()
    elif len(start) > 2:
        print("ERROR in function file_entry: Too many arguments in start date in input file.")
        sys.exit()
    elif len(start) < 1:
        print("ERROR in function file_entry: Required argument 'Start date' is missing in input file.")
        sys.exit()

    # Ensemsble number
    ens = int(args[3])
    if not args[3]:
        print("ERROR in function file_entry: Required argument 'Number of ensembles' is missing in input file.")
        sys.exit()

    # Optional
    end = args[4]
    analysis = args[5]
    total = args[6]
    plot = args[7]
    monthly = args[8]
    grid = args[9]
    sample = args[10]
    mask = args[11]
    output = args[12]
    covary = args[13]
    hist = args[14]
    lb = args[15]
    save_ext = args[16]
    func = args[17]
    calc_areas = args[18]
    index = args[19]

    # Make sure ens > 0
    if ens < 1:
        print("ERROR in function file_entry: Number of ensembles should be > 0.")
        sys.exit()

    # Check end dates
    if end:
        end = end.split(',')
        if len(end) == 2:
            end[0] = end[0].strip()
            end[1] = end[1].strip()
        elif len(end) > 2:
            print("ERROR in function file_entry: Too many arguments in end date in input file.")
            sys.exit()


    # Check for analysis
    if analysis:
        analysis = analysis.split(',')
        for i in range(len(analysis)):
            analysis[i] = analysis[i].strip()

    # Check for plot
    if plot:
        plot = [int(plot.strip())]

    # Check for grid and sample
    if grid:
        grid = grid.split(',')
        # If not split, then it is a file
        if len(grid) == 1:
            grid[0] = grid[0].strip()
            # Check if txt or nc file or linear or rotate
            check_sample_grid_one_arg(grid, 'file_entry')
        elif len(grid) == 2:
            try:
                grid[0] = float(grid[0].strip())
                grid[1] = float(grid[1].strip())
            except Exception:
                print("ERROR in function file_entry: Grid argument in input file may be missing a comma.")
                sys.exit()
        else:
            print("ERROR in function file_entry: Number of arguments given for grid point is invalid. ")
            print("   - Latitude and Longitude pair or file name expected.")
            sys.exit()
    elif sample:
        sample = sample.split(',')
        # If not split, then it is a file
        if len(sample) == 1:
            sample[0] = sample[0].strip()
            # Check if txt or nc file or linear or rotate
            check_sample_grid_one_arg(sample, 'file_entry')
        elif len(sample) == 2:
            try:
                sample[0] = float(sample[0].strip())
                sample[1] = float(sample[1].strip())
            except Exception:
                print("ERROR in function file_entry: Sample argument in input file may be missing a comma.")
                sys.exit()
        else:
            print("ERROR in function file_entry: Number of arguments given for sample point is invalid. ")
            print("   - Latitude and Longitude pair or file name expected.")
            sys.exit()

    # Histogram option
    if len(hist) == 0 or not hist:
        hist = 'fd'

    # Check for longitude centre
    if lb:
        try:
            lb = [float(lb.strip())]
        except Exception:
            print("ERROR in function file_entry: Longitude centre may not be a float.")
            sys.exit()

    # Get variables and put in list
    varbs = args[2].split(',')
    for i in range(len(varbs)):
        varbs[i] = varbs[i].strip()

    # Check for user function
    if func:
        func = func.split(',')
        try:
            func[0] = func[0].strip()
            func[1] = func[1].strip()
        except Exception:
            print("ERROR in function file_entry: User function argument in input file may be missing a comma.")
            sys.exit()

    # Check total, if ensemble is just one, then turn it off
    if ens == 1 and total:
        total = False

    # Check indices
    if index:
        index = index.strip()
        check_valid_indices(index)

    # Check boolean options
    # total, monthly, output, covary, save_ext, calc_areas
    if not isinstance(total, bool):
        print("ERROR in function file_entry: Total ensemble stats argument in input file is invalid. "
              "It must be set to True/False or left empty (False).")
        sys.exit()
    if not isinstance(monthly, bool):
        print("ERROR in function file_entry: Monthly argument in input file is invalid. "
              "It must be set to True/False or left empty (False).")
        sys.exit()
    if not isinstance(output, bool):
        print("ERROR in function file_entry: Save Output argument in input file is invalid. "
              "It must be set to True/False or left empty (False).")
        sys.exit()
    if not isinstance(covary, bool):
        print("ERROR in function file_entry: Covary argument in input file is invalid. "
              "It must be set to True/False or left empty (False).")
        sys.exit()
    if not isinstance(save_ext, bool):
        print("ERROR in function file_entry: Save extract data argument in input file is invalid. "
              "It must be set to True/False or left empty (False).")
        sys.exit()
    if not isinstance(calc_areas, bool):
        print("ERROR in function file_entry: Calculate areas argument in input file is invalid. "
              "It must be set to True/False or left empty (False).")
        sys.exit()

    return algae_type, start, varbs, ens, end, analysis, total, plot, monthly, grid, sample, mask, output, covary, hist, lb, save_ext, func, calc_areas, index


def user_entry():
    """
    Get user input from command line or from input file and run full program.
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
    parser.add_argument('-pf', '--prefix', nargs='+', required=True, help="<Required> This is the prefix of the file - in the filenames format section, this is the START OF FILENAME.")
    parser.add_argument('start_date', nargs='+', help="""Start date of analysis 
    Can be in the following formats:
    ----------------------------------
    YYYY-MM-DD : e.g. 2020-04-12
    YYYY-MM    : e.g. 2020-04
    YYYY       : e.g. 2020 
    - If day is not given, the 1st of the given month will be used i.e 2020-04 => 2020-04-01
    - If day and month is not given, 1st Jan will be used as the start date i.e 2020 => 2020-01-01""")
    parser.add_argument('end_date', nargs='*', help=""" <Not required> End date of analysis - format is the same as start_date
    -----------------------------------end_date not given-------------------------------------
    - If only start year is given, the end_date is automatically set to the 31 Dec of start year
    - If start year and month is given, then end_date is set to the end of the start month
       -----------------------------------end_date given-------------------------------------
    - If day is not given, the end of the given month will be used i.e 2020-04 => 2020-04-30
    - If day and month is not given, 31 Dec will be used as the end date i.e 2020 => 2020-12-31""")
    parser.add_argument('-v', '--vars', nargs='+', metavar="variables", help="<Required> Variables of data to analyse",
                        required=True)
    parser.add_argument('-p', '--plot', nargs=1, metavar=("ensemble_number"), help="""Plot map, histogram and timeseries graphs
    E.g. --plot 1
    The ensemble to plot must be included. """)
    parser.add_argument('-m', '--monthly', action="store_true", help="Data in file is stored in monthly increments.")
    group = parser.add_mutually_exclusive_group()
    group.add_argument('-g', '--grid', nargs='+', type=float, metavar=("(lat, lon) or filename or linear/rotate"),
                       help="""
                       Grid Point: Latitude, Longitude
                        Uses grid point that latitude and longitude lies in.
                        Other commands:
                        - You can define a list of grid points in a .txt file e.g check INPUT/sample_points.txt
                           - Grid Point: sample_points.txt
                        - You can regrid to a grid (using nearest neighbour interpolation) defined in a NETCDF file:
                           - Grid Point: example_file.nc
                        Cannot be used in conjunction with sample point.
                       """)
    group.add_argument('-s', '--sample', nargs='+', type=float, metavar=("(lat, lon) or filename or linear/rotate"),
                       help="""
                       Sample Point: Latitude, Longitude
                        Uses sample point given by latitude and longitude using interpolation.
                        Other commands:
                        - You can define a list of sample points in a .txt file e.g check INPUT/sample_points.txt
                           - Sample Point: sample_points.txt
                        - You can regrid to a grid (using linear interpolation) defined in a NETCDF file:
                           - Sample Point: example_file.nc
                        Cannot be used in conjunction with grid point.
                       """)
    group.add_argument('-lc', '--lon_centre', nargs=1, type=float, help="Longitude to centre map on.")
    parser.add_argument('-mk', '--mask', nargs=1, metavar="filename", help="Uses masking grid given as a file "
                                                                           "(contains boolean array to be imposed on "
                                                                           "the global grid).")
    parser.add_argument('-o', '--output', action="store_true", help="If plot option selected, save data output of histogram and timeseries "
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
    parser.add_argument('-se', '--save_extract', action="store_true", help="Save extracted (iris.cube) data in pkl file.")
    parser.add_argument('-u', '--user', nargs=2, metavar=('file_name', 'function_name'),
                        help="""Use function written by the user and stored in user_function folder for analysis. 
                        file_name : name of file that contains function in user_function folder
                        function_name : name of function to call 
                        Note: user functions are expected to only take in a cube as an argument. An example of a function 
                        can be found in user_function/example_function.py
                        """)
    parser.add_argument('-a', '--analysis', nargs='+', help="""Analysis performed on data set.
    If not specified, then all analysis listed below will be performed.
    Types of analysis:
    - mean
    - std (Standard deviation)
    - rms (Root mean squared error)
    - median
    You can also select a combination of analysis to perform e.g. -a mean rms """)
    parser.add_argument('-ca', '--areas', action="store_true", help="Calculate areas of grid boxes of latitude and"
                                                                    " longitude and saves to NetCDF file areas.nc in results folder")
    parser.add_argument('-t', '--total', action="store_true",
                        help="""Total ensemble stats: True/False : The analysis will be performed over the whole ensemble given.
                        - If set True, all the ensembles will be averaged as a collection.
                        - If set False, the ensembles will be averaged individually.""")
    parser.add_argument('-i', '--index', metavar=('index'),
                        help="""Calculate index given 
                            The control run is the FIRST file prefix set and the corresponding start/end date. 
                            The future run is the SECOND file prefix set and the corresponding second start/end date
                            Types of inidices that can be calculated:          
                            enso : The Oceanic Niño Index (ONI) 
                            nino12 : Niño 1+2 Index
                            nino4 : Niño 4 Index
                            tni : The Trans-Niño Index (TNI)
                            iod : Indian Ocean Dipole (IOD) Mode Index 
                            amo : Atlantic Multidecadal Oscillation (AMO) Index
                            pdo : Pacific Decadal Oscillation (PDO) Index 
                            ao : Arctic Oscillation (AO; Northern Annular Mode) Index 
                            aao : Antarctic Oscillation (AAO; Southern Annular Mode) Index 
                            nao : North Atlantic Oscillation (NAO) Index
                            """)
    # Log output
    old_stdout = sys.stdout
    log_file = open("message.log", "w")
    sys.stdout = log_file

    # Initialise the variables
    algae_type, varbs, start, end, ens, monthly, lat, lon, grid, sample, mask, output, covary, hist, plot, \
        lon_centre, save_ext, func = None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None
    analysis, calc_areas, total, points_sample_grid = None, None, None, None
    argv = None
    loaded_data, second_date_given, start2, end2 = None, False, None, None
    saved, ens_files, abs_files, full_saved, dim_renames = None, None, None, None, None
    args_dict = {}

    # If no arguments are given, use input file
    if len(sys.argv) == 1:
        algae_type, start, varbs, ens, end, analysis, total, plot, monthly, grid, sample, mask, output, covary, hist, lon_centre, save_ext, func, calc_areas, index = file_entry()
    elif len(sys.argv) == 2 and (sys.argv[1] == '-ex' or sys.argv[1] == '--example'):
        algae_type, start, varbs, ens, end, analysis, total, plot, monthly, grid, sample, mask, output, covary, hist, lon_centre, save_ext, func, calc_areas, index = file_entry(example=True)
    elif len(sys.argv) == 2 and sys.argv[1][-3:] == 'pkl':  # Get pickle file to open
        loaded_data = sys.argv[1]
        # Check that it is a pickle file
        if loaded_data is None or loaded_data[-3:] != 'pkl':
            print("ERROR in function user_entry : Pickle file must be used.")
            sys.exit()
    else:
        # Arguments
        args = parser.parse_args()

        algae_type = args.prefix
        start = args.start_date
        varbs = args.vars
        ens = args.ens[0]
        end = args.end_date
        analysis = args.analysis
        total = args.total
        plot = args.plot
        monthly = args.monthly
        grid = args.grid
        sample = args.sample
        mask = args.mask
        output = args.output
        covary = args.covary
        hist = args.hist
        lon_centre = args.lon_centre
        save_ext = args.save_extract
        func = args.user
        calc_areas = args.areas
        index = args.index

    # If not loading data from pickle file, then get from command line
    if loaded_data is None:
        # Get command line arguments
        argv = 'python main.py ' + algae_type[0]
        if len(algae_type) == 2:
            argv = ' ' + algae_type[1]
        argv = ' ' + start[0]
        if len(start) == 2:
            argv = ' ' +  start[1]
        if end:
            argv = argv + ' ' + end[0]
            if len(end) == 2:
                argv = argv + ' ' + end[1]
        av = ' '.join(varbs)
        argv = argv + ' -v ' + av + ' -e ' + str(ens)

        if end and len(start) < len(end):
            print("ERROR in function user_entry: Start dates are required.")
            sys.exit()

        if len(algae_type) > 2:
            print("ERROR in function user_entry: Too many arguemnts given for 'Prefix' argument.")
            sys.exit()

        # All dates
        day_s, mon_s, yr_s, day_e, mon_e, yr_e = None, None, None, None, None, None
        day_s2, mon_s2, yr_s2, day_e2, mon_e2, yr_e2 = None, None, None, None, None, None
        # Get split start date
        if len(start) == 1:
            day_s, mon_s, yr_s = get_date(start[0])
            if not end:  # If end date not given, use the end of start year
                if StartBools.just_start_year:
                    end = str(yr_s)
                elif StartBools.just_start_year_month:
                    end = str(yr_s) + "-" + str(mon_s)
            else:
                end = end[0]
            # Get split end date
            day_e, mon_e, yr_e = get_date(end, start=False)

        # If extra year is given
        if len(start) == 2:
            second_date_given = True
            StartBools.just_start_year, StartBools.just_start_year_month = False, False
            day_s, mon_s, yr_s = get_date(start[0])
            if not end:
                if StartBools.just_start_year:
                    fst_end = str(yr_s)
                elif StartBools.just_start_year_month:
                    fst_end = str(yr_s) + "-" + str(mon_s)
            else:
                fst_end = end[0]
            # Get split end date
            day_e, mon_e, yr_e = get_date(fst_end, start=False)
            # Get next start and end date
            day_s2, mon_s2, yr_s2 = get_date(start[1])
            if not end or (len(end) == 1):
                if StartBools.just_start_year:
                    end = str(yr_s2)
                elif StartBools.just_start_year_month:
                    end = str(yr_s2) + "-" + str(mon_s2)
            else:
                end = end[1]
            # Get split end date
            day_e2, mon_e2, yr_e2 = get_date(end, start=False)
        elif len(start) > 2:
            print("ERROR in function user_entry: Too many arguemnts given for 'Start date' argument.")
            sys.exit()

        # Print user input
        print("Arguments:")
        if len(algae_type) == 1:
            print("- file prefix: ", algae_type[0])
        if len(algae_type) == 2:
            print("- first file prefix: ", algae_type[0])
            print("- second file prefix: ", algae_type[1])
        print("- variables: ", varbs)
        print("- start date: " + str(yr_s) + "-" + str(mon_s) + "-" + str(day_s))
        print("- end date: " + str(yr_e) + "-" + str(mon_e) + "-" + str(day_e))
        if second_date_given:
            print("- second start date: " + str(yr_s2) + "-" + str(mon_s2) + "-" + str(day_s2))
            print("- second end date: " + str(yr_e2) + "-" + str(mon_e2) + "-" + str(day_e2))

        # Check that dates are in valid order
        is_valid = check_valid_order([day_s, mon_s, yr_s], [day_e, mon_e, yr_e])
        if not is_valid:
            print("ERROR in function user_entry: Invalid start and end date")
            print("  - The end date is earlier than the start date")
            sys.exit()
        if second_date_given:
            is_valid = check_valid_order([day_s2, mon_s2, yr_s2], [day_e2, mon_e2, yr_e2])
            if not is_valid:
                print("ERROR in function user_entry: Invalid second start and second end date")
                print("  - The end date is earlier than the start date")
                sys.exit()
        print("Number of ensembles:", ens)

        if analysis:
            print("Analysis: ", analysis)
            a_ = ' '.join(analysis)
            argv = argv + ' -a ' + a_
            check_analysis(analysis)
        if total:
            print("Total ensemble stats option selected.")
            argv = argv + ' -t'
        if plot:
            print("Plotting option selected.")
            argv = argv + ' -p ' + str(plot[0])
        else:
            plot = None
        if monthly:
            print("Monthly date expected.")
            argv = argv + ' -m'

        if grid:
            if len(grid) == 2:
                lat, lon = grid[0], grid[1]
                print("Grid point option selected.")
                argv = argv + ' -g ' + str(grid[0]) + ' ' + str(grid[1])
            elif len(grid) == 1:
                # Check if txt or nc file or linear or rotate
                check_sample_grid_one_arg(grid, 'user_entry')
                points_sample_grid = grid[0]
                print("Grid point option selected.")
                argv = argv + ' -g ' + str(grid[0])
            else:
                print("ERROR in function user_entry: Grid point argument has invalid number of arguments.")
                sys.exit()
        elif sample:
            if len(sample) == 2:
                lat, lon = sample[0], sample[1]
                print("Sample point option selected.")
                argv = argv + ' -s ' + str(sample[0]) + ' ' + str(sample[1])
            elif len(sample) == 1:
                # Check if txt or nc file or linear or rotate
                check_sample_grid_one_arg(sample, 'user_entry')
                points_sample_grid = sample[0]
                print("Sample point option selected.")
                argv = argv + ' -s ' + str(sample[0])
            else:
                print("ERROR in function user_entry: Sample point argument has invalid number of arguments.")
                sys.exit()

        if mask:
            if isinstance(mask, list):
                mask = mask[0]
            print("Masking grid option selected.")
            argv = argv + ' -mk ' + mask
        elif not mask:
            mask = None
        if output:
            print("Save analysis data output selected.")
            argv = argv + ' -o'
        if covary:
            print("Co-varying option selected.")
            argv = argv + ' -cv'
            check_variables_covary(varbs)

        if not hist:
            hist = 'fd'
        elif hist:
            argv = argv + ' -h ' + hist

        print("Histogram bin selection option:", hist)

        if func:
            print("User function given: " + str(func[0]) + ", " + str(func[1]))
            argv = argv + ' -u ' + func[0] + ' ' + func[1]

        if calc_areas:
            print("Calculate areas option selected.")
            argv = argv + ' -ca'

        # Check index is given with second date
        if index and not second_date_given:
            print("ERROR in function user_entry: Index must be given with a second start date set.")
            sys.exit()

        if index:
            print("Index option selected: " + index)
            argv = argv + ' -i'

        if lon_centre:
            lon_centre = lon_centre[0]
            print("Longitude centering option selected.")
            argv = argv + ' -lc ' + str(lon_centre)
        elif not lon_centre:
            lon_centre = None

        # Call functions to perform analysis
        start = [day_s, mon_s, yr_s]
        end = [day_e, mon_e, yr_e]
        if second_date_given:
            start2 = [day_s2, mon_s2, yr_s2]
            end2 = [day_e2, mon_e2, yr_e2]


        # Calculate indices
        if index:  # Self contained action
            calculate_index(algae_type, index, varbs, start, end, start2, end2, monthly=monthly, test=True)
            sys.exit()

        # Extract data from files
        extract = Extract(algae_type[0], varbs, start, end, ens, monthly=monthly, lat=lat, lon=lon, grid=grid,
                                                               points_sample_grid=points_sample_grid,
                                                               lon_centre=lon_centre, maskfile=mask,
                                                               calc_areas=calc_areas)
        saved, ens_files, abs_files, full_saved = extract.extract_data()

        if second_date_given:
            at = None
            if len(algae_type) == 2:
                at = algae_type[1]
            else:
                at = algae_type[0]
            extract = Extract(at, varbs, start2, end2, ens, monthly=monthly, lat=lat, lon=lon, grid=grid,
                                                                   points_sample_grid=points_sample_grid,
                                                                   lon_centre=lon_centre, maskfile=mask,
                                                                   calc_areas=calc_areas)
            saved2, ens_files2, abs_files2, full_saved2 = extract.extract_data()

        # Put all values in dictionary
        args_dict['algae_type'] = algae_type
        args_dict['varbs'] = varbs
        args_dict['start'] = start
        args_dict['end'] = end
        args_dict['start2'] = start2
        args_dict['end2'] = end2
        args_dict['analysis'] = analysis
        args_dict['total'] = total
        args_dict['cov'] = covary
        args_dict['ens'] = ens
        args_dict['monthly'] = monthly
        args_dict['lat'] = lat
        args_dict['lon'] = lon
        args_dict['grid'] = grid
        args_dict['lon_centre'] = lon_centre
        args_dict['save_out'] = output
        args_dict['hist'] = hist
        args_dict['plot'] = plot
        args_dict['argv'] = argv
        args_dict['mask'] = mask
        args_dict['func'] = func
        args_dict['points_sample_grid'] = points_sample_grid

        # Save data to pickle file
        if save_ext:
            # Add to arguments
            argv = argv + ' -se'
            args_dict['argv'] = argv
            # Save data to file
            save_extract_data_to_file(saved, full_saved, ens_files, abs_files, args_dict)
            if second_date_given:
                save_extract_data_to_file(saved2, full_saved2, ens_files2, abs_files2, args_dict, num=2)


    else:
        saved, ens_files, abs_files, args_dict, full_saved = load_extract_data_from_file(loaded_data)

    # COMPUTE ANALYSIS
    anlys = Analysis(saved)
    # user analysis
    func, ens_stats, func_name, analysis_str, nan_indices = None, None, None, None, None
    ens_stats2 = None
    if args_dict['func']:
        func = args_dict['func']
        file_name, func_name = func[0], func[1]
        ens_stats = anlys.compute_user_analysis(file_name, func_name)
    else:
        if second_date_given:
            ens_stats, analysis_str, nan_indices = anlys.calc_stats_difference(saved2, args_dict['analysis'],
                                                                         total=args_dict['total'])
        else:
            ens_stats, analysis_str, nan_indices = anlys.compute_stats_analysis(args_dict['analysis'],
                                                                          total=args_dict['total'])


    # Warning for mask and sample/grid
    if args_dict['mask'] is not None and args_dict['lat'] is not None:
        print("WARNING: Please ensure that sample/grid point is in the masked region.")

    # PLOTTING
    plot, save_out = args_dict['plot'], args_dict['save_out']
    if plot is not None or save_out:
        plot_ens_num = int(plot[0]) if plot is not None else 1

        # Plot histogram
        create_histogram(saved, ens_stats, args_dict['start'], args_dict['end'], args_dict['varbs'], sel=args_dict['hist'],
                         save_out=args_dict['save_out'], ens_num=plot_ens_num, cov=args_dict['cov'], mask=args_dict['mask'],
                         total=args_dict['total'], analysis_str=analysis_str, nan_indices=nan_indices, plot=plot,
                         second_date_given=second_date_given, start_date2=args_dict['start2'], end_date2=args_dict['end2'])

        # Only plot timeseries and map if plot is enabled
        if plot is not None:
            # Only plot map of analysis if using analysis: mean, median, std or rms and NOT grid/sample point
            if analysis_str:
                if func is None:
                    plot_map_analysis(ens_stats, args_dict['varbs'], save_out=args_dict['save_out'], ens_num=plot_ens_num,
                             analysis_str=analysis_str, total=args_dict['total'],
                             second_date_given=second_date_given)
                else:
                    print("WARNING: Map not plotted as user function is used.")
            else:
                plot_map(saved, args_dict['varbs'], save_out=args_dict['save_out'], ens_num=plot_ens_num,
                         total=args_dict['total'], second_date_given=second_date_given)

            # Plot time series and boxplot

            create_timeseries(saved, args_dict['start'], args_dict['end'], args_dict['varbs'],
                              save_out=args_dict['save_out'], ens_num=plot_ens_num, func_name=func_name,
                              second_date_given=second_date_given)

    # WRITE ANALYSIS TO NETCDF FILE
    if save_out:
        wo = WriteOutput(ens_files, abs_files, ens_stats, analysis_str, args_dict['varbs'],
                                       args_dict['start'], args_dict['end'], args_dict['argv'], saved, full_saved,
                                       total=args_dict['total'], lon_centre=args_dict['lon_centre'],
                                        mask=args_dict['mask'], lon=args_dict['lon'], lat=args_dict['lat'],
                                        grid=args_dict['grid'], user_func=func_name,
                                      points_sample_grid=args_dict['points_sample_grid'],
                                      second_date_given=second_date_given, test=True)
        wo.write_analysis_to_netcdf_file()

    print("PROGRAM SUCCESSFUL - TERMINAL FINISHED.")
    # End logging
    sys.stdout = old_stdout
    log_file.close()

    # Print to terminal when finished
    print_end_statement()

    # Show graphs
    if plot is not None:
        plt.show()


