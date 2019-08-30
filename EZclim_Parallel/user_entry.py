"""
Created by Adanna Akwataghibe (Github: AdannaAkwats)
"""
import argparse
from calendar import monthrange
from Extract import *
from Analysis import *
from WriteOutput import *
from plots import *
from utils import check_valid_order, check_analysis, check_variables_covary, print_end_statement
from calculate_indices import *
from file_entry import file_entry
from ProgressBar import *

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
    parser.add_argument('-end', '--end_date', nargs='*', help=""" <Not required> End date of analysis - format is the same as start_date
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
    parser.add_argument('-ht', '--hist', nargs='*',
                        metavar="number_of_bins_in_histogram", help=" Options for bin size selection. If not set, the "
                                                                    "default value = fd (Freedman "
                                                                    "Diaconis Estimator). The list of the potential "
    "options are listed in: \n"
    "https://docs.scipy.org/doc/numpy/reference/generated/numpy.histogram_bin_edges.html#numpy.histogram_bin_edges")
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
    parser.add_argument('-sp', '--spatial', action="store_true", help="Calculates averages spatially.")
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
    log_file = open("output.log", "w")
    sys.stdout = log_file

    # Init progress bar
    progress = ProgressBar(description='Climate Modelling software output', n_iter=5)

    # Initialise the variables
    algae_type, start, varbs, ens, end, analysis, spatial, total = None, None, None, None, None, None, None, None
    plot, monthly, grid, sample, mask, output, covary, hist = None, None, None, None, None, None, None, None
    lon_centre, func, calc_areas, index, lat, lon, points_sample_grid = None, None, None, None, None, None, None
    second_date_given, start2, end2 = False, None, None

    # If no arguments are given, use input file
    if len(sys.argv) == 1:
        algae_type, start, varbs, ens, end, analysis, spatial, total, plot, monthly, grid, sample, mask, output, covary, hist, lon_centre, func, calc_areas, index = file_entry()
    elif len(sys.argv) == 2 and (sys.argv[1] == '-ex' or sys.argv[1] == '--example'):
        algae_type, start, varbs, ens, end, analysis, spatial, total, plot, monthly, grid, sample, mask, output, covary, hist, lon_centre, func, calc_areas, index = file_entry(example=True)
    else:
        # Arguments
        args = parser.parse_args()

        algae_type = args.prefix
        start = args.start_date
        varbs = args.vars
        ens = args.ens[0]
        end = args.end_date
        analysis = args.analysis
        spatial = args.spatial
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
        func = args.user
        calc_areas = args.areas
        index = args.index

    # Update progress after getting input from user
    progress.update()

    # Get command line arguments
    argv = 'python main.py'
    argv = argv + ' ' + start[0]
    if len(start) == 2:
        argv = argv + start[1]
    argv = argv + ' -pf ' + algae_type[0]
    if len(algae_type) == 2:
        argv = argv + ' ' + algae_type[1]
    if end:
        argv = argv + ' -end ' + end[0]
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

    if spatial and not analysis:
        print("Error in function user_entry: Spatial argument cannot be set when no analysis is selected.")
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

    # 2 end years must be given with 2 start years
    if len(start) == 2 and len(end) != 2:
        print("ERROR in function user_entry: Both end dates must be given with both start dates.")
        sys.exit()

    # If extra year is given
    if len(start) == 2:
        second_date_given = True
        # Get first start date
        StartBools.just_start_year, StartBools.just_start_year_month = False, False
        day_s, mon_s, yr_s = get_date(start[0])

        # Get first end date
        fst_end = end[0]
        day_e, mon_e, yr_e = get_date(fst_end, start=False)

        # Get next start
        day_s2, mon_s2, yr_s2 = get_date(start[1])

        # Get next end date
        end = end[1]
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
    if spatial:
        print("Spatial analysis option selected.")
        argv = argv + ' -sp'
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
        hist = ['fd']
    elif hist:
        argv = argv + ' -ht ' + str(hist[0])
        if len(hist) == 2:
            argv = argv + ' ' + str(hist[1])
        elif len(hist) > 2:
            print("ERROR in function user_entry: Histogram argument has invalid number of arguments.")
            sys.exit()

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

    # Update progress after preocessing input from user
    progress.update()


    # Calculate indices
    if index:  # Self contained action
        calculate_index(algae_type, index, varbs, start, end, start2, end2, monthly=monthly, test=True)
        # Update progress after calculating index
        progress.update()
        progress.finish()
        sys.exit()

    # EXTRACT DATA FROM FILES
    extract = Extract(algae_type[0], varbs, start, end, ens, monthly=monthly, lat=lat, lon=lon, grid=grid,
                                                           points_sample_grid=points_sample_grid,
                                                           lon_centre=lon_centre, maskfile=mask,
                                                           calc_areas=calc_areas)
    saved, ens_files, abs_files, full_saved, dim_coords = extract.extract_data()

    saved2, ens_files2, abs_files2, full_saved2 = None, None, None, None
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
        saved2, ens_files2, abs_files2, full_saved2, _ = extract.extract_data()

    # Update progress after extracting data
    progress.update()

    # COMPUTE ANALYSIS
    anlys = Analysis(saved)
    ens_stats, func_name, analysis_str, nan_indices = None, None, None, None
    spat_calcs, spat_calcs2 = None, None
    ens_stats2 = None
    if func:  # user analysis
        file_name, func_name = func[0], func[1]
        ens_stats = anlys.compute_user_analysis(file_name, func_name)
    else:
        if second_date_given:
            ens_stats, spat_calcs, spat_calcs2, analysis_str, nan_indices= anlys.calc_stats_difference(saved2, analysis, total=total,
                                                                               spatial=spatial, dim_coords=dim_coords)
        else:
            ens_stats, analysis_str, nan_indices = anlys.compute_stats_analysis(analysis, total=total,
                                                                                spatial=spatial, dim_coords=dim_coords)

    # Warning for mask and sample/grid
    if mask is not None and lat is not None:
        print("WARNING: Please ensure that sample/grid point is in the masked region.")

    # Update progress after computing analysis
    progress.update()


    # PLOTTING
    try:
        if plot is not None or output:
            plot_ens_num = int(plot[0]) if plot is not None else 1

            # Plot histogram
            create_histogram(saved, ens_stats, start, end, varbs, sel=hist, monthly=monthly,
                             save_out=output, ens_num=plot_ens_num, cov=covary, mask=mask,
                             total=total, analysis_str=analysis_str, nan_indices=nan_indices, plot=plot,
                             second_date_given=second_date_given, start_date2=start2, end_date2=end2, spatial=spatial)

            # Only plot timeseries and map if plot is enabled
            if plot is not None:
                # Only plot map of analysis if using analysis: mean, median, std or rms and NOT grid/sample point
                if analysis_str:
                    if func is None or not func:
                        plot_map_analysis(ens_stats, varbs, save_out=output, ens_num=plot_ens_num,
                                          analysis_str=analysis_str, total=total,
                                          second_date_given=second_date_given)
                    else:
                        print("WARNING: Map not plotted as user function is used.")
                else:
                    plot_map(saved, varbs, save_out=output, ens_num=plot_ens_num, total=total,
                             second_date_given=second_date_given)

            # Plot time series and boxplot
            if analysis_str:
                create_timeseries_analysis(ens_stats, start, end, varbs, analysis_str, monthly=monthly,
                                           save_out=output, ens_num=plot_ens_num,
                                           second_date_given=second_date_given, total=total, spatial=spatial,
                                           calcs=spat_calcs, calcs2=spat_calcs2, plot=plot)
            else:
                create_timeseries(saved, start, end, varbs,
                                  save_out=output, ens_num=plot_ens_num, func_name=func_name, monthly=monthly,
                                  second_date_given=second_date_given, plot=plot)
            # Update progress after plotting
            progress.update()

    except Exception as err:
        print("Exception thrown in function user_entry when plotting: " + str(err))

    # WRITE ANALYSIS TO NETCDF FILE
    if output:
        wo = WriteOutput(ens_files, abs_files, ens_stats, analysis_str, varbs,
                         start, end, argv, saved, full_saved,
                         total=total, lon_centre=lon_centre,
                         mask=mask, lon=lon, lat=lat,
                         grid=grid, user_func=func_name,
                         points_sample_grid=points_sample_grid,
                         second_date_given=second_date_given, test=True)
        wo.write_analysis_to_netcdf_file()
        # Update progress after writing output
        progress.update()

    print("PROGRAM SUCCESSFUL - TERMINAL FINISHED.")
    # End logging
    sys.stdout = old_stdout
    log_file.close()

    # Print to terminal when finished
    print_end_statement()

    progress.finish()

    # Show graphs
    if plot is not None:
        plt.show()


