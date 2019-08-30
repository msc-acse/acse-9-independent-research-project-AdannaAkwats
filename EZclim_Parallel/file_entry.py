"""
Created by Adanna Akwataghibe (Github: AdannaAkwats)
"""
import os
import directories
from utils import *


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
    spatial = args[6]
    total = args[7]
    plot = args[8]
    monthly = args[9]
    grid = args[10]
    sample = args[11]
    mask = args[12]
    output = args[13]
    covary = args[14]
    hist = args[15]
    lb = args[16]
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
    if not hist:
        hist = ['fd']
    else:
        hist = hist.split(',')
        if len(hist) == 2:
            try:
                hist[0] = hist[0].strip()
                hist[1] = hist[1].strip()
            except Exception:
                print("ERROR in function file_entry: Histogram argument in input file may be missing a comma.")
                print(" Note that for a 2D histogram, if two bin inputs are given, both must be integers.")
                sys.exit()
        elif len(hist) > 2:
            print("ERROR in function file_entry: Number of arguments given for histogram is invalid. ")
            sys.exit()

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
    # spatial, total, monthly, output, covary, save_ext, calc_areas
    if not isinstance(spatial, bool):
        print("ERROR in function file_entry: Spatial argument in input file is invalid. "
              "It must be set to True/False or left empty (False).")
        sys.exit()
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
    if not isinstance(calc_areas, bool):
        print("ERROR in function file_entry: Calculate areas argument in input file is invalid. "
              "It must be set to True/False or left empty (False).")
        sys.exit()

    return algae_type, start, varbs, ens, end, analysis, spatial, total, plot, monthly, grid, sample, mask, output, covary, hist, lb, func, calc_areas, index
