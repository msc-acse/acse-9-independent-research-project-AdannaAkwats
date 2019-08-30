"""
Created by Adanna Akwataghibe (Github: AdannaAkwats)
"""
import xarray as xr
import directories
from utils import get_files_time_period
import os
import sys

"""
Self contained Python script used to calculate indices
"""

def calculate_index(file_prefixes, index, varbs, start, end, start2, end2, monthly=False, test=False):
    """
    Calculate index given by index
    :param file_prefixes: list of prefixes
    :param index: index calculated, string
    :param varbs variables, list
    :param start: start date
    :param end: end date
    :param start2: second start date
    :param end2: second end date
    :param monthly: data is stored in monthly increments (time = 12) else assumed (time = 365)
    :return: index file
    """

    # Init files, files2
    files, files2 = None, None

    # Get files and min and maximum year
    files, _,_ = get_files_time_period(file_prefixes[0], start[2], end[2])

    if len(file_prefixes) == 2:
        files2, _,_ = get_files_time_period(file_prefixes[1], start2[2], end2[2])
    else:
        files2, _,_ = get_files_time_period(file_prefixes[0], start2[2], end2[2])

    files = [os.path.abspath(os.path.join(directories.DATA, file)) for file in files]
    files2 = [os.path.abspath(os.path.join(directories.DATA, file)) for file in files2]
    if test:
        files = [f.replace("Adanna Akwataghibe", "Adanna") for f in files]
        files2 = [f.replace("Adanna Akwataghibe", "Adanna") for f in files2]

    # Create output file
    f_name = file_prefixes[0] + '_' + str(start[2]) + '_' + str(end[2]) + '.nc'
    combined_control = os.path.abspath(os.path.join(directories.ANALYSIS, f_name))

    f2_name = file_prefixes[0] + '_' + str(start2[2]) + '_' + str(end2[2]) + '.nc'
    if len(file_prefixes) == 2:
        f2_name = file_prefixes[1] + '_' + str(start2[2]) + '_' + str(end2[2]) + '.nc'
    combined_future = os.path.abspath(os.path.join(directories.ANALYSIS, f2_name))

    if test:
        combined_control = combined_control.replace("Adanna Akwataghibe", "Adanna")
        combined_future = combined_future.replace("Adanna Akwataghibe", "Adanna")


    # Combine files into one
    if len(files) > 1:
        try:
        # Combine each ensemble file into one
            times_append = xr.open_mfdataset(files)
            times_append.to_netcdf(path=combined_control, mode='w')
        except Exception as er:
            print("ERROR in function calculate_index: failed trying to combine multiple input files into one file.")
            print(er)
            print(" - Files given: " + str(files))
            sys.exit()
    else:
        combined_control = files[0]

    if len(files2) > 1:
        try:
            # Combine each ensemble file into one
            times_append = xr.open_mfdataset(files2)
            times_append.to_netcdf(path=combined_future, mode='w')
        except Exception as er:
            print("ERROR in function calculate_index: failed trying to combine multiple input files into one file.")
            print(er)
            print(" - Files given: " + str(files2))
            sys.exit()
    else:
        combined_future = files2[0]

    # Get time string
    timely = "monthly"
    if not monthly:
        timely = "daily"

    # Create output file
    output_file = os.path.abspath(os.path.join(directories.ANALYSIS, index.upper() + '.nc'))
    if test:
        output_file = output_file.replace("Adanna Akwataghibe", "Adanna")
    of = os.path.join(directories.ANALYSIS, '')

    print("CDO ./calculate_indices.sh script running")
    os.system("./calculate_indices.sh %s %s %s %s %s %s" % (index,combined_future,combined_control,timely,output_file,of))
    print("CDO script done")

    # Print to terminal when finished
    print("PROGRAM FINISHED.")