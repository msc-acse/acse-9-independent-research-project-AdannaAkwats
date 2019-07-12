from nco import Nco
import time
from utils import *
from netCDF4 import Dataset


def write_means_to_netcdf_file(ens_files, abs_files, ens_means, variables, start_date, end_date, argv, test=False):
    """
    Write means computed in netcdf files
    :param ens_files: initial files arranged in ensemble order
    :param abs_files: absolute path of ens_files
    :param ens_means: ensemble means calculated calling function compute_average
    :param variables: list of variables
    :param start_date: start date list in [day, month, year] format
    :param end_date: end date list in [day, month, year] format
    :param argv: string containing command line arguments used
    :param test: if test is true, make some changes specific to files on my pc
    :return: None, files created in folder analysis/ensemble_means
    """

    # Initialise Nco
    nco = Nco()

    # Get data path
    data_path = directories.CLIMATE_DATA

    # Start and end date string
    start_end_str = str(start_date[0]) + "-" + str(start_date[1]) + "-" + str(start_date[2]) + " and " + \
                    str(end_date[0]) + "-" +str(end_date[1]) + "-" + str(end_date[2])

    # Get folder to store ensemble means
    results = directories.ANALYSIS
    mean_folder = os.path.abspath(os.path.join(results, directories.MEANS))
    if test:
        mean_folder = mean_folder.replace("Adanna Akwataghibe", "Adanna")

    # Go through ensembles, merge files to get output and write to output
    for i in range(len(ens_means)):
        # Get first file name in specific ensemble and add last year to name - use as output file name
        output_file = ""
        if ens_files[i][0].endswith(".nc"):
            # append start and end date to file name
            output_file = ens_files[i][0][:-3] + '_' + str(end_date[2]) + '.nc'
            output_file = os.path.basename(os.path.normpath(output_file))

        output_file = os.path.join(mean_folder, output_file)

        # Merge files in ensemble in output_file
        nco.ncecat(input=abs_files[i], output=output_file)

        # Write means to file
        with Dataset(abs_files[i][0], 'r') as src, Dataset(output_file, 'a') as dest:
            for var in variables:
                # create dataset identical to original variable in file
                mean_var_name = var + '_mean'
                datatype = src.variables[var].datatype
                # Get dimensions without time
                dims = src.variables[var].dimensions[1:]
                mean_var = dest.createVariable(mean_var_name, datatype, dims)
                # save means in variable
                mean_var[:] = ens_means[i][var]
                mean_var.setncatts(src[var].__dict__)
                mean_var.long_name = mean_var.long_name + ' averaged between ' + start_end_str

            # Write to description and history of file
            desc_str = "Added averages of variables " + ', '.join(variables) + " within time period " + \
                                       start_end_str

            if 'description' in dest.ncattrs():
                dest.description = desc_str + ' \n' + dest.description
            else:
                dest.description = desc_str
            dest.history = time.ctime(time.time()) + ': Commands used to produce file: ' + argv + ' \n' + \
                               time.ctime(time.time()) + ': Functions used:  extract_data, compute_average,' \
                                                         ' write_means_to_netcdf_file' + ' \n' + dest.history

    print("Mean ensemble files created in " + os.path.join(directories.ANALYSIS, directories.MEANS) + " folder.")


