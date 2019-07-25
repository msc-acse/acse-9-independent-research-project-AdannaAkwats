from nco import Nco
import time
from utils import *
from netCDF4 import Dataset
import pickle


def write_means_to_netcdf_file(ens_files, abs_files, ens_means, analysis_str, variables, start_date, end_date, argv, test=False):
    """
    Write means computed in netcdf files
    :param ens_files: initial files arranged in ensemble order
    :param abs_files: absolute path of ens_files
    :param ens_means: ensemble means calculated calling function compute_average
    :param analysis_str: type of analysis computed in ens_means: 'mean', 'std', 'median', 'rms' or 'all' of them
    :param variables: list of variables
    :param start_date: start date list in [day, month, year] format
    :param end_date: end date list in [day, month, year] format
    :param argv: string containing command line arguments used
    :param test: if test is true, make some changes specific to files on my pc
    :return: None, files created in folder analysis/ensemble_means
    """
    # Assertions
    assert ens_files is not None and abs_files is not None and variables is not None
    assert check_list_date(start_date) and check_list_date(end_date)

    # For all stats
    analysis_strs = ['mean', 'std', 'median', 'rms']

    # Initialise Nco
    nco = Nco()

    # Get data path
    data_path = directories.CLIMATE_DATA

    # Start and end date string
    start_end_str = str(start_date[2]) + "-" + str(start_date[1]) + "-" + str(start_date[0]) + " and " + \
                    str(end_date[2]) + "-" +str(end_date[1]) + "-" + str(end_date[0])

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
                loop = 1
                if analysis_str == 'all':
                    loop = len(analysis_strs)
                for a in range(loop):
                    # create dataset identical to original variable in file
                    mean_var_name = var + '_' + analysis_str
                    if analysis_str == 'all':
                        mean_var_name = var + '_' + analysis_strs[a]
                    datatype = src.variables[var].datatype
                    # Get dimensions without time
                    dims = src.variables[var].dimensions[1:]
                    mean_var = dest.createVariable(mean_var_name, datatype, dims)
                    # save means in variable
                    if analysis_str == 'all':
                        mean_var[:] = ens_means[i][a][var]
                    else:
                        mean_var[:] = ens_means[i][var]
                    mean_var.setncatts(src[var].__dict__)
                    mean_var.long_name = mean_var.long_name + ' averaged between ' + start_end_str

            # Write to description and history of file
            desc_str = "Added " + analysis_str + " of variables " + ', '.join(variables) + " within time period " + \
                                       start_end_str
            if analysis_str == 'all':
                desc_str = "Added " + "mean, standard deviation, RMS and median " + " of variables " + ', '.join(variables) + " within time period " + \
                           start_end_str

            if 'description' in dest.ncattrs():
                dest.description = desc_str + ' \n' + dest.description
            else:
                dest.description = desc_str
            dest.history = time.ctime(time.time()) + ': Commands used to produce file: ' + argv + ' \n' + \
                               time.ctime(time.time()) + ': Functions used:  extract_data, compute_stats_analysis,' \
                                                         ' write_means_to_netcdf_file' + ' \n' + dest.history

    print("Ensemble files created in " + os.path.join(directories.ANALYSIS, directories.MEANS) + " folder.")


def save_extract_data_to_file(list_ens, full_list_ens, ens_files, abs_files, args_dict):
    """
    Save given objects to pickle file in current directory
    :param list_ens: List of ensemble dicts
    :param full_list_ens: List of ensemble dicts of full map (only if lat and lon are given)
    :param ens_files: List of file used
    :param abs_files: Absolute path of ens_files
    :param args_dict: dictionary of all arguments in command line
    """
    afile = open(r'extracted_data.pkl', 'wb')
    pickle.dump(list_ens, afile)
    pickle.dump(ens_files, afile)
    pickle.dump(abs_files, afile)
    pickle.dump(args_dict, afile)
    if args_dict['lat'] is not None:
        pickle.dump(full_list_ens, afile)
    afile.close()


def load_extract_data_from_file(filename):
    """
    Load objects from pickle file
    :param filename: name of pickle file
    :return: - list of ensemble (iris cube) data
             - list of ensemble files
             - list of ensemble files in absolute path
             - dictionary of arguments in command line
             - list of ensemble data (full map) - this is not None if lat/lon are not None
    """
    afile = open(filename, "rb")
    list_ens = pickle.load(afile)
    ens_files = pickle.load(afile)
    abs_files = pickle.load(afile)
    args_dict = pickle.load(afile)
    full_list_ens = None
    if args_dict['lat'] is not None:
        full_list_ens = pickle.load(afile)
    afile.close()

    return list_ens, ens_files, abs_files, args_dict, full_list_ens



