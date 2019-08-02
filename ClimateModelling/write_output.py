from nco import Nco
import time
from utils import *
from netCDF4 import Dataset
import pickle


def write_means_to_netcdf_file(ens_files, abs_files, ens_means, analysis_str, variables, start_date, end_date, argv,
                               test=False):
    """
    Write means computed in netcdf files
    :param ens_files: initial files arranged in ensemble order
    :param abs_files: absolute path of ens_files
    :param ens_means: ensemble means calculated calling function compute_stats_analysis
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
                    dims = tuple([j for j in src.variables[var].dimensions if j.lower() != 'time' and j.lower() != 't'])
                    mean_var = dest.createVariable(mean_var_name, datatype, dims)
                    # save means in variable
                    if analysis_str == 'all':
                        cube = ens_means[i][a][var]
                        mean_var[:] = cube.data
                    else:
                        cube = ens_means[i][var]
                        mean_var[:] = cube.data
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


def write_user_analysis_to_netcdf_file(ens_files, abs_files, user_analysis, func_name, variables, start_date, end_date, argv, test=False):
    """
    Write user analysis computed in netcdf files
    :param ens_files: initial files arranged in ensemble order
    :param abs_files: absolute path of ens_files
    :param user_analysis: ensemble analysis calculated calling function compute_user_analysis
    :param func_name: user-defined function name
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

    # Initialise Nco
    nco = Nco()

    # Start and end date string
    start_end_str = str(start_date[2]) + "-" + str(start_date[1]) + "-" + str(start_date[0]) + " and " + \
                    str(end_date[2]) + "-" + str(end_date[1]) + "-" + str(end_date[0])

    # Get folder to store ensemble means
    results = directories.ANALYSIS
    mean_folder = os.path.abspath(os.path.join(results, directories.MEANS))
    if test:
        mean_folder = mean_folder.replace("Adanna Akwataghibe", "Adanna")

    # Go through ensembles, merge files to get output and write to output
    for i in range(len(user_analysis)):
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
                # Save user analysis if given else use means
                if user_analysis:
                    # Get cube of variable
                    cube = user_analysis[i][var]
                    # Get dimensions of cube and file
                    coords = cube.dim_coords
                    # Get coords names
                    coords_names = [j.name() for j in coords]
                    dims, dims_var = src.dimensions, src.variables[var].dimensions

                    # Get datatype of variable
                    datatype = src.variables[var].datatype

                    # Save new dimensions
                    new_dims = []
                    for c in range(len(coords_names)):
                        dim = dims[coords_names[c]]
                        dim_name = None
                        # If shape of dimensions differ, then create new dimension with new shape
                        if len(coords[c].points) != len(dims[coords_names[c]]):
                            dim_name = coords_names[c] + '_' + 'new'
                            dest.createDimension(dim_name, len(coords[c].points))
                        new_dims.append(dim_name)

                    # Create variable and put data
                    new_dims = tuple(new_dims)
                    user_var_name = var + '_' + func_name
                    user_var = dest.createVariable(user_var_name, datatype, new_dims)
                    user_var[:] = cube.data
                    user_var.setncatts(src[var].__dict__)
                    user_var.long_name = user_var.long_name + ' between ' + start_end_str + ' using user-defined function ' + func_name

            # Write to description and history of file
            desc_str = "Added analysis using user-defined function " + func_name + " of variables " + ', '.join(variables) + " within time period " + \
                       start_end_str
            if 'description' in dest.ncattrs():
                dest.description = desc_str + ' \n' + dest.description
            else:
                dest.description = desc_str
            dest.history = time.ctime(time.time()) + ': Commands used to produce file: ' + argv + ' \n' + \
                           time.ctime(time.time()) + ': Functions used:  extract_data, compute_user_analysis,' \
                                                     ' write_user_analysis_to_netcdf_file' + ' \n' + dest.history

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

    print("function save_extract_data_from_file: Extracted data successfully saved to extracted_data.pkl.")


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

    print("function load_extract_data_from_file: Extracted data successfully loaded from " + filename + ".")

    return list_ens, ens_files, abs_files, args_dict, full_list_ens



