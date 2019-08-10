import time
from utils import *
from netCDF4 import Dataset, MFDataset
import pickle
import xarray as xr
from nco import Nco


def write_means_to_netcdf_file(ens_files, abs_files, ens_means, analysis_str, variables, start_date, end_date, argv,
                               saved, full_saved, total=False, lon_centre=None, mask=None, lon=None, lat=None,
                                grid=None, test=False):
    """
    Write analysis computed in netcdf files
    :param ens_files: initial files arranged in ensemble order
    :param abs_files: absolute path of ens_files
    :param ens_means: ensemble means calculated calling function compute_stats_analysis
    :param analysis_str: type of analysis computed in ens_means: 'mean', 'std', 'median', 'rms' or 'all' of them
    :param variables: list of variables
    :param start_date: start date list in [day, month, year] format
    :param end_date: end date list in [day, month, year] format
    :param argv: string containing command line arguments used
    :param saved: iris cubes containing data
    :param full_saved: iris cubes containing data if sample/grid point are used
    :param total: set if all ensembles have been calculated together, instead of separately, boolean
    :param lon_centre: longitude center , float
    :param mask: file that contains mask data, string
    :param lon: longitude, set if grid or sample point, floats
    :param lat: latitude, set if grid or sample point, floats
    :param grid: set if grid point is given
    :param test: if test is true, make some changes specific to files on my pc
    :return: None, files created in folder analysis/ensemble_means
    """
    # Assertions
    assert ens_files is not None and abs_files is not None and variables is not None and analysis_str is not None
    assert check_list_date(start_date) and check_list_date(end_date)

    if total:
        return write_total(ens_files, abs_files, ens_means, analysis_str, variables, start_date, end_date, argv,
                               saved, full_saved, lon_centre=lon_centre, mask=mask, lon=lon, lat=lat,
                                grid=grid, test=test)

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
        else:
            print("ERROR in function write_means_to_netcdf_file: Non-NetCDF file discovered " + str(ens_files[i][0]))
            sys.exit()
        output_file = os.path.join(mean_folder, output_file)

        if lon_centre is not None:  # Add lc_'centre' to the name
            output_file = output_file[:-3] + '_lc_' + str(lon_centre) + '.nc'
        if mask is not None:  # Add masked to the name
            output_file = output_file[:-3] + '_masked' + '.nc'

        sample = False
        if lat and lon and not grid:
            sample = True

        if sample:
            output_file = output_file[:-3] + '_s_' + str(lat) + '_' + str(lon) + '.nc'
        elif grid:
            output_file = output_file[:-3] + '_g_' + str(lat) + '_' + str(lon) + '.nc'

        # Combine each ensemble file into one
        times_append = xr.open_mfdataset(abs_files[i])
        times_append.to_netcdf(output_file)

        # Get global atrributes of file to use later
        glob_attrs = times_append.attrs

        # Get previous names
        list_dim_names = list(times_append.dims.keys())
        time_name, lat_name, lon_name = None, None, None
        for dd in list_dim_names:
            if dd.lower() == 't' or 'time' in dd.lower():
                time_name = dd
            elif dd[0].lower() == 'y' or 'lat' in dd.lower() or 'latitude' in dd.lower():
                lat_name = dd
            elif dd[0].lower() == 'x' or 'lon' in dd.lower() or 'longitude' in dd.lower():
                lon_name = dd

        # If sample/ grid point , get original grid
        if full_saved is not None:
            saved = full_saved

        # Turn all cubes to xarray and write them in file
        for var in variables:
            # Save original grid to file
            cube_saved = saved[i][var]
            # Replace to original names
            coord = cube_saved.coord("longitude")
            coord.rename(lon_name)
            coord = cube_saved.coord("latitude")
            coord.rename(lat_name)
            coord = cube_saved.coord("time")
            coord.rename(time_name)
            # Convert cube to xarray
            converted_saved = xr.DataArray.from_iris(cube_saved)
            try:
                if mask is not None:  # Make a new file with masked data
                    converted_saved.to_netcdf(output_file, mode='w')
                else:
                    converted_saved.to_netcdf(output_file, mode='a')
            except ValueError as err:
                print("ERROR in function write_means_to_netcdf_file: " + str(err))

            # Save all analysis to file
            for a in range(len(analysis_str)):
                # Save analysis grid to file
                cube = ens_means[i][a][var]
                # Rename to original names
                coord = cube.coord("longitude")
                coord.rename(lon_name)
                coord = cube.coord("latitude")
                coord.rename(lat_name)
                coord = cube.coord("time")
                coord.rename(time_name)
                # Convert to xarray
                converted = xr.DataArray.from_iris(cube)
                # Add analysis name to variable
                new_name = converted.name + '_' + analysis_str[a]
                converted = converted.rename(new_name)
                # Add new long name of file
                try:
                    new_long_name = converted.attrs['long_name'] + ' averaged between ' + start_end_str
                except KeyError:
                    new_long_name = var + ' averaged between ' + start_end_str
                converted.attrs['long_name'] = new_long_name
                # Append analysis to file
                converted.to_netcdf(output_file, mode='a')

        # Global attributes
        dest = Dataset(output_file, 'a')
        # Write to description and history of file
        desc_str = "Added " + ', '.join(analysis_str) + " of variables " + ', '.join(variables) + " within time period " + \
                   start_end_str
        if grid:
            desc_str = "Added " + ', '.join(analysis_str) + " of variables " + ', '.join(
                variables) + " using grid point (" + str(lat) + "," + str(lon) + ")" + " within time period " + \
                       start_end_str
        elif sample:
            desc_str = "Added " + ', '.join(analysis_str) + " of variables " + ', '.join(
                variables) + " using sample point (" + str(lat) + "," + str(lon) + ")" + " within time period " + \
                       start_end_str
        hist_str = time.ctime(time.time()) + ': Commands used to produce file: ' + argv + ' \n ' + \
                           time.ctime(time.time()) + ': Functions used:  extract_data, compute_stats_analysis,' \
                                                     ' write_means_to_netcdf_file'
        if 'description' in dest.ncattrs():
            dest.description = desc_str + ' \n ' + dest.description
        else:
            dest.description = desc_str
        if 'history' in dest.ncattrs():
            dest.history = hist_str + ' \n ' + dest.history
        else:
            dest.history = hist_str

        # Add global attributes from original files
        if mask is not None:
            for dest_key, dest_val in glob_attrs.items():
                if dest_key == 'description':
                    dest.description = dest.description + ' \n ' + dest_val
                elif dest_key == 'history':
                    dest.history = dest.history + ' \n ' + dest_val
                elif dest_key == 'filename':
                    dest.filename = dest_val
                elif dest_key == 'title':
                    dest.title = dest_val
                elif dest_key == 'grid_type':
                    dest.grid_type = dest_val
                elif dest_key == 'grid_tile':
                    dest.grid_tile = dest_val
                elif dest_key == 'source':
                    dest.source = dest_val
        dest.close()

    print("Ensemble files created in " + os.path.join(directories.ANALYSIS, directories.MEANS) + " folder.")


def write_total(ens_files, abs_files, ens_means, analysis_str, variables, start_date, end_date, argv,
                               saved, full_saved, lon_centre=None, mask=None, lon=None, lat=None,
                                grid=None, test=False):
    """
    Write analysis computed in netcdf files if ensembles calculated together
    :param ens_files: initial files arranged in ensemble order
    :param abs_files: absolute path of ens_files
    :param ens_means: ensemble means calculated calling function compute_stats_analysis
    :param analysis_str: type of analysis computed in ens_means: 'mean', 'std', 'median', 'rms' or 'all' of them
    :param variables: list of variables
    :param start_date: start date list in [day, month, year] format
    :param end_date: end date list in [day, month, year] format
    :param argv: string containing command line arguments used
    :param saved: iris cubes containing data
    :param full_saved: iris cubes containing data if sample/grid point are used
    :param total: set if all ensembles have been calculated together, instead of separately, boolean
    :param lon_centre: longitude center , float
    :param mask: file that contains mask data, string
    :param lon: longitude, set if grid or sample point, floats
    :param lat: latitude, set if grid or sample point, floats
    :param grid: set if grid point is given
    :param test: if test is true, make some changes specific to files on my pc
    :return: None, files created in folder analysis/ensemble_means
    """
    # Assertions
    assert ens_files is not None and abs_files is not None and variables is not None and analysis_str is not None
    assert check_list_date(start_date) and check_list_date(end_date)

    # Start and end date string
    start_end_str = str(start_date[2]) + "-" + str(start_date[1]) + "-" + str(start_date[0]) + " and " + \
                    str(end_date[2]) + "-" +str(end_date[1]) + "-" + str(end_date[0])

    # Get folder to store ensemble means
    results = directories.ANALYSIS
    mean_folder = os.path.abspath(os.path.join(results, directories.MEANS))
    if test:
        mean_folder = mean_folder.replace("Adanna Akwataghibe", "Adanna")

    ens_members = []

    for i in range(len(ens_means)):
        if ens_files[i][0].endswith(".nc"):
            # Combine each ensemble file into one
            times_append = xr.open_mfdataset(abs_files[i])
            ens_members.append(times_append)

    appended = xr.concat(ens_members, 'ensemble_member')

    # Get first file name in specific ensemble and add last year to name - use as output file name
    if ens_files[0][0].endswith(".nc"):
        # append start and end date to file name
        file_prefix = ens_files[0][0].split('_ens')[0]
        output_file = file_prefix + '_all_ens' + '_' + str(start_date[2]) + '_' + str(end_date[2]) + '.nc'
        output_file = os.path.basename(os.path.normpath(output_file))
    else:
        print("ERROR in function write_means_to_netcdf_file: Non-NetCDF file discovered " + str(ens_files[i][0]))
        sys.exit()
    output_file = os.path.join(mean_folder, output_file)
    if lon_centre is not None:  # Add lc_'centre' to the name
        output_file = output_file[:-3] + '_lc_' + str(lon_centre) + '.nc'
    if mask is not None:  # Add masked to the name
        output_file = output_file[:-3] + '_masked' + '.nc'

    sample = False
    if lat and lon and not grid:
        sample = True

    if sample:
        output_file = output_file[:-3] + '_s_' + str(lat) + '_' + str(lon) + '.nc'
    elif grid:
        output_file = output_file[:-3] + '_g_' + str(lat) + '_' + str(lon) + '.nc'

    appended.to_netcdf(output_file)

    # Get global atrributes of file to use later
    glob_attrs = appended.attrs

    # Get previous names
    list_dim_names = list(appended.dims.keys())
    time_name, lat_name, lon_name = None, None, None
    for dd in list_dim_names:
        if dd.lower() == 't' or 'time' in dd.lower():
            time_name = dd
        elif dd[0].lower() == 'y' or 'lat' in dd.lower() or 'latitude' in dd.lower():
            lat_name = dd
        elif dd[0].lower() == 'x' or 'lon' in dd.lower() or 'longitude' in dd.lower():
            lon_name = dd

    # If sample/ grid point , get original grid
    if full_saved is not None:
        saved = full_saved

    changed_names = {"longitude": lon_name, "latitude": lat_name, "time":time_name}

    # Turn all cubes to xarray and write them in file
    for var in variables:
        # Save original grid to file (or masked/lon centre)
        ens_arrays = []
        for i in range(len(ens_means)):
            cube = saved[i][var]
            ens_arrays.append(xr.DataArray.from_iris(cube))
        xr_saved = xr.concat(ens_arrays, 'ensemble_member')
        # Replace to original names
        xr_saved = xr_saved.rename(changed_names)
        print(xr_saved)
        try:
            if mask is not None:  # Make a new file with masked data
                xr_saved.to_netcdf(output_file, mode='w')
            else:
                xr_saved.to_netcdf(output_file, mode='a')
        except ValueError as err:
            print("ERROR in function write_means_to_netcdf_file: " + str(err))

        # Save all analysis to file
        for a in range(len(analysis_str)):
            # Save analysis grid to file
            cube = ens_means[a][var]
            # Rename to original names
            coord = cube.coord("longitude")
            coord.rename(lon_name)
            coord = cube.coord("latitude")
            coord.rename(lat_name)
            coord = cube.coord("time")
            coord.rename(time_name)
            # Convert to xarray
            converted = xr.DataArray.from_iris(cube)
            # Add analysis name to variable
            new_name = converted.name + '_' + analysis_str[a]
            converted = converted.rename(new_name)
            # Add new long name of file
            new_long_name = converted.attrs['long_name'] + ' averaged between ' + start_end_str
            converted.attrs['long_name'] = new_long_name
            # Append analysis to file
            converted.to_netcdf(output_file, mode='a')

    # Global attributes
    dest = Dataset(output_file, 'a')
    # Write to description and history of file
    desc_str = "Added " + ', '.join(analysis_str) + " of variables " + ', '.join(variables) + " within time period " + \
               start_end_str
    if grid:
        desc_str = "Added " + ', '.join(analysis_str) + " of variables " + ', '.join(
            variables) + " using grid point (" + str(lat) + "," + str(lon) + ")" + " within time period " + \
                   start_end_str
    elif sample:
        desc_str = "Added " + ', '.join(analysis_str) + " of variables " + ', '.join(
            variables) + " using sample point (" + str(lat) + "," + str(lon) + ")" + " within time period " + \
                   start_end_str
    hist_str = time.ctime(time.time()) + ': Commands used to produce file: ' + argv + ' \n ' + \
                       time.ctime(time.time()) + ': Functions used:  extract_data, compute_stats_analysis,' \
                                                 ' write_means_to_netcdf_file'
    if 'description' in dest.ncattrs():
        dest.description = desc_str + ' \n ' + dest.description
    else:
        dest.description = desc_str
    if 'history' in dest.ncattrs():
        dest.history = hist_str + ' \n ' + dest.history
    else:
        dest.history = hist_str

    # Add global attributes from original files
    if mask is not None:
        for dest_key, dest_val in glob_attrs.items():
            if dest_key == 'description':
                dest.description = dest.description + ' \n ' + dest_val
            elif dest_key == 'history':
                dest.history = dest.history + ' \n ' + dest_val
            elif dest_key == 'filename':
                dest.filename = dest_val
            elif dest_key == 'title':
                dest.title = dest_val
            elif dest_key == 'grid_type':
                dest.grid_type = dest_val
            elif dest_key == 'grid_tile':
                dest.grid_tile = dest_val
            elif dest_key == 'source':
                dest.source = dest_val
    dest.close()


def write_user_analysis_to_netcdf_file(ens_files, abs_files, user_analysis, func_name, variables, start_date, end_date,
                                       argv, test=False):
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



