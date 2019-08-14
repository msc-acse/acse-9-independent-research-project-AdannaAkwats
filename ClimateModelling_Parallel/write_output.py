import time
from utils import *
from netCDF4 import Dataset, MFDataset
import pickle
import xarray as xr
from multiprocessing import Pool
from functools import partial
import parallel_settings
import sys
import write_output_serial

def write_parallel(abs_files, ens_means, analysis_str, variables, start_date, end_date, argv,
                               saved, full_saved, lon_centre, mask, lon, lat,
                                grid, mean_folder, start_end_str, user_func, tup_ens):
    """
    Write to files using processes
    :param abs_files: absolute path of ens_files
    :param ens_means: ensemble means calculated calling function compute_stats_analysis
    :param analysis_str: type of analysis computed in ens_means: 'mean', 'std', 'median', 'rms' or 'all' of them
    :param variables: list of variables
    :param start_date: start date list in [day, month, year] format
    :param end_date: end date list in [day, month, year] format
    :param argv: string containing command line arguments used
    :param saved: iris cubes containing data
    :param full_saved: iris cubes containing data if sample/grid point are used
    :param lon_centre: longitude center , float
    :param mask: file that contains mask data, string
    :param lon: longitude, set if grid or sample point, floats
    :param lat: latitude, set if grid or sample point, floats
    :param grid: set if grid point is given
    :param test: if test is true, make some changes specific to files on my pc
    :param mean_folder: path of where to store averages in reults
    :param start_end_str: string with start and end date
    :param user_func: user function name
    :param tup_ens: tuple containing ensemble files and corresponding index in original list
    :return: None, files created in folder analysis/ensemble_means
    """

    i, ens = tup_ens

    # Get first file name in specific ensemble and add last year to name - use as output file name
    output_file = ""
    if ens[0].endswith(".nc"):
        # append start and end date to file name
        # Check if file name has 2 times in it
        if not get_file_two_years(ens[0]):
            output_file = re.sub(str(start_date[2]), str(start_date[2]) + '_' + str(end_date[2]), ens[0])
        else:
            output_file = re.sub(r'_(\d+)_(\d+)', '_' + str(start_date[2]) + '_' + str(end_date[2]), ens[0])
        output_file = os.path.basename(os.path.normpath(output_file))
    else:
        print("ERROR in function write_analysis_to_netcdf_file: Non-NetCDF file discovered " + str(ens[0]))
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
    if len(full_saved[0]) != 0:
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
            converted_saved.to_netcdf(output_file, mode='a')
        except Exception:
            converted_saved.to_netcdf(output_file, mode='w')

        # If user function given, change analysis string
        if user_func:
            cube = ens_means[i][var]
            # Convert to xarray
            converted = xr.DataArray.from_iris(cube)
            # Add analysis name to variable
            new_name = converted.name + '_' + user_func
            converted = converted.rename(new_name)
            # Add new long name of file
            try:
                new_long_name = converted.attrs['long_name'] + ' averaged between ' + start_end_str
            except KeyError:
                new_long_name = var + ' averaged between ' + start_end_str
            converted.attrs['long_name'] = new_long_name
            # Append analysis to file
            converted.to_netcdf(output_file, mode='a')

        else:
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
    as_ = ', '.join(analysis_str)
    if user_func:
        as_ = user_func
    # Write to description and history of file
    desc_str = "Added " + as_ + " of variables " + ', '.join(variables) + " within time period " + \
               start_end_str
    if grid:
        desc_str = "Added " + as_ + " of variables " + ', '.join(
            variables) + " using grid point (" + str(lat) + "," + str(lon) + ")" + " within time period " + \
                   start_end_str
    elif sample:
        desc_str = "Added " + as_ + " of variables " + ', '.join(
            variables) + " using sample point (" + str(lat) + "," + str(lon) + ")" + " within time period " + \
                   start_end_str
    hist_str = time.ctime(time.time()) + ': Commands used to produce file: ' + argv + ' \n ' + \
               time.ctime(time.time()) + ': Functions used:  extract_data, compute_stats_analysis,' \
                                         ' write_analysis_to_netcdf_file'
    if user_func:
        hist_str = time.ctime(time.time()) + ': Commands used to produce file: ' + argv + ' \n ' + \
                   time.ctime(time.time()) + ': Functions used:  extract_data, compute_user_analysis,' \
                                             ' write_analysis_to_netcdf_file'

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


def write_analysis_to_netcdf_file(ens_files, abs_files, ens_means, analysis_str, variables, start_date, end_date, argv,
                               saved, full_saved, total=False, lon_centre=None, mask=None, lon=None, lat=None,
                                grid=None, user_func=False, test=False):
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
    :param user_func: user function name
    :param test: if test is true, make some changes specific to files on my pc
    :return: None, files created in folder analysis/ensemble_means
    """
    # Assertions
    assert ens_files is not None
    assert abs_files is not None
    assert variables is not None
    assert analysis_str is not None
    assert check_list_date(start_date) and check_list_date(end_date)

    if total:
        return write_total(ens_files, abs_files, ens_means, analysis_str, variables, start_date, end_date, argv,
                               saved, full_saved, lon_centre=lon_centre, mask=mask, lon=lon, lat=lat,
                                grid=grid, test=test)

    else:
        # Start and end date string
        start_end_str = str(start_date[2]) + "-" + str(start_date[1]) + "-" + str(start_date[0]) + " and " + \
                        str(end_date[2]) + "-" +str(end_date[1]) + "-" + str(end_date[0])

        # Get folder to store ensemble means
        results = directories.ANALYSIS
        mean_folder = os.path.abspath(os.path.join(results, directories.MEANS))
        if test:
            mean_folder = mean_folder.replace("Adanna Akwataghibe", "Adanna")

        try:
            pool = Pool(processes=parallel_settings.NUM_PROCESSORS)
            func = partial(write_parallel, abs_files, ens_means, analysis_str, variables, start_date, end_date, argv,
                           saved, full_saved, lon_centre, mask, lon, lat,
                           grid, mean_folder, start_end_str, user_func)
            tup_ens = (list(range(len(ens_means))), ens_files)
            pool.map(func, zip(list(range(len(ens_means))), ens_files))
            pool.close()
            pool.join()
        except Exception:  # If size is too big to send then switch to serial
            write_output_serial.write_analysis_to_netcdf_file(ens_files, abs_files, ens_means, analysis_str, variables,
                                                              start_date, end_date, argv,saved, full_saved, total=total,
                                                              lon_centre=lon_centre, mask=mask, lon=lon, lat=lat,
                                                              grid=grid, user_func=user_func, test=test)

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
    assert ens_files is not None
    assert abs_files is not None
    assert variables is not None
    assert analysis_str is not None
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
        print("ERROR in function write_analysis_to_netcdf_file: Non-NetCDF file discovered " + str(ens_files[i][0]))
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
    if len(full_saved[0]) != 0:
        saved = full_saved

    changed_names = {"longitude": lon_name, "latitude": lat_name, "time": time_name}

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
        try:
            xr_saved.to_netcdf(output_file, mode='a')
        except Exception:
            xr_saved.to_netcdf(output_file, mode='w')

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
    as_ = ', '.join(analysis_str)
    desc_str = "Added " + as_ + " of variables " + ', '.join(variables) + " within time period " + \
               start_end_str
    if grid:
        desc_str = "Added " + as_ + " of variables " + ', '.join(
            variables) + " using grid point (" + str(lat) + "," + str(lon) + ")" + " within time period " + \
                   start_end_str
    elif sample:
        desc_str = "Added " + as_ + " of variables " + ', '.join(
            variables) + " using sample point (" + str(lat) + "," + str(lon) + ")" + " within time period " + \
                   start_end_str
    hist_str = time.ctime(time.time()) + ': Commands used to produce file: ' + argv + ' \n ' + \
                       time.ctime(time.time()) + ': Functions used:  extract_data, compute_stats_analysis,' \
                                                 ' write_analysis_to_netcdf_file'
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



