"""
Created by Adanna Akwataghibe (Github: AdannaAkwats)
"""
import time
from utils import *
from netCDF4 import Dataset
import xarray as xr


class WriteOutput:
    def __init__(self, ens_files, abs_files, ens_means, analysis_str, variables, start_date, end_date, argv,
                               saved, full_saved, total=False, lon_centre=None, mask=None, lon=None, lat=None,
                                grid=None, user_func=False, points_sample_grid=None, second_date_given=False,test=False):
        self.ens_files = ens_files
        self.abs_files = abs_files
        self.ens_means = ens_means
        self.analysis_str = analysis_str
        self.variables = variables
        self.start_date = start_date
        self.end_date = end_date
        self.argv = argv
        self.saved = saved
        self.full_saved = full_saved
        self.total = total
        self.lon_centre = lon_centre
        self.mask = mask
        self.lon = lon
        self.lat = lat
        self.grid = grid
        self.user_func = user_func
        self.points_sample_grid = points_sample_grid
        self.second_date_given = second_date_given
        self.test = test


    def write_analysis_to_netcdf_file(self):
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
        :param points_sample_grid: file (txt or nc) that contains sample or grid points
        :param second_date_given: if set, multi model averages calculated
        :param start_date2: second model start date
        :param end_date2: second model end date
        :param test: if test is true, make some changes specific to files on my pc
        :return: None, files created in folder analysis/ensemble_means
        """
        # Assertions
        assert self.ens_files is not None
        assert self.abs_files is not None
        assert self.variables is not None
        assert check_list_date(self.start_date) and check_list_date(self.end_date)

        if self.total:
            return self.write_total()

        # Start and end date string
        start_end_str = str(self.start_date[2]) + "-" + str(self.start_date[1]) + "-" + str(self.start_date[0]) + " and " + \
                        str(self.end_date[2]) + "-" +str(self.end_date[1]) + "-" + str(self.end_date[0])

        # Get folder to store ensemble means
        results = directories.ANALYSIS
        mean_folder = os.path.abspath(os.path.join(results, directories.MEANS))
        if self.test:
            mean_folder = mean_folder.replace("Adanna Akwataghibe", "Adanna")

        if not self.analysis_str and self.user_func is None:
            self.ens_means = self.ens_files

        # Go through ensembles, merge files to get output and write to output
        for i in range(len(self.ens_means)):
            # Get first file name in specific ensemble and add last year to name - use as output file name
            output_file = ""
            if self.ens_files[i][0].endswith(".nc"):
                # Check if file name has 2 times in it
                if not get_file_two_years(self.ens_files[i][0]):
                    output_file = re.sub(str(self.start_date[2]), str(self.start_date[2]) + '_' + str(self.end_date[2]), self.ens_files[i][0])
                    if self.second_date_given:
                        output_file = re.sub(str(self.start_date[2]), 'multi_model', self.ens_files[i][0])
                else:
                    output_file = re.sub(r'_(\d+)_(\d+)', '_' + str(self.start_date[2]) + '_' + str(self.end_date[2]), self.ens_files[i][0])
                    output_file = re.sub(r'_(\d+)-(\d+)', '_' + str(self.start_date[2]) + '-' + str(self.end_date[2]), self.ens_files[i][0])
                    if self.second_date_given:
                        output_file = re.sub(r'_(\d+)_(\d+)', '_multi_model', self.ens_files[i][0])
                        output_file = re.sub(r'_(\d+)-(\d+)', '_multi_model', self.ens_files[i][0])
                # append start and end date to file name
                output_file = os.path.basename(os.path.normpath(output_file))
            else:
                print("ERROR in function write_analysis_to_netcdf_file: Non-NetCDF file discovered " + str(self.ens_files[i][0]))
                sys.exit()
            output_file = os.path.join(mean_folder, output_file)

            if self.lon_centre is not None:  # Add lc_'centre' to the name
                output_file = output_file[:-3] + '_lc_' + str(self.lon_centre) + '.nc'
            if self.mask is not None:  # Add masked to the name
                output_file = output_file[:-3] + '_masked' + '.nc'
            if self.user_func is not None:  # Add user function to the name
                output_file = output_file[:-3] + '_' + self.user_func + '.nc'

            sample, nc_true, txt_true = False, False, False
            if self.lat and self.lon and not self.grid:
                sample = True
            if self.points_sample_grid is not None and not self.grid:
                sample = True

            if sample:
                if self.points_sample_grid is not None:
                    output_file = output_file[:-3] + '_s_' + 'regrid' + '.nc'
                else:
                    output_file = output_file[:-3] + '_s_' + str(self.lat) + '_' + str(self.lon) + '.nc'
            elif self.grid:
                if self.points_sample_grid is not None:
                    output_file = output_file[:-3] + '_g_' + 'regrid' + '.nc'
                else:
                    output_file = output_file[:-3] + '_g_' + str(self.lat) + '_' + str(self.lon) + '.nc'

            # Combine each ensemble file into one
            times_append = xr.open_mfdataset(self.abs_files[i])
            try:
                times_append.to_netcdf(path=output_file, mode='w')
            except Exception as err:
                print("ERROR in write_analysis_to_netcdf_file: " + str(err))

            # Get global atrributes of file to use later
            glob_attrs = times_append.attrs

            # Get previous names
            list_dim_names = list(times_append.dims.keys())
            time_name, lat_name, lon_name = None, None, None
            for dd in list_dim_names:
                not_bnds = 'bound' not in dd and 'bnd' not in dd
                if not_bnds:
                    if dd.lower() == 't' or dd.lower() == 'time':
                        time_name = dd
                    elif dd[0].lower() == 'y' or 'lat' in dd.lower() or 'latitude' in dd.lower():
                        lat_name = dd
                    elif dd[0].lower() == 'x' or 'lon' in dd.lower() or 'longitude' in dd.lower():
                        lon_name = dd

            # If sample/ grid point , get original grid
            if self.full_saved[0] is not None:
                self.saved = self.full_saved

            # Turn all cubes to xarray and write them in file
            for var in self.variables:
                # Save original grid to file
                cube_saved = self.saved[i][var]
                # Replace to original names
                try:
                    coord = cube_saved.coord("longitude")
                    coord.rename(lon_name)
                    coord = cube_saved.coord("latitude")
                    coord.rename(lat_name)
                except Exception:
                    pass
                try:
                    coord = cube_saved.coord("time")
                    coord.rename(time_name)
                except Exception:
                    pass
                # Convert cube to xarray
                converted_saved = xr.DataArray.from_iris(cube_saved)
                try:
                    converted_saved.to_netcdf(path=output_file, mode='a')
                except Exception:
                    converted_saved.to_netcdf(path=output_file, mode='w')

                if self.analysis_str:
                    # Save all analysis to file
                    for a in range(len(self.analysis_str)):
                        # Save analysis grid to file
                        cube = self.ens_means[i][a][var]
                        d = np.asarray(cube.data)
                        cube.data = d
                        # Rename to original names
                        try:
                            coord = cube.coord("longitude")
                            coord.rename(lon_name)
                            coord = cube.coord("latitude")
                            coord.rename(lat_name)
                        except Exception:
                            pass
                        try:  # Time may not always be here if we take the averages
                            coord = cube.coord("time")
                            coord.rename(time_name)
                        except Exception:
                            pass
                        # Convert to xarray
                        converted = xr.DataArray.from_iris(cube)
                        # Add analysis name to variable
                        new_name = converted.name + '_' + self.analysis_str[a]
                        converted = converted.rename(new_name)
                        # Add new long name of file
                        try:
                            new_long_name = converted.attrs['long_name'] + ' averaged between ' + start_end_str
                            if self.second_date_given:
                                new_long_name = converted.attrs['long_name'] + ' multi model average'
                        except KeyError:
                            new_long_name = var + ' averaged between ' + start_end_str
                            if self.second_date_given:
                                new_long_name = converted.attrs['long_name'] + ' multi model average'
                        converted.attrs['long_name'] = new_long_name

                        if self.second_date_given:  # Write to new file
                            converted.to_netcdf(output_file, mode='w')
                        else:
                            # Append analysis to file
                            try:
                                converted.to_netcdf(path=output_file, mode='a')
                            except Exception:
                                converted.to_netcdf(path=output_file, mode='w')

            # If user function given, change analysis string
            if self.user_func:
                cube = self.ens_means[i]
                name = cube.name()
                # Convert to xarray
                converted = xr.DataArray.from_iris(cube)
                # Add analysis name to variable
                new_name = converted.name + '_' + self.user_func
                converted = converted.rename(new_name)
                # Add new long name of file
                try:
                    new_long_name = converted.attrs['long_name'] + ' averaged between ' + start_end_str
                except KeyError:
                    new_long_name = name + ' averaged between ' + start_end_str
                converted.attrs['long_name'] = new_long_name
                # Append analysis to file
                try:
                    converted.to_netcdf(path=output_file, mode='a')
                except Exception:
                    converted.to_netcdf(path=output_file, mode='w')

            # Global attributes
            dest = Dataset(output_file, 'a')
            if not self.analysis_str:
                self.analysis_str = ['analysis']
            as_ = ', '.join(self.analysis_str)
            if self.user_func:
                as_ = self.user_func
            # Write to description and history of file
            desc_str = "Added " + as_ + " of variables " + ', '.join(self.variables) + " within time period " + \
                       start_end_str
            if self.grid:
                desc_str = "Added " + as_ + " of variables " + ', '.join(
                    self.variables) + " using grid point (" + str(self.lat) + "," + str(self.lon) + ")" + " within time period " + \
                           start_end_str
            elif sample:
                desc_str = "Added " + as_ + " of variables " + ', '.join(
                    self.variables) + " using sample point (" + str(self.lat) + "," + str(self.lon) + ")" + " within time period " + \
                           start_end_str
            if self.points_sample_grid is not None:
                desc_str = "Added " + as_ + " of variables " + ', '.join(
                    self.variables) + " regridded using points in " + self.points_sample_grid + " within time period " + \
                           start_end_str
            hist_str = time.ctime(time.time()) + ': Commands used to produce file: ' + self.argv + ' \n ' + \
                               time.ctime(time.time()) + ': Functions used:  extract_data, compute_stats_analysis,' \
                                                         ' write_analysis_to_netcdf_file'
            if self.user_func:
                hist_str = time.ctime(time.time()) + ': Commands used to produce file: ' + self.argv + ' \n ' + \
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
            if self.mask is not None:
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

            # If no analysis, set back to false for next ensemble
            if self.analysis_str == ['analysis']:
                self.analysis_str = False

        print("Ensemble files created in " + os.path.join(directories.ANALYSIS, directories.MEANS) + " folder.")


    def write_total(self):
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
        :param points_sample_grid: file (txt or nc) that contains sample or grid points
        :param second_date_given: if set, multi model averages calculated
        :param start_date2: second model start date
        :param end_date2: second model end date
        :param test: if test is true, make some changes specific to files on my pc
        :return: None, files created in folder analysis/ensemble_means
        """
        # Assertions
        assert self.ens_files is not None
        assert self.abs_files is not None
        assert self.variables is not None
        assert self.analysis_str is not None
        assert check_list_date(self.start_date) and check_list_date(self.end_date)

        # Start and end date string
        start_end_str = str(self.start_date[2]) + "-" + str(self.start_date[1]) + "-" + str(self.start_date[0]) + " and " + \
                        str(self.end_date[2]) + "-" +str(self.end_date[1]) + "-" + str(self.end_date[0])

        # Get folder to store ensemble means
        results = directories.ANALYSIS
        mean_folder = os.path.abspath(os.path.join(results, directories.MEANS))
        if self.test:
            mean_folder = mean_folder.replace("Adanna Akwataghibe", "Adanna")

        ens_members = []

        for i in range(len(self.ens_files)):
            if self.ens_files[i][0].endswith(".nc"):
                # Combine each ensemble file into one
                times_append = xr.open_mfdataset(self.abs_files[i])
                ens_members.append(times_append)

        appended = xr.concat(ens_members, 'ensemble_member')

        # Get first file name in specific ensemble and add last year to name - use as output file name
        if self.ens_files[0][0].endswith(".nc"):
            # append start and end date to file name
            file_prefix = self.ens_files[0][0].split('_ens')[0]
            output_file = file_prefix + '_all_ens' + '_' + str(self.start_date[2]) + '_' + str(self.end_date[2]) + '.nc'
            output_file = os.path.basename(os.path.normpath(output_file))
        else:
            print("ERROR in function write_analysis_to_netcdf_file: Non-NetCDF file discovered " + str(self.ens_files[i][0]))
            sys.exit()
        output_file = os.path.join(mean_folder, output_file)
        if self.lon_centre is not None:  # Add lc_'centre' to the name
            output_file = output_file[:-3] + '_lc_' + str(self.lon_centre) + '.nc'
        if self.mask is not None:  # Add masked to the name
            output_file = output_file[:-3] + '_masked' + '.nc'

        sample = False
        if self.lat and self.lon and not self.grid:
            sample = True
        if self.points_sample_grid is not None and not self.grid:
            sample = True

        if sample:
            if self.points_sample_grid is not None:
                output_file = output_file[:-3] + '_s_' + 'regrid' + '.nc'
            else:
                output_file = output_file[:-3] + '_s_' + str(self.lat) + '_' + str(self.lon) + '.nc'
        elif self.grid:
            if self.points_sample_grid is not None:
                output_file = output_file[:-3] + '_g_' + 'regrid' + '.nc'
            else:
                output_file = output_file[:-3] + '_g_' + str(self.lat) + '_' + str(self.lon) + '.nc'

        appended.to_netcdf(path=output_file, mode='w')

        # Get global atrributes of file to use later
        glob_attrs = appended.attrs

        # Get previous names
        list_dim_names = list(appended.dims.keys())
        time_name, lat_name, lon_name = None, None, None
        for dd in list_dim_names:
            not_bnds = 'bound' not in dd and 'bnd' not in dd
            if not_bnds:
                if dd.lower() == 't' or dd.lower() == 'time':
                    time_name = dd
                elif dd[0].lower() == 'y' or 'lat' in dd.lower() or 'latitude' in dd.lower():
                    lat_name = dd
                elif dd[0].lower() == 'x' or 'lon' in dd.lower() or 'longitude' in dd.lower():
                    lon_name = dd

        # If sample/ grid point , get original grid
        if self.full_saved[0] is not None:
            self.saved = self.full_saved

        changed_names = {"longitude": lon_name, "latitude": lat_name, "time":time_name}

        # Turn all cubes to xarray and write them in file
        for var in self.variables:
            # Save original grid to file (or masked/lon centre)
            ens_arrays = []
            for i in range(len(self.ens_files)):
                cube = self.saved[i][var]
                ens_arrays.append(xr.DataArray.from_iris(cube))
            xr_saved = xr.concat(ens_arrays, 'ensemble_member')
            # Replace to original names
            xr_saved = xr_saved.rename(changed_names)
            try:
                xr_saved.to_netcdf(path=output_file, mode='a')
            except Exception:
                xr_saved.to_netcdf(path=output_file, mode='w')

            # If no analysis
            if not self.analysis_str:
                for var in self.variables:
                    # Save original grid to file
                    cube_saved = self.ens_means[0][var]
                    # Replace to original names
                    try:
                        coord = cube_saved.coord("longitude")
                        coord.rename(lon_name)
                        coord = cube_saved.coord("latitude")
                        coord.rename(lat_name)
                    except Exception:
                        pass
                    try:
                        coord = cube_saved.coord("time")
                        coord.rename(time_name)
                    except Exception:
                        pass
                    # Convert cube to xarray
                    converted_saved = xr.DataArray.from_iris(cube_saved)
                    try:
                        converted_saved.to_netcdf(path=output_file, mode='a')
                    except Exception:
                        converted_saved.to_netcdf(path=output_file, mode='w')

            else:
                # Save all analysis to file
                for a in range(len(self.analysis_str)):
                    # Save analysis grid to file
                    cube = self.ens_means[a][var]
                    d = np.asarray(cube.data)
                    cube.data = d
                    # Rename to original names
                    try:
                        coord = cube.coord("longitude")
                        coord.rename(lon_name)
                        coord = cube.coord("latitude")
                        coord.rename(lat_name)
                    except Exception:
                        pass
                    try:
                        coord = cube.coord("time")
                        coord.rename(time_name)
                    except Exception:
                        pass
                    # Convert to xarray
                    converted = xr.DataArray.from_iris(cube)
                    # Add analysis name to variable
                    new_name = converted.name + '_' + self.analysis_str[a]
                    converted = converted.rename(new_name)
                    # Add new long name of file
                    try:
                        new_long_name = converted.attrs['long_name'] + ' averaged between ' + start_end_str
                        if self.second_date_given:
                            new_long_name = converted.attrs['long_name'] + ' multi model average'
                    except KeyError:
                        new_long_name = var + ' averaged between ' + start_end_str
                        if self.second_date_given:
                            new_long_name = converted.attrs['long_name'] + ' multi model average'
                    converted.attrs['long_name'] = new_long_name
                    if self.second_date_given:  # Write to new file
                        converted.to_netcdf(path=output_file, mode='w')
                    else:
                        # Append analysis to file
                        try:
                            converted.to_netcdf(path=output_file, mode='a')
                        except Exception:
                            converted.to_netcdf(path=output_file, mode='w')

        # Global attributes
        dest = Dataset(output_file, 'a')
        # Write to description and history of file
        if not self.analysis_str:
            self.analysis_str = ['analysis']
        as_ = ', '.join(self.analysis_str)
        desc_str = "Added " + as_ + " of variables " + ', '.join(self.variables) + " within time period " + \
                   start_end_str
        if self.grid:
            desc_str = "Added " + as_ + " of variables " + ', '.join(
                self.variables) + " using grid point (" + str(self.lat) + "," + str(self.lon) + ")" + " within time period " + \
                       start_end_str
        elif sample:
            desc_str = "Added " + as_ + " of variables " + ', '.join(
                self.variables) + " using sample point (" + str(self.lat) + "," + str(self.lon) + ")" + " within time period " + \
                       start_end_str
        if self.points_sample_grid is not None:
            desc_str = "Added " + as_ + " of variables " + ', '.join(
                self.variables) + " regridded using points in " + self.points_sample_grid + " within time period " + \
                       start_end_str
        hist_str = time.ctime(time.time()) + ': Commands used to produce file: ' + self.argv + ' \n ' + \
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
        if self.mask is not None:
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


