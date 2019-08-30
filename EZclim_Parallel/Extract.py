"""
Created by Adanna Akwataghibe (Github: AdannaAkwats)
"""
from utils import *
import directories
import xarray as xr
from netCDF4 import Dataset
import time

def warn(*args, **kwargs):
    pass
import warnings

warnings.warn = warn
import iris
import warnings
from functools import reduce
import operator
import math
from shapely.geometry import Polygon, Point
from multiprocessing import Pool
from functools import partial
import parallel_settings
warnings.filterwarnings("ignore")



class Extract:
    def __init__(self, algae_type, variables, start_date, end_date, num_ens, monthly=False, lat=None, lon=None, grid=None,
                 points_sample_grid=None, lon_centre=None, maskfile=None, calc_areas=False, test=True):
        self.algae_type = algae_type
        self.variables = variables
        self.start_date = start_date
        self.end_date = end_date
        self.num_ens = num_ens
        self.monthly = monthly
        self.lat = lat
        self.lon = lon
        self.grid = grid
        self.points_sample_grid = points_sample_grid
        self.lon_centre = lon_centre
        self.maskfile = maskfile
        self.calc_areas = calc_areas
        self.test = test


    def contains_points(self, polygon, xp, yp):
        """
        Construct boolean array that shows if coords xp and yp are in polygon
        :param polygon: shapely.geometry.Polygon
        :param xp: list of x coordinates
        :param yp: list of y coordinates
        :return: boolean array
        """

        return np.array([Point(x, y).intersects(polygon) for x, y in zip(xp, yp)],
                        dtype=np.bool)


    def get_mask(self, maskfile, cube_data, lons, lats):
        """
        Construct mask (boolean) array from mask file
        :param maskfile: file that contains data for mask
        :param cube_data: 3D/4D array
        :param lons: list of longitudes
        :param lats: list of latitudes
        :return: extracted data using mask
        """

        # Assertions
        assert maskfile is not None
        assert lons is not None
        assert lats is not None

        # Get polygons from mask file
        maskfile = os.path.join(directories.INPUT, maskfile)
        polygons, level = get_polygons(maskfile)

        # if only one polygon
        if not is_nested_list(polygons):
            # Cast nested list as nested list of tuples
            polygons = [tuple(l) for l in polygons]
            # Change to 3D polygon
            polygons = [polygons]

        # Initialise mask array to list of Falses
        full_mask = [False for _ in range(len(lats) * len(lons))]

        # Go through each polygon in the list
        for i in range(len(polygons)):
            # Sort coordinates of polygons to be in clockwise direction
            center = tuple(map(operator.truediv, reduce(lambda x, y: map(operator.add, x, y), polygons[i]),
                               [len(polygons[i])] * 2))
            sorted_poly = sorted(polygons[i], key=lambda coord: (-135 - math.degrees(
                math.atan2(*tuple(map(operator.sub, coord, center))[::-1]))) % 360)

            # Construct shapely polygon and check for points inside polygon
            poly = Polygon(sorted_poly)
            x, y = np.meshgrid(lons, lats)

            # Get mask
            mask = self.contains_points(poly, x.ravel(), y.ravel())
            full_mask = np.logical_or(full_mask, mask)

        # Reshape mask from 1D to 2D array
        m = full_mask.reshape(len(lats), len(lons))

        # Increase 2D array to 3D array
        tiles = None
        cube_shape = cube_data.shape
        if level is None:
            level = [1, cube_shape[0]]

        # Check level is within limits
        if len(level) == 1:
            if level[0] > cube_data.shape[0] or level[0] < 0:
                print("ERROR in function get_mask: level given is not within depth range 0 <= level <= "
                      + str(cube_data.shape[0]))
                sys.exit()
        elif len(level) == 2:
            if level[1] > cube_data.shape[0] or level[0] < 0:
                print("ERROR in function get_mask: level given is not within depth range 0 <= level <= "
                      + str(cube_data.shape[0]))
                sys.exit()

        if len(cube_shape) == 3:
            if len(level) == 2:
                tiles = (level[1] - level[0] + 1, 1, 1)
            else:
                tiles = (1, 1, 1)
        elif len(cube_shape) == 4:  # depth is included
            if len(level) == 2:
                tiles = (level[1] - level[0] + 1, cube_shape[1], 1, 1)
            else:
                tiles = (1, cube_shape[1], 1, 1)

        mask_arr = np.tile(m, tiles)

        print("function get_mask: Mask successfully constructed from mask file.")

        return mask_arr, level


    def calculate_areas(self, cube, lat_name, lon_name):
        """
        Calculate areas of grid boxes of the cube given
        :param cube: iris.cube
        :param lat_name: name of latitude coordinate, string
        :param lon_name: name of longitude cootdinate, string
        :return: saves the areas in a netcdf file
        """
        cube.coord(lat_name).guess_bounds()
        cube.coord(lon_name).guess_bounds()

        areas_arr = None
        if len(cube.data.shape) == 3:
            areas_arr = (np.pi / 180) * np.asarray(iris.analysis.cartography.area_weights(cube[0]))
        elif len(cube.data.shape) == 4:
            areas_arr = (np.pi / 180) * np.asarray(iris.analysis.cartography.area_weights(cube[0, 0]))

        # Write areas to netcdf
        file_area = os.path.join(directories.ANALYSIS, "areas.nc")
        dataset = Dataset(file_area, 'w', format='NETCDF4_CLASSIC')

        # Create dimensions
        lats = cube.coord(lat_name).points
        lons = cube.coord(lon_name).points
        dataset.createDimension('lat', len(lats))
        dataset.createDimension('lon', len(lons))

        # Create variables
        latitudes = dataset.createVariable('latitude', np.float32, ('lat',))
        longitudes = dataset.createVariable('longitude', np.float32, ('lon',))
        areas = dataset.createVariable('areas', np.float32, ('lat', 'lon'))

        # Global Attributes
        dataset.description = 'Calculated areas of grid boxes'
        dataset.history = 'Created by Adanna Akwataghibe (aa14415@ic.ac.uk) on ' + time.ctime(time.time())

        # Variable Attributes
        latitudes.units = 'degree_north'
        longitudes.units = 'degree_east'
        areas.units = 'metres squared'

        # Writing data
        latitudes[:] = lats
        longitudes[:] = lons
        areas[:] = areas_arr

        # Close dataset
        dataset.close()

        print("Areas file created in " + directories.ANALYSIS + " folder.")


    def regrid_from_file(self, points_sample_grid):
        """
        Get latitude and longitude points from netcdf file
        :param points_sample_grid: netcdf file name
        :return: latitude, longitude, arrays
        """

        # Get cube from file
        cubes =  iris.load(os.path.join(directories.DATA, points_sample_grid))
        cube = cubes[0]

        # Get latiude and longitude
        lat_name, lon_name = None, None
        # Get names of coordinates
        for c in cube.coords():
            dd = c.name()
            not_bnds = 'bound' not in dd and 'bnd' not in dd
            if not_bnds:
                if dd[0].lower() == 'y' or 'lat' in dd.lower() or 'latitude' in dd.lower():
                    lat_name = dd
                if dd[0].lower() == 'x' or 'lon' in dd.lower() or 'longitude' in dd.lower():
                    lon_name = dd

        # Get points
        return cube.coord(lat_name).points, cube.coord(lon_name).points



    def extract_parallel(self, till_start, till_end, dr, dim_coords, set_day, min_yr, num_leap_years, const_lon_name, const_lat_name, const_time_name, lon_name,
                         lat_name, time_name, mask_set, level, mask_arr, sample, sp, nc_true, regrid_lats,
                         regrid_lons, lons_points, lats_points, a_ens_files):
        """
        Helper function for parallelisation
        """
        datasets = xr.open_mfdataset(a_ens_files)

        saved, orig_saved = {}, None
        if sample or self.grid:
            orig_saved = {}

        for var in self.variables:
            # Check if variable is in the dataset
            if var not in datasets:  # throw an error and stop
                print("ERROR in function extract_file : Variable %s not found in files" % var)
                sys.exit()

            # Select time period in dataset
            try:
                if 'time' in datasets[var].coords:
                    if not set_day:
                        if not self.monthly:
                            days = datasets.coords['time'].dt.dayofyear
                            num_leap_years = len(days[days == 366])
                            # Get exact indices of start and end date
                        till_start, till_end = get_diff_start_end(self.start_date, self.end_date,
                                                                  min_yr=min_yr, monthly=self.monthly,
                                                                  num_leap_year_input=num_leap_years)
                        set_day = True
                    time_selected = datasets[var][dict(time=slice(till_start, till_end))]
                elif 't' in datasets[var].coords:
                    if not set_day:
                        if not self.monthly:
                            days = datasets.coords['t'].dt.dayofyear
                            num_leap_years = len(days[days == 366])
                            # Get exact indices of start and end date
                        till_start, till_end = get_diff_start_end(self.start_date, self.end_date, min_yr=min_yr,
                                                                  monthly=self.monthly,
                                                                  num_leap_year_input=num_leap_years)
                        set_day = True
                    time_selected = datasets[var][dict(t=slice(till_start, till_end))]
            except Exception as err:
                print("Exception thrown in function extract_data", str(err))
                time_selected = datasets[var]
            # Convert to cube
            cube = time_selected.to_iris()

            if not dr:
                # Get names of longitude and latitude variables
                for c in cube.coords():
                    dd = c.var_name
                    not_bnds = 'bound' not in dd and 'bnd' not in dd
                    if not_bnds:
                        if dd.lower() == 't' or dd.lower() == 'time':
                            time_name = c.name()
                            dim_coords['time'] = time_name
                        elif dd[0].lower() == 'y' or 'lat' in dd.lower() or 'latitude' in dd.lower():
                            lat_name = c.name()
                            dim_coords['latitude'] = lat_name
                        elif dd[0].lower() == 'x' or 'lon' in dd.lower() or 'longitude' in dd.lower():
                            lon_name = c.name()
                            dim_coords['longitude'] = lon_name
                        else:  # 3D data
                            depth_name = c.name()
                            dim_coords['depth'] = depth_name
                dr = True

            # Rename coords
            try:
                if lon_name is not None and lat_name is not None:
                    coord = cube.coord(lon_name)
                    coord.rename(const_lon_name)
                    coord = cube.coord(lat_name)
                    coord.rename(const_lat_name)
            except Exception:
                print("WARNING in function extract_data: No longitude and latitude coordinate found.")
            try:
                if time_name is not None:
                    coord = cube.coord(time_name)
                    coord.rename(const_time_name)
            except Exception:
                print("WARNING in function extract_data: No time coordinate found.")

            # Calculate areas
            if self.calc_areas:
                self.calculate_areas(cube, const_lat_name, const_lon_name)
                self.calc_areas = False

            # Get longitude and latitude
            lons, lats = None, None
            try:
                lons, lats = cube.coord(const_lon_name).points, cube.coord(const_lat_name).points
            except Exception:
                pass

            # Centre to new longitude centre
            if self.lon_centre is not None:
                # Move longitude centre
                lon_low = self.lon_centre - 180
                lon_high = self.lon_centre + 180.000000000001
                lons = iris.analysis.cartography.wrap_lons(lons, lon_low, lon_high - lon_low)
                count, _ = shift_by_index(lons, self.lon_centre)
                lons.sort()

                # Get axis number
                coords = cube.dim_coords
                axis = 0
                for j in coords:
                    if j.name() == const_lon_name:
                        break
                    axis += 1

                # Move map centre and replace data
                shifted_cube = np.roll(cube.data, count, axis=axis)
                cube.data = shifted_cube
                dim_coord = cube.coord(const_lon_name)
                centred_coord = iris.coords.DimCoord(lons, standard_name=dim_coord.standard_name,
                                                     units=dim_coord.units,
                                                     long_name=dim_coord.long_name, var_name=dim_coord.var_name,
                                                     attributes=dim_coord.attributes, bounds=dim_coord.bounds)
                # Replace coordinate
                cube.replace_coord(centred_coord)

            # Reduced cube
            reduced_cube = None

            # Masking
            if self.maskfile is not None:
                # Swap dimensions if not (depth, time, lat, lon)
                c_n = [c.name() for c in cube.coords(dim_coords=True)]
                indices = []
                if len(c_n) == 4:
                    for dim in c_n:
                        if dim == depth_name:
                            indices.append(0)
                        if dim == time_name:
                            indices.append(1)
                        if dim == const_lat_name:
                            indices.append(2)
                        if dim == const_lon_name:
                            indices.append(3)
                elif len(c_n) == 3:
                    for dim in c_n:
                        if dim == time_name:
                            indices.append(0)
                        if dim == const_lat_name:
                            indices.append(1)
                        if dim == const_lon_name:
                            indices.append(2)

                if indices != [0, 1, 2, 3]:
                    cube.transpose(indices)

                # Get mask array and update cube
                if not mask_set:
                    mask_arr, level = self.get_mask(self.maskfile, cube.data, lons, lats)
                    mask_set = True

                if len(level) == 2:
                    # Get specific cube.data
                    level_low, level_high = level[0] - 1, level[1]
                    mask_cube = np.ma.array(cube.data[level_low:level_high], mask=~mask_arr)
                else:  # only 1 level given
                    mask_cube = np.ma.array(cube.data[level[0] - 1], mask=~mask_arr)

                if cube.data.shape != mask_cube.shape:
                    if len(level) == 2:
                        reduced_cube = cube[level_low:level_high]
                    else:
                        reduced_cube = cube[level[0] - 1]
                    reduced_cube.data = mask_cube
                else:
                    cube.data = mask_cube
                    # Transpose back
                    cube.transpose(indices)

            if sample or self.grid:
                # Save interpolated values
                found_grid = None

                # Check lat and lon are in grid
                if sp and (lons_points is not None or lats_points is not None):
                    if not np.all(min(lons) <= lons_points) and np.all(lons_points <= max(lons)) and \
                            not np.all(min(lats) <= lats_points) and np.all(lats_points <= max(lats)):
                        print("ERROR: extract_data function: latitudes and longitudes are outside grid.")
                        sys.exit()
                elif not sp:
                    if not (min(lons) <= self.lon <= max(lons)) or not (min(lats) <= self.lat <= max(lats)):
                        print("ERROR: extract_data function: latitude and longitude are outside grid.")
                        sys.exit()

                # Construct lat and lon to give to iris cube interpolate function
                # If more points are given
                samples = None
                if sp and (lons_points is not None or lats_points is not None):
                    samples = [(const_lat_name, lats_points), (const_lon_name, lons_points)]
                elif sp and nc_true:
                    samples = [(const_lat_name, regrid_lats), (const_lon_name, regrid_lons)]
                else:
                    samples = [(const_lat_name, self.lat), (const_lon_name, self.lon)]

                if self.grid:
                    # Uses nearest neighbour interpolation and takes into account spherical distance
                    if reduced_cube is None:
                        found_grid = cube.interpolate(samples, iris.analysis.Nearest())
                    else:
                        found_grid = reduced_cube.interpolate(samples, iris.analysis.Nearest())

                    # If in NaN grid space, tell user
                    if np.isnan(found_grid.data).any():
                        warnings.warn("NaN values found in interpolated grids.")
                if sample:
                    if reduced_cube is None:
                        found_grid = cube.interpolate(samples, iris.analysis.Linear())
                    else:
                        found_grid = reduced_cube.interpolate(samples, iris.analysis.Linear())

                saved[var] = found_grid

                # Save original cube
                if reduced_cube is None:
                    orig_saved[var] = cube
                else:
                    orig_saved[var] = reduced_cube

            else:
                if reduced_cube is None:
                    saved[var] = cube
                else:
                    saved[var] = reduced_cube

        return saved, orig_saved, dim_coords


    def extract_data(self):
        """
        Extracts the data given by the user and stores them
        :param algae_type: name of prefix of filename to look into
        :param variables: list of variables to extract from files e.g. ['temp', 'sal']
        :param start_date: start date given to user
        :param end_date: end date given to user
        :param num_ens: number of ensembles, int
        :param monthly: data is stored in monthly increments (time = 12) else assumed (time = 365)
        :param lat: latitude, set if grid or sample point, floats
        :param lon: longitude, set if grid or sample point, floats
        :param grid: set if grid point is given
        :param points_sample_grid: file (txt or nc) that contains sample or grid points
        :param lon_centre: longitude center , float
        :param maskfile: file that contains mask data, string
        :param calc_areas: set if areas of grid boxes should be calculated, boolean
        :return: dictionary storing arrays or list of arrays:
                e.g. if only one file inputted and variables = ['temp', 'sal'], then
                    dict = {'temp': [...], 'sal': [...]
                    if multiple files inputted and variables = ['temp', 'sal'], then
                    dict = {'temp': [ [..], [..], ..], 'sal': [ [..], [..], ..]
                time_name: name of time variable in file
                ens_files: the data nc files that will be used when writing output in netcdf file
                abs_files: same as ens_files except files are in absolute path
                orig_saved: dict above with full map - this is filled if lat/lon are not None
        """
        # Assertions
        assert self.variables is not None
        assert check_list_date(self.start_date) and check_list_date(self.end_date)

        # If points of sample or grid, then open file
        sp, lons_points, lats_points, nc_true = False, None, None, False
        if self.points_sample_grid is not None:
            lons_points, lats_points, nc_true = get_sample_grid_points(self.points_sample_grid)
            sp = True

        # Check if sample point
        sample = False
        if self.lat and self.lon and not self.grid:
            sample = True
        if sp and not self.grid:
            sample = True

        # Only once, dim_coords
        dr, dim_coords = False, {}

        # Get day, month and year
        day_s, mon_s, yr_s = self.start_date[0], self.start_date[1], self.start_date[2]
        day_e, mon_e, yr_e = self.end_date[0], self.end_date[1], self.end_date[2]

        # Get path
        path = directories.DATA

        # Get files and min and maximum year
        files, min_yr, max_yr = get_files_time_period(self.algae_type, yr_s, yr_e)

        # Save list of dictionaries - each dict in the list is an ensemble
        saved = [{} for _ in range(self.num_ens)]
        # Save original (without sample or grid point)
        orig_saved = None
        if sample or self.grid:
            orig_saved = [{} for _ in range(self.num_ens)]

        # Get files in absolute path saved in ensemble groups
        files = [os.path.join(path, file) for file in files]
        ens_files = [[] for _ in range(self.num_ens)]
        abs_files = [[] for _ in range(self.num_ens)]

        for i in range(len(files)):
            ens_indx = ens_to_indx(get_ens_num(files[i]), self.num_ens)
            if ens_indx == -1:  # Ensemble number > Number of ensembles
                continue
            # save in ens_files
            ens_files[ens_indx].append(files[i])
            # Get absolute path
            joined = os.path.abspath(files[i])

            if self.test:  # make changes specific to my PC
                joined = joined.replace("Adanna Akwataghibe", "Adanna")
            # save in ens_files
            abs_files[ens_indx].append(joined)

        # Get corrct start and end positions
        till_start, till_end, set_day, num_leap_years = None, None, False, None

        # Save names of latitude, longitude, depth or time
        time_name, lat_name, lon_name = None, None, None
        # Set if data is 4D
        depth_name = None
        # Change name of latitude and longitude to regular names
        const_lon_name, const_lat_name, const_time_name = "longitude", "latitude", "time"
        # If mask set, then compute mask array
        # Used to compute mask array and longitude centre only once
        mask_set, mask_arr, level = False, None, None

        # If sample/grid point and regrid file
        regrid_lats, regrid_lons = None, None
        if sp and nc_true:
            regrid_lats, regrid_lons = self.regrid_from_file(self.points_sample_grid)

        pool = Pool(processes=parallel_settings.NUM_PROCESSORS)
        # Go through variables in each ensemble

        func = partial(self.extract_parallel, till_start, till_end, dr, dim_coords, set_day, min_yr, num_leap_years, const_lon_name, const_lat_name, const_time_name,
                       lon_name, lat_name, time_name, mask_set, level, mask_arr, sample, sp, nc_true, regrid_lats,
                         regrid_lons, lons_points, lats_points)

        r = pool.map(func, ens_files)
        pool.close()
        pool.join()

        # Get saved and orig_saved
        saved = [tup[0] for tup in r]
        orig_saved = [tup[1] for tup in r]
        dim_coords = r[0][2]

        print("function extract_data: Data successfully extracted from files.")

        return saved, ens_files, abs_files, orig_saved, dim_coords


