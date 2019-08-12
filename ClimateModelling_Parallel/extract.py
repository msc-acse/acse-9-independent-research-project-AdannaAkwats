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
import cartopy.crs as ccrs
from functools import reduce
import operator
import math
from shapely.geometry import Polygon, Point
from multiprocessing import Pool
from functools import partial
import parallel_settings


def contains_points(polygon, xp, yp):
    """
    Construct boolean array that shows if coords xp and yp are in polygon
    :param polygon: shapely.geometry.Polygon
    :param xp: list of x coordinates
    :param yp: list of y coordinates
    :return: boolean array
    """

    return np.array([Point(x, y).intersects(polygon) for x, y in zip(xp, yp)],
                    dtype=np.bool)


def get_mask(maskfile, cube_data, lons, lats):
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
        mask = contains_points(poly, x.ravel(), y.ravel())
        full_mask = np.logical_or(full_mask,  mask)

    # Reshape mask from 1D to 2D array
    m = full_mask.reshape(len(lats), len(lons))

    # Increase 2D array to 3D array
    tiles = None
    cube_shape = cube_data.shape
    if level is None:
        level = cube_shape[0]
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


def calculate_areas(cube, lat_name, lon_name):
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
        areas_arr = iris.analysis.cartography.area_weights(cube[0])
    elif len(cube.data.shape) == 4:
        areas_arr = iris.analysis.cartography.area_weights(cube[0, 0])

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


def extract_parallel(variables, till_start, till_end, dr, const_lon_name, const_lat_name, const_time_name,
                 calc_areas, lon_centre, lon_name, lat_name, time_name, maskfile, mask_set, level, mask_arr,
                 sample, grid, lat, lon, target_xy, a_ens_files):
    datasets = xr.open_mfdataset(a_ens_files)

    saved, orig_saved = {}, {}
    for var in variables:
        # Check if variable is in the dataset
        if var not in datasets:  # throw an error and stop
            print("ERROR in function extract_file : Variable %s not found in files" % var)
            sys.exit()

        # Select time period in dataset
        time_selected = datasets[var][dict(time=slice(till_start, till_end))]
        # Convert to cube
        cube = time_selected.to_iris()

        if not dr:
            # Get names of longitude and latitude variables
            # Get list of coordinate names
            coord_names = [coord.name() for coord in cube.coords()]
            # Get names of coordinates
            for dd in list(coord_names):
                if dd.lower() == 't' or dd.lower() == 'time':
                    time_name = dd
                elif dd[0].lower() == 'y' or 'lat' in dd.lower() or 'latitude' in dd.lower():
                    lat_name = dd
                elif dd[0].lower() == 'x' or 'lon' in dd.lower() or 'longitude' in dd.lower():
                    lon_name = dd
                # Check if 3D data
                else:
                    if not (dd == 'bounds' or dd == 'bnds'):
                        depth_name = dd
            dr = True

        # Rename coords
        coord = cube.coord(lon_name)
        coord.rename(const_lon_name)
        coord = cube.coord(lat_name)
        coord.rename(const_lat_name)
        coord = cube.coord(time_name)
        coord.rename(const_time_name)

        # Calculate areas
        if calc_areas:
            calculate_areas(cube, const_lat_name, const_lon_name)
            calc_areas = False

        # Get longitude and latitude
        lons, lats = cube.coord(const_lon_name).points, cube.coord(const_lat_name).points

        # Centre to new longitude centre
        if lon_centre is not None:
            # Move longitude centre
            lon_low = lon_centre - 180
            lon_high = lon_centre + 180.1
            lons = iris.analysis.cartography.wrap_lons(lons, lon_low, lon_high - lon_low)
            count, _ = shift_by_index(lons, lon_centre)
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
        if maskfile is not None:
            # Get mask array and update cube
            if not mask_set:
                mask_arr, level = get_mask(maskfile, cube.data, lons, lats)
                mask_set = True

            # Check level is within limits
            depth_c = cube.coord_dims(depth_name)
            print(depth_c)

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

        if sample or grid:
            # Save interpolated values
            found_grid = None

            # Check lat and lon are in grid
            if not (min(lons) <= lon <= max(lons)) or not (min(lats) <= lat <= max(lats)):
                print("ERROR: extract_data function: latitude and longitude are outside grid.")
                sys.exit()

            # Construct lat and lon to give to iris cube interpolate function
            samples = [(const_lat_name, target_xy[1]), (const_lon_name, target_xy[0])]
            if grid:
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

            # TODO: what if variable = Nan ???
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

    return saved, orig_saved


def extract_data(algae_type, variables, start_date, end_date, num_ens, monthly=False, lat=None, lon=None,
                 grid=None, lon_centre=None, maskfile=None, calc_areas=False, test=True):
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
    assert variables is not None
    assert check_list_date(start_date) and check_list_date(end_date)

    # Check if sample point
    sample = False
    if lat and lon and not grid:
        sample = True

    # Only once
    dr = False

    # Get day, month and year
    day_s, mon_s, yr_s = start_date[0], start_date[1], start_date[2]
    day_e, mon_e, yr_e = end_date[0], end_date[1], end_date[2]

    # Get path
    path = directories.CLIMATE_DATA

    # Get files and min and maximum year
    files, min_yr, max_yr = get_files_time_period(algae_type, yr_s, yr_e)

    # Get files in absolute path saved in ensemble groups
    files = [os.path.join(path, file) for file in files]
    ens_files = [[] for _ in range(num_ens)]
    abs_files = [[] for _ in range(num_ens)]

    for i in range(len(files)):
        ens_indx = ens_to_indx(get_ens_num(files[i]), num_ens)
        if ens_indx == -1:  # Ensemble number > Number of ensembles
            continue
        # save in ens_files
        ens_files[ens_indx].append(files[i])
        # Get absolute path
        joined = os.path.abspath(files[i])

        if test:  # make changes specific to my PC
            joined = joined.replace("Adanna Akwataghibe", "Adanna")
        # save in ens_files
        abs_files[ens_indx].append(joined)

    # Get exact indices of start and end date
    till_start, till_end = get_diff_start_end(start_date, end_date, min_yr=min_yr, monthly=monthly)

    # Save names of latitude, longitude, depth or time
    time_name, lat_name, lon_name = None, None, None
    # Set if data is 4D
    depth_name = None
    # Change name of latitude and longitude to regular names
    const_lon_name, const_lat_name, const_time_name = "longitude", "latitude", "time"
    # If mask set, then compute mask array
    # Used to compute mask array and longitude centre only once
    mask_set, mask_arr, level = False, None, None

    # Transform the lon lat point into rotated coordinates
    target_xy = None
    if sample or grid:
        pp_cs = iris.coord_systems.GeogCS(iris.fileformats.pp.EARTH_RADIUS)
        rot_pole = pp_cs.as_cartopy_crs()
        target_xy = rot_pole.transform_point(lon, lat, ccrs.Geodetic())

    pool = Pool(processes=parallel_settings.NUM_PROCESSORS)
    # Go through variables in each ensemble

    func = partial(extract_parallel, variables, till_start, till_end, dr, const_lon_name, const_lat_name, const_time_name,
                 calc_areas, lon_centre, lon_name, lat_name, time_name, maskfile, mask_set, level, mask_arr,
                 sample, grid, lat, lon, target_xy)

    r = pool.map(func, ens_files)
    pool.close()
    pool.join()

    # Get saved and orig_saved
    saved = [tup[0] for tup in r]
    orig_saved = [tup[1] for tup in r]

    print("function extract_data: Data successfully extracted from files.")

    return saved, ens_files, abs_files, orig_saved


