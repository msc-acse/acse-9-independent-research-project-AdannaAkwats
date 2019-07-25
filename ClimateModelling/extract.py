from scipy import interpolate
from utils import *
from netCDF4 import Dataset
import directories
import xarray as xr
import iris
import warnings
import cartopy.crs as ccrs
from PIL import Image, ImageDraw
from functools import reduce
import operator
import math


def extract_data(algae_type, variables, start_date, end_date, num_ens, monthly=False, lat=None, lon=None,
                 grid=None, lon_bounds=None, test=True):
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
    :param lon_bounds: longitude center range, tuple
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

    # Get day, month and year
    day_s, mon_s, yr_s = start_date[0], start_date[1], start_date[2]
    day_e, mon_e, yr_e = end_date[0], end_date[1], end_date[2]

    # Get path
    path = directories.CLIMATE_DATA

    # Get files and min and maximum year
    files, min_yr, max_yr = get_files_time_period(algae_type, yr_s, yr_e)

    # Save list of dictionaries - each dict in the list is an ensemble
    saved = [{} for _ in range(num_ens)]
    # Save original (without sample or grid point)
    orig_saved = None
    if sample or grid:
        orig_saved = [{} for _ in range(num_ens)]

    # Get files in absolute path saved in ensemble groups
    files = [os.path.join(path, file) for file in files]
    ens_files = [[] for _ in range(num_ens)]
    abs_files = [[] for _ in range(num_ens)]

    for i in range(len(files)):
        ens_indx = ens_to_indx(get_ens_num(files[i]))
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

    # Transform the lon lat point into rotated coordinates
    target_xy = None
    if sample or grid:
        pp_cs = iris.coord_systems.GeogCS(iris.fileformats.pp.EARTH_RADIUS)
        rot_pole = pp_cs.as_cartopy_crs()
        target_xy = rot_pole.transform_point(lon, lat, ccrs.Geodetic())

    # Go through variables in each ensemble
    for i in range(num_ens):
        datasets = xr.open_mfdataset(ens_files[i])
        for var in variables:
            # Check if variable is in the dataset
            if var not in datasets:  # throw an error and stop
                print("Error in function extract_file : Variable %s not found in files" % var)
                sys.exit()

            # Select time period in dataset
            time_selected = datasets[var][dict(time=slice(till_start, till_end))]
            # Convert to cube
            cube = time_selected.to_iris()

            # Get names of longitude and latitude variables
            # Get list of coordinate names
            coord_names = [coord.name() for coord in cube.coords()]
            # Get names of coordinates
            for dd in list(coord_names):
                if dd == 't' or dd == 'time':
                    time_name = dd
                elif dd[0].lower() == 'y' or 'lat' in dd.lower() or 'latitude' in dd.lower():
                    lat_name = dd
                elif dd[0].lower() == 'x' or 'lon' in dd.lower() or 'longitude' in dd.lower():
                    lon_name = dd
                # Check if 3D data
                else:
                    if not (dd == 'bounds' or dd == 'bnds'):
                        depth_name = dd

            # Change name of latitude and longitude to regular names
            const_lon_name, const_lat_name, const_time_name = "longitude", "latitude", "time"
            coord = cube.coord(lon_name)
            coord.rename(const_lon_name)
            coord = cube.coord(lat_name)
            coord.rename(const_lat_name)
            coord = cube.coord(time_name)
            coord.rename(const_time_name)

            # Update lon and lat name
            lon_name, lat_name = const_lon_name, const_lat_name

            # Centre to new longitude bounds
            lons = cube.coord(lon_name).points
            if lon_bounds:
                lons = iris.analysis.cartography.wrap_lons(lons, lon_bounds[0], lon_bounds[1] - lon_bounds[0])
            if sample or grid:
                # Save interpolated values
                found_grid = None

                # Check lat and lon are in grid
                lats = cube.coord(lat_name).points
                if not (min(lons) <= lon <= max(lons)) or not (min(lats) <= lat <= max(lats)):
                    print("Error: extract_data function: latitude and longitude are outside grid.")
                    sys.exit()

                # Construct lat and lon to give tp iris cube interpolate function
                samples = [(lat_name, target_xy[1]), (lon_name, target_xy[0])]
                if grid:
                    # Uses nearest neighbour interpolation and takes into account spherical distance
                    found_grid = cube.interpolate(samples, iris.analysis.Nearest())
                    # found_grid = res.data

                    # If in NaN grid space, tell user
                    if np.isnan(found_grid.data).any():
                        warnings.warn("NaN values found in interpolated grids.")
                if sample:
                    found_grid = cube.interpolate(samples, iris.analysis.Linear())

                # TODO: what if variable = Nan ???
                saved[i][var] = found_grid

                # If longitude centre given, then change current coord data to centred data
                if lon_bounds:
                    dim_coord = cube.coord(lon_name)
                    centred_coord = iris.coords.AuxCoord(lons, standard_name=dim_coord.standard_name,
                                                         units=dim_coord.units,
                                                         long_name=dim_coord.long_name, var_name=dim_coord.var_name,
                                                         attributes=dim_coord.attributes, bounds=dim_coord.bounds)
                    # Replace coordinate
                    cube.replace_coord(centred_coord)
                # Save original cube
                orig_saved[i][var] = cube

            else:
                # If longitude centre given, then change current coord data to centred data
                if lon_bounds:
                    dim_coord = cube.coord(lon_name)
                    centred_coord = iris.coords.AuxCoord(lons, standard_name=dim_coord.standard_name,
                                                             units=dim_coord.units,
                                                             long_name=dim_coord.long_name, var_name=dim_coord.var_name,
                                                             attributes=dim_coord.attributes, bounds=dim_coord.bounds)
                    # Replace coordinate
                    cube.replace_coord(centred_coord)

                saved[i][var] = cube

    print("function extract_data: Data successfully extracted from files.")

    return saved, ens_files, abs_files, orig_saved


def get_mask(list_ens, maskfile):
    """
    Apply mask on arrays of ensemble data
    :param list_ens: the list of ensembles (dicts) containing the data of the climate variables
    :param maskfile: file containing the polygons array of mask to go over grid, string
    :param plot: shows example plot of the mask over a the first variable and a random time period
    :return: 2D mask array
    """
    # Assertions
    assert list_ens is not None

    if maskfile is None or not maskfile:
        return None, None

    # Get polygons from mask file
    polygons = get_polygons(maskfile)
    points_x, points_y = [], []

    # if only one polygon
    if not is_nested_list(polygons):
        # Cast nested list as nested list of tuples
        polygons = [tuple(l) for l in polygons]
        # Change to 3D polygon
        polygons = [polygons]

    # if multiple polygons
    # Go through each polygon in the list
    for i in range(len(polygons)):
        # Sort coordinates of polygons to be in clockwise direction
        center = tuple(map(operator.truediv, reduce(lambda x, y: map(operator.add, x, y), polygons[i]),
                           [len(polygons[i])] * 2))
        sorted_poly = sorted(polygons[i], key=lambda coord: (-135 - math.degrees(
            math.atan2(*tuple(map(operator.sub, coord, center))[::-1]))) % 360)
        polygons[i] = sorted_poly

        xs, ys = [], []
        for j in range(len(sorted_poly)):
            xs.append(sorted_poly[j][0])
            ys.append(sorted_poly[j][1])

        points_x.append(xs)
        points_y.append(ys)

    print("function get_mask: Mask successfully constructed from mask file.")

    return points_x, points_y


def get_mask2(list_ens, maskfile):
    """
    Apply mask on arrays of ensemble data
    :param list_ens: the list of ensembles (dicts) containing the data of the climate variables
    :param maskfile: file containing the polygons array of mask to go over grid, string
    :param plot: shows example plot of the mask over a the first variable and a random time period
    :return: 2D mask array
    """
    # Assertions
    assert list_ens is not None

    if maskfile is None or not maskfile:
        return None

    # Get latitude and longitude - use ad height and width
    first_key = list(list_ens[0])[0]
    dims = list_ens[0][first_key].data[0].shape
    width, height = dims[1], dims[0]
    # Get polygons from mask file
    polygons = get_polygons(maskfile)
    shape = np.asarray(polygons).shape
    # Cast nested list as nested list of tuples
    polygons = [tuple(l) for l in polygons]

    # if only one polygon
    if len(shape) == 2:
        # Change to 3D polygon
        polygons = [polygons]
        # Update shape accordingly
        shape = np.asarray(polygons).shape

    # if multiple polygons
    if len(shape) == 3:
        # Create new image with latitude and longitude as height and width respectively
        img = Image.new('L', (width, height), 0)
        # Go through each polygon in the list
        for i in range(len(polygons)):
            # Sort coordinates of polygons to be in clockwise direction
            center = tuple(map(operator.truediv, reduce(lambda x, y: map(operator.add, x, y), polygons[i]),
                               [len(polygons[i])] * 2))
            sorted_poly = sorted(polygons[i], key=lambda coord: (-135 - math.degrees(
                math.atan2(*tuple(map(operator.sub, coord, center))[::-1]))) % 360)
            polygons[i] = sorted_poly
            # Draw polygon on image
            ImageDraw.Draw(img).polygon(polygons[i], fill=i+1)
        # Get mask from image
        mask = np.array(img)

        import matplotlib.pyplot as plt
        plt.imshow(mask)

        print("function get_mask: Mask successfully constructed from mask file.")

        return mask

    else:
        print("Error in function apply_mask : Arrays is not 2D or 3D.")
        sys.exit()


def extract_data2(algae_type, variables, start_date, end_date, num_ens, monthly=False, lat=None, lon=None,
                 grid=False, mask=None):
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
    :param mask: set if mask file is given, file containing the boolean array of mask to go over grid, string
    :return: dictionary storing arrays or list of arrays:
            e.g. if only one file inputted and variables = ['temp', 'sal'], then
                dict = {'temp': [...], 'sal': [...]
                if multiple files inputted and variables = ['temp', 'sal'], then
                dict = {'temp': [ [..], [..], ..], 'sal': [ [..], [..], ..]
            units: units of variables
            files: the data nc files that will be used if function write write_averages_to_netcdf_file
    """
    # Make sure data structures are not empty
    assert algae_type is not None
    assert variables is not None

    # Check if sample point
    sample = False
    if lat and lon and not grid:
        sample = True

    # Get day, month and year
    day_s, mon_s, yr_s = start_date[0], start_date[1], start_date[2]
    day_e, mon_e, yr_e = end_date[0], end_date[1], end_date[2]

    # Get path
    path = directories.CLIMATE_DATA

    # Get files and min and maximum year
    files, min_yr, max_yr = get_files_time_period(algae_type, yr_s, yr_e)

    # Save list of dictionaries - each dict in the list is an ensemble
    saved = [{} for _ in range(num_ens)]

    # Save the units of the variables to use later
    save_units = True  # only save in the first for loop
    units = {}

    save_mask = True  # only save mask in the first for loop
    mask_arr = None

    nan_values = {}

    for file in files:
        # For each relevant file, make dataset and get variable data
        dataset = Dataset(os.path.join(path, file), 'r')

        # Get file ensemble number
        ens_num = get_ens_num(file)
        # Get corresponding index in list
        indx = ens_to_indx(ens_num)

        # Grid point: Get size and name of dimensions
        time_size, lat_size, lon_size = None, None, None
        time_name, lat_name, lon_name = 'time', None, None

        # If grid and sample point, save names and size of time, latitude and longitude
        if grid or sample or mask:
            for dd in dataset.dimensions:
                if dd == time_name:
                    time_size = dataset.dimensions[dd].size
                if dd[0].lower() == 'y' or dd[:3].lower() == 'lat':
                    lat_size = dataset.dimensions[dd].size
                    lat_name = dd
                if dd[0].lower() == 'x' or dd[:3].lower() == 'lon':
                    lon_size = dataset.dimensions[dd].size
                    lon_name = dd

        # If mask, then save mask array for only the first loop
        # Check if mask
        if mask and save_mask:
            # open file and save a np array
            try:
                mask_arr = np.loadtxt(mask, usecols=range(lon_size), dtype=np.int)
                save_mask = False
            except IndexError:
                print("Error: extract_data function: mask file does not have correct latitude and longitude.")
                sys.exit()

        # Save the data for each variable
        for var in variables:
            if save_units:  # Save for only the first file
                nan_val = dataset.variables[var].missing_value
                nan_values[var] = nan_val

            ds = np.array(dataset.variables[var])

            # If we have grid or sample point
            if grid or sample:
                # Check that dimensions match : Note this only works with 3d data
                if (time_size, lat_size, lon_size) == ds.shape:
                    # Get lat and lon array
                    lat_arr, lon_arr = np.array(dataset.variables[lat_name]), np.array(dataset.variables[lon_name])
                    if grid:
                        # Get index of closest value to lat and lon in arrays
                        lat_indx, lon_indx = find_nearest(lat_arr, lat), find_nearest(lon_arr, lon)
                        # Get specific grid point in variable
                        ds = ds[:, lat_indx, lon_indx]
                    if sample:
                        ds_ = []
                        for j in range(ds.shape[0]):
                            f = interpolate.interp2d(lon_arr, lat_arr, ds[j])
                            ds_.append(f(lon, lat))
                        # Replace ds
                        ds = np.asarray(ds_).flatten()
                else:
                    print("Error: extract_data function: Dimensions do not match with variables.")
                    sys.exit()

            if save_units:  # Save the units for only the first file
                unit = dataset.variables[var].units
                units[var] = unit

            # Check if variable name is already in dict, if so add to the list in dict
            if var in saved[indx]:
                cur_d = saved[indx].get(var)
                # Concatenate list saved and new list
                ds = np.concatenate((cur_d, ds))

            if mask:
                # use mask to select relevant point in grid
                for i in range(ds.shape[0]):
                    d = np.ma.array(ds[i], mask=mask_arr, fill_value=nan_values[var])
                    ds[i] = d.filled()

            # Save variable name and data in the dict
            saved[indx][var] = ds

            # Close datset
            dataset.close()

        # Do not save units anymore, since we have all units now
        save_units = False

    # Get specific time frame
    till_start, till_end = get_diff_start_end(start_date, end_date, min_yr=min_yr, monthly=monthly)
    # For multiple years in one file, the dicts in saved should have shape of max_yr - min_yr + 1
    # If days are reduced then select time frame
    if (till_end - till_start) != saved[0][variables[0]].shape[0]:
        for var in variables:
            for indx in range(num_ens):
                saved[indx][var] = saved[indx][var][till_start:till_end, :, :]

    return saved, units, files, nan_values


