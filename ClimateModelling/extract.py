import sys
import time
from scipy import interpolate
from utils import *
import matplotlib.pyplot as plt
from netCDF4 import Dataset
from nco import Nco
import directories


def extract_data(algae_type, variables, start_date, end_date, num_ens, monthly=False, lat=None, lon=None,
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


def create_histogram(list_ens, units, start_date, end_date, nan_values, monthly=False, save_out=None,
                     cov=None, sel=None, plot=False):
    """
    Analysis the data given - in this case it computes the histogram (assumes grid/sample point)
    :param list_ens: the list of ensembles (dicts) containing the data of the climate variables
    :param units: the units matching to each variable
    :param start_date: start date given to user
    :param end_date: end date given to user
    :param nan_values: the missing values in the variables data array
    :param monthly: data is stored in monthly increments (time = 12) else assumed (time = 365)
    :param save_out: if set, then save output of histogram/ rimeseries
    :param cov: if set, then perform covariance analysis
    :param sel: selection option for bin siz, default is fd - Freedman Diaconis Estimator
    :param plot: if plot is true, then shows plot of histogram
    :return: None
    """

    if not sel:
        sel = 'fd'
    fig, axs = plt.subplots(len(list_ens), len(list_ens[0]), squeeze=False)
    time_str = "daily"
    if monthly:
        time_str = "monthly"
    fig.suptitle("Variables " + str(list(list_ens[0])) + " measured " + time_str + " between " + str(start_date[0]) +
                 "-" + str(start_date[1]) + "-" + str(start_date[2]) + " and " + str(end_date[0]) + "-" +
                 str(end_date[1]) + "-" + str(end_date[2]) + " using the E2S2M climate model")
    a, e = 0, 0
    for dict_ in list_ens:
        axs[e, a].set_title("Ensemble " + str(e))
        for d in dict_:
            ens = dict_[d].flatten()
            indices = np.argwhere(np.isclose(ens, nan_values[d]))
            ens = np.delete(ens, indices)
            hist, bin_edges = np.histogram(ens, bins=sel)
            print(ens)
            if plot and not cov:
                axs[e, a].hist(ens, bins=sel)
                axs[e, a].set_ylabel("Frequency")
                axs[e, a].set_xlabel(d + ' (' + units[d] + ')')
                # a += 1

            if plot and cov:  # If covariance between 2 variables, plot a 2d histogram
                axs[a].hist2d(ens, bins=sel)
        e += 1

    plt.show()

    return None


def create_timeseries(list_ens, units, start_date, end_date, monthly=False, save_out=None, cov=None):
    """
    Analysis the data given - in this case it computes the timeseries (assumes grid/sample point)
    :param list_ens: the list of ensembles (dicts) containing the data of the climate variables
    :param units: the units matching to each variable
    :param start_date and end_date: extract data within this time frame
    :param monthly: data is stored in monthly increments (time = 12) else assumed (time = 365)
    :param save_out: if set, then save output of histogram/ rimeseries
    :param cov: if set, then perform covariance analysis
    :return: None
    """
    return None


def compute_average(list_ens, nan_values):
    """
    Analysis the data given - in this case it computes the mean
    :param list_ens: the list of ensembles (dicts) containing the data of the climate variables
    :param nan_values: missing values in data set
    :return:
        ens_means: list of averages of the different ensembles
    """

    # Holds the means for each ensemble
    ens_means = []
    for dict_ in list_ens:
        # Save the mean of each variable in a dict of list
        means = {}
        # Calculate the mean of each variable in the dictionary given
        for d in dict_:
            # Select the parts of the data within timeframe
            mean = np.mean(dict_[d], axis=0)
            # Replace values close to nan values to actual nan values
            if mean.shape:
                mean[np.isclose(mean, nan_values[d], rtol=1)] = nan_values[d]
            # Save mean for variable
            means[d] = mean
        ens_means.append(means)

    return ens_means


def write_means_to_netcdf_file(files, ens_means, variables, start_date, end_date, argv, test=False):
    """
    Write means computed in netcdf files
    :param files: initial files
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

    # Start and end date string
    start_end_str = str(start_date[0]) + "-" + str(start_date[1]) + "-" + str(start_date[2]) + " and " + \
                    str(end_date[0]) + "-" +str(end_date[1]) + "-" + str(end_date[2])

    # Get path
    path = directories.CLIMATE_DATA

    # Get normal and files in absolute path saved in ensemble groups
    ens_files = [[] for _ in range(len(ens_means))]
    abs_files = [[] for _ in range(len(ens_means))]

    # Get absolute path of each file
    for i in range(len(files)):
        # Get file ensemble number index
        ens_indx = ens_to_indx(get_ens_num(files[i]))
        # save in ens_files
        ens_files[ens_indx].append(files[i])
        # Get absolute path
        joined = os.path.abspath(os.path.join(path, files[i]))

        if test:
            joined = joined.replace("Adanna Akwataghibe", "Adanna")
        # save in ens_files
        abs_files[ens_indx].append(joined)

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
            output_file = ens_files[i][0][:-3] + '_' + str(end_date[2]) + '.nc'

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
                mean_var[:] = ens_means[i][var][:]
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

    print("Mean ensemble files created in " + os.path.join(directories.ANALYSIS, directories.MEANS) + "folder.")


def plot_graph(file):
    """
    Plot the data given in file
    :param file: file from analysis
    :return: graph plot in output + saves as png file
    """

    # Open file
    dataset = Dataset(file, 'r')

    d = np.array(dataset.variables['mean_air_temperature'])

    fig, ax = plt.subplots()
    ax.imshow(d)

    png_file = file.rstrip('nc') + 'png'
    fig.savefig(png_file)

    print("Image is saved in the " + directories.ANALYSIS + " folder as a png file.")

    return None
