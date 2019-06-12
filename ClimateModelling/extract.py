import os
import time
import numpy as np
from math import ceil
import re
from utils import *
import matplotlib.pyplot as plt
from netCDF4 import Dataset
import directories


def extract_data(algae_type, variables, start_date, end_date, monthly=False, num_ens=1, lat=None, lon=None, mask=None):
    """
    Extracts the data given by the user and stores them
    :param algae_type: name of prefix of filename to look into
    :param variables: list of variables to extract from files e.g. ['temp', 'sal']
    :param start_date and end_date: extract data within this time frame
    :param monthly: data is stored in monthly increments (time = 12) else assumed (time = 365)
    :param num_ens: number of ensembles, int
    :param lat, lon: latitude and longitude, floats, used if sample or grid point is selected
    :param mask: file containing the boolean array of mask to go over grid, string
    :return: dictionary storing arrays or list of arrays
            e.g. if only one file inputted and variables = ['temp', 'sal'], then
                dict = {'temp': [...], 'sal': [...]
                if multiple files inputted and variables = ['temp', 'sal'], then
                dict = {'temp': [ [..], [..], ..], 'sal': [ [..], [..], ..]
    """

    # Get day, month and year
    day_s, mon_s, yr_s = start_date[0], start_date[1], start_date[2]
    day_e, mon_e, yr_e = end_date[0], end_date[1], end_date[2]

    # Get specific time frame
    till_start, till_end = get_diff_start_end(start_date, end_date, monthly=monthly)

    # Get path and folder
    path = directories.CLIMATE_DATA + '/'
    folder = os.listdir(path)

    # Files should be automatically ordered by year assuming that the format of files is what we expect
    files = []

    # List of years to extract
    years = list(range(yr_s, yr_e + 1))

    # Save lowest and highest year in data for later - only used if multiple years are in the same file
    min_yr = yr_s
    max_yr = yr_e
    mult_yr_seen = False
    # Go through the files in the folder and get the relevant files within the timeframe
    for file in folder:
        if os.path.isfile(os.path.join(path, file)) and file.startswith(algae_type):
            # If file with just one year in it
            if not get_file_two_years(file):
                for year in years:
                    if str(year) in file:
                        files.append(file)
            else: # file has multiple years in it
                mult_yr_seen = True
                fst_yr, snd_yr = get_file_two_years(file)
                # Get files that have data within the years
                if overlaps(fst_yr, snd_yr, yr_s, yr_e):
                    files.append(file)
                    if fst_yr < min_yr:
                        min_yr = fst_yr
                    if snd_yr > max_yr:
                        max_yr = snd_yr

    # Save list of dictionaries - each dict in the list is an ensemble
    saved = [{} for i in range(num_ens)]

    # Save the units of the variables to use later
    save_units = True  # only save in the first for loop
    units = []

    for file in files:
        # For each relevant file, make dataset and get variable data
        dataset = Dataset(os.path.join(path, file), 'r')

        # Get file ensemble number
        ens_num = get_ens_num(file)
        # Get corresponding index in list
        indx = ens_to_indx(ens_num)

        # Save the data for each variable
        for var in variables:
            ds = np.array(dataset.variables[var])

            if save_units: # Save the units for only the first file
                unit = dataset.variables[var].units
                units.append(unit)
            # Check if variable name is already in dict, if so add to the list in dict
            if var in saved[indx]:
                cur_d = saved[indx].get(var)
                # Concatenate list saved and new list
                ds = np.concatenate((cur_d, ds))

            # Save variable name and data in the dict
            saved[indx][var] = ds

            # Dont save units anymore, since we have all units now
            save_units = False

    # TODO get specific time for multiple years
    
    # Get more specific time frame
    # If days are reduced then select time frame
    if (till_end - till_start) != saved[0][variables[0]].shape[0]:
        for var in variables:
            for indx in range(num_ens):
                saved[indx][var] = saved[indx][var][till_start:till_end, :, :]
    return saved, units


def analysis(list_ens, units, start_date, end_date, monthly=False):
    """
    Analysis the data given - in this case it computes the mean
    :param list_ens: the list of ensembles (dicts) containing the data of the climate variables
    :param units: the units matching to each variable
    :param start_date and end_date: extract data within this time frame
    :param monthly: data is stored in monthly increments (time = 12) else assumed (time = 365)
    :return: netcdf file that contains the mean data
    """

    # Holds the means for each ensemble
    ens_means = []
    for dict in list_ens:
        # Save the mean of each variable in a list
        means = []
        # Calculate the mean of each variable in the dictionary given
        for d in dict:
            # Select the parts of the data within timeframe
            # print(dict[d])
            # mean = np.mean(dict[d], axis=0)
            # means.append(mean)
            pass

        ens_means.append(means)

    # Write to netcdf file, pass in means and variable names
    # write_to_netcdf_file(ens_means, list(list_ens[0]), units)


    return None


def write_to_netcdf_file(means, variables, units, time_period="TEMP TODO"):
    """
    Writes to a netcdf file and saves the mean data
    :param means: data created by analysis function
    :param variables: variables of climate model
    :param units: units of climate model
    :param time_period: TODO
    :return: saved netcdf file
    """
    # Get latituide and longitude size by means dimensions
    lat_size = means[0].shape[0]
    lon_size = means[0].shape[1]

    # Create netcdf file to write results in
    path = directories.ANALYSIS + '/analysis.nc'
    dataset = Dataset(path, 'w', format='NETCDF4_CLASSIC')
    dataset.createDimension('lat', lat_size)
    dataset.createDimension('lon', lon_size)
    # dataset.createDimension('time', None)

    # times = dataset.createVariable('time', np.float64, 'time')
    # latitudes = dataset.createVariable('latitude', np.float32, 'lat')
    # longitudes = dataset.createVariable('longitude', np.float32, 'lon')

    for i in range(len(variables)):
        name = 'mean_'+variables[i]
        mean = dataset.createVariable(name, np.float32, ('lat', 'lon'))
        mean.units = units[i]
        mean[:] = means

    # Global attributes
    desc_string = 'Contains the average mean of '
    for i in range(len(variables)):
        desc_string += variables[i]
        if i == len(variables) - 2: # Penultimate element
            desc_string += ' and '
        if i < len(variables) - 2: # Last element
            desc_string += ', '

    dataset.description =  desc_string
    dataset.history = 'Created ' + time.ctime(time.time()) + 'by Adanna Akwataghibe aa14415@ic.ac.uk'
    # dataset.source = 'netCDF4 python module tutorial'

    # Variable Attributes
    # latitudes.units = 'degree_north'
    # longitudes.units = 'degree_east'
    # times.units = 'hours since 0001-01-01 00:00:00'
    # times.calendar = 'gregorian'

    print("New nc file has been saved in the " + directories.ANALYSIS + " folder.")


def plot(file):
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

    print("Image is saved in the " +  directories.ANALYSIS + " folder as a png file.")

    return None

def extract_data2():

    # For the future, dont open specific files like this
    # open all files in data folder for example
    dataset = Dataset('data/sresa1b_ncar_ccsm3-example.nc', 'r')

    print(dataset.dimensions.keys())
    print(dataset.dimensions['lat'])

    print(dataset.dimensions['lon'])

    print(dataset.variables.keys())

    print(dataset.variables['pr'])

    # print(dataset.variables['pr'])

    d = np.array(dataset.variables['pr'])

    # fig, ax = plt.subplots()
    plt.imshow(d[0])
    plt.show()

