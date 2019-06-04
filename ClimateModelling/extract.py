import os
import time
import numpy as np
import matplotlib.pyplot as plt
from netCDF4 import Dataset


def extract_data(variables, time_period, path="data"):
    """
    Extracts the data given by the user and stores them
    :inputs:
        variables - list of variables to extract from files e.g. ['temp', 'sal']
        time_period - extract data within this time frame
        path - folder where data is stored, this is "data" in our case, user can
                change this if their data is saved in another folder. Note that the
                folder must be in the same directory as this python file
    :return: dictionary storing arrays or list of arrays
            e.g. if only one file inputted and variables = ['temp', 'sal'], then
                dict = {'temp': [...], 'sal': [...]
                if multiple files inputted and variables = ['temp', 'sal'], then
                dict = {'temp': [ [..], [..], ..], 'sal': [ [..], [..], ..]
    """

    # Save dict
    saved = {}

    # Go through each file in the folder 'data'
    filenames = os.listdir(path)

   # TODO: do you go through only one ensemble? e.g only ens101 or multiple?
    for i in range(len(filenames)):
        # Select relevant files
        # TODO

        # For each relevant file, make dataset and get variable data
        dataset = Dataset(path + "/" + filenames[i], 'r')
        print(dataset.variables['air_temperature'].units)

        # Save the data for each variable
        for var in variables:
            ds = np.array(dataset.variables[var])
            if var in saved:
                cur_d = saved.get(var)
                # Concatenate list saved and new list
                ds = np.concatenate((cur_d, ds))

            saved[var] = ds

        break # TODO: remove
            # TODO: Remove later
            # plt.imshow(d[0])
            # plt.show()

    return saved


def analysis(dict):

    # Save the mean of each variable in a list
    means = []

    # Calculate the mean of each variable in the dictionary given
    for d in dict:
        print(saved[d].shape)
        mean = np.mean(saved[d], axis=0)
        means.append(mean)


    # Write to netcdf file, pass in means and variable names
    # write_to_netcdf_file(means, list(dict))

    return None


def write_to_netcdf_file(means, variables):
    # Get latituide and longitude size by means dimensions
    lat_size = means[0].shape[0]
    lon_size = means[0].shape[1]

    # Create netcdf file to write results in
    dataset = Dataset('results/mean.nc', 'w', format='NETCDF4_CLASSIC')
    dataset.createDimension('lat', lat_size)
    dataset.createDimension('lon', lon_size)
    dataset.createDimension('time', None)

    times = dataset.createVariable('time', np.float64, 'time')
    latitudes = dataset.createVariable('latitude', np.float32, 'lat')
    longitudes = dataset.createVariable('longitude', np.float32, 'lon')

    for i in range(len(variables)):
        name = 'mean_'+variables[i]
        mean = dataset.createVariable(name, np.float32, ('time', 'lon', 'lat'))

    # Global attributes
    desc_string = 'Contains the average mean of '
    for i in range(len(variables)):
        desc_string += variables[i]
        if i == len(variables) - 2: # Penultimate element
            desc_string += ' and '
        if i < len(variables) - 2: # Last element
            desc_string += ', '

    dataset.description =  desc_string
    dataset.history = 'Created ' + time.ctime(time.time())
    # dataset.source = 'netCDF4 python module tutorial'

    # Variable Attributes
    latitudes.units = 'degree_north'
    longitudes.units = 'degree_east'
    temp.units = 'K'
    times.units = 'hours since 0001-01-01 00:00:00'
    times.calendar = 'gregorian'


def plot(file):

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


# Checking that everything works
vars = ['air_temperature']
saved = extract_data(vars, None)
analysis(saved)
# TODO: what time peroid will be given - how will it be given (1950 - 1952 or 2 years)

