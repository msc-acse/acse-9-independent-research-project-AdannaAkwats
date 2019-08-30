"""
Created by Adanna Akwataghibe (Github: AdannaAkwats)
"""
from netCDF4 import Dataset
import time
import numpy as np
from numpy.random import uniform, weibull
import directories
import os

"""
Script that creates example NetCDF files to use when testing
"""


def ex(year, monthly=False):
    """
    Create file with daily data for specific year >= 2000, one variable in file
    :param year: year of file
    :param monthly: if set to true, then creates with monthly data
    """
    assert(year >= 2000)
    # time - 365 daily , 2000 - 2005
    file_name = 'ex1_ens101_' + str(year) + '.nc'
    if monthly:
        file_name = 'ex1monthly_ens101_' + str(year) + '.nc'
    file_name = os.path.join(directories.DATA, file_name)
    dataset = Dataset(file_name, 'w', format='NETCDF4_CLASSIC')

    # Create datasets
    dataset.createDimension('latitude', 10)
    dataset.createDimension('longitude', 20)
    dataset.createDimension('time', None)

    times = dataset.createVariable('time', np.float64, ('time',))
    latitudes = dataset.createVariable('latitude', np.float64, ('latitude',))
    longitudes = dataset.createVariable('longitude', np.float64, ('longitude',))

    # Create the 3D variable
    temp = dataset.createVariable('temp', np.float32, ('time', 'latitude', 'longitude'), fill_value=-1.e+20)

    # Global Attributes
    dataset.description = 'Example NetCDF file'
    dataset.history = 'Created by Adanna Akwataghibe ' + time.ctime(time.time())
    dataset.source = 'netCDF4 python module tutorial'

    # Variable Attributes
    latitudes.units = 'degree_north'
    longitudes.units = 'degree_east'
    longitudes.long_name = 'longitude'
    latitudes.long_name = 'latitude'
    times.calendar = 'NOLEAP'
    times.calendar_type = 'NOLEAP'
    times.units = 'days since 2000-01-01 00:00:00'
    times.long_name = 'time'
    temp.units = 'K'
    temp.long_name = 'temperature'

    # Fill in lons and lats
    lats = np.linspace(-90, 90, 10)
    lons = np.linspace(-180, 180, 20)
    latitudes[:] = lats
    longitudes[:] = lons

    nlats = len(dataset.dimensions['latitude'])
    nlons = len(dataset.dimensions['longitude'])

    if not monthly:
        temp[:, :, :] = uniform(size=(365, nlats, nlons))
    else:
        temp[:, :, :] = uniform(size=(12, nlats, nlons))

    if not monthly:
        times[:] = np.asarray(range(365)) + 365 * (year - 2000)
    else:
        t = np.cumsum([0, 31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30])
        times[:] = t + 365 * (year - 2000)

    print("New file " + file_name + " created in " + directories.DATA + " directory.")

    dataset.close()


def ex_level(year, swap=False):
    """
    Create file with daily data for specific year >= 2000, one variable in file, added level to dimensions
    :param year: year of file
    :param swap: if swap, then change time and level dim
    """
    assert (year >= 2000)
    # time - 365 daily , 2000 - 2005
    file_name = 'ex1level_ens101_' + str(year) + '.nc'
    file_name = os.path.join(directories.DATA, file_name)
    dataset = Dataset(file_name, 'w', format='NETCDF4_CLASSIC')

    # Create datasets
    dataset.createDimension('level', 5)
    dataset.createDimension('latitude', 10)
    dataset.createDimension('longitude', 20)
    dataset.createDimension('time', None)

    levels = dataset.createVariable('level', np.float64, ('level',))
    times = dataset.createVariable('time', np.float64, ('time',))
    latitudes = dataset.createVariable('latitude', np.float64, ('latitude',))
    longitudes = dataset.createVariable('longitude', np.float64, ('longitude',))

    # Create the 3D variable
    if not swap:
        temp = dataset.createVariable('temp', np.float32, ('level', 'time', 'latitude', 'longitude'), fill_value=-1.e+20)
    else:
        temp = dataset.createVariable('temp', np.float32, ('time', 'level', 'latitude', 'longitude'), fill_value=-1.e+20)

    # Global Attributes
    dataset.description = 'Example NetCDF file'
    dataset.history = 'Created by Adanna Akwataghibe ' + time.ctime(time.time())
    dataset.source = 'netCDF4 python module tutorial'

    # Variable Attributes
    latitudes.units = 'degree_north'
    longitudes.units = 'degree_east'
    longitudes.long_name = 'longitude'
    latitudes.long_name = 'latitude'
    times.calendar = 'NOLEAP'
    times.calendar_type = 'NOLEAP'
    times.units = 'days since 2000-01-01 00:00:00'
    times.long_name = 'time'
    temp.units = 'K'
    temp.long_name = 'temperature'
    levels.units = "hPa"

    # Fill in lons and lats
    lats = np.linspace(-90, 90, 10)
    lons = np.linspace(-180, 180, 20)
    latitudes[:] = lats
    longitudes[:] = lons
    levels[:] = [100, 75, 50, 25, 0]

    nlats = len(dataset.dimensions['latitude'])
    nlons = len(dataset.dimensions['longitude'])

    if not swap:
        temp[0:5, :, :, :] = uniform(size=(5, 365, nlats, nlons))
    else:
        temp[:, :, :, :] = uniform(size=(365, 5, nlats, nlons))

    times[:] = np.asarray(range(365)) + 365 * (year - 2000)

    print("New file " + file_name + " created in " + directories.DATA + " directory.")

    dataset.close()



def ex2(year):
    """
    Create file with daily data for specific year >= 2000, two variables in file
    :param year: year of file
    """
    assert (year >= 2000)
    # time - 365 daily , from 2000, 2 variables
    file_name = 'ex2_ens101_' + str(year) + '.nc'
    file_name = os.path.join(directories.DATA, file_name)
    dataset = Dataset(file_name, 'w', format='NETCDF4_CLASSIC')

    # Create datasets
    dataset.createDimension('latitude', 10)
    dataset.createDimension('longitude', 20)
    dataset.createDimension('time', None)

    times = dataset.createVariable('time', np.float64, ('time',))
    latitudes = dataset.createVariable('latitude', np.float64, ('latitude',))
    longitudes = dataset.createVariable('longitude', np.float64, ('longitude',))

    # Create the 3D variable
    temp = dataset.createVariable('temp', np.float32, ('time', 'latitude', 'longitude'), fill_value=-1.e+20)
    sal = dataset.createVariable('sal', np.float32, ('time', 'latitude', 'longitude'), fill_value=-1.e+20)


    # Global Attributes
    dataset.description = 'Example NetCDF file'
    dataset.history = 'Created by Adanna Akwataghibe ' + time.ctime(time.time())
    dataset.source = 'netCDF4 python module tutorial'

    # Variable Attributes
    latitudes.units = 'degree_north'
    longitudes.units = 'degree_east'
    longitudes.long_name = 'longitude'
    latitudes.long_name = 'latitude'
    times.calendar = 'NOLEAP'
    times.calendar_type = 'NOLEAP'
    times.units = 'days since 2000-01-01 00:00:00'
    times.long_name = 'time'
    temp.units = 'K'
    temp.long_name = 'temperature'
    sal.long_name = 'salinity'

    # Fill in lons and lats
    lats = np.linspace(-90, 90, 10)
    lons = np.linspace(-180, 180, 20)
    latitudes[:] = lats
    longitudes[:] = lons

    nlats = len(dataset.dimensions['latitude'])
    nlons = len(dataset.dimensions['longitude'])
    temp[:, :, :] = uniform(size=(365, nlats, nlons))
    sal[:, :, :] = weibull(5, size=(365, nlats, nlons))

    times[:] = np.asarray(range(365)) + 365 * (year - 2000)

    print("New file " + file_name + " created in " + directories.DATA + " directory.")

    dataset.close()


def ex2_level(year):
    """
    Create file with daily data for specific year >= 2000, two variables in file, added level in dimensions
    :param year: year of file
    """
    assert (year >= 2000)
    # time - 365 daily , from 2000, 2 variables
    file_name = 'ex2level_ens101_' + str(year) + '.nc'
    file_name = os.path.join(directories.DATA, file_name)
    dataset = Dataset(file_name, 'w', format='NETCDF4_CLASSIC')

    # Create datasets
    dataset.createDimension('level', 5)
    dataset.createDimension('latitude', 10)
    dataset.createDimension('longitude', 20)
    dataset.createDimension('time', None)

    levels = dataset.createVariable('level', np.float64, ('level',))
    times = dataset.createVariable('time', np.float64, ('time',))
    latitudes = dataset.createVariable('latitude', np.float64, ('latitude',))
    longitudes = dataset.createVariable('longitude', np.float64, ('longitude',))

    # Create the 3D variable
    temp = dataset.createVariable('temp', np.float32, ('level', 'time', 'latitude', 'longitude'), fill_value=-1.e+20)
    sal = dataset.createVariable('sal', np.float32, ('level', 'time', 'latitude', 'longitude'), fill_value=-1.e+20)


    # Global Attributes
    dataset.description = 'Example NetCDF file'
    dataset.history = 'Created by Adanna Akwataghibe ' + time.ctime(time.time())
    dataset.source = 'netCDF4 python module tutorial'

    # Variable Attributes
    latitudes.units = 'degree_north'
    longitudes.units = 'degree_east'
    longitudes.long_name = 'longitude'
    latitudes.long_name = 'latitude'
    times.calendar = 'NOLEAP'
    times.calendar_type = 'NOLEAP'
    times.units = 'days since 2000-01-01 00:00:00'
    times.long_name = 'time'
    temp.missing_value = -1.e+20
    temp.units = 'K'
    temp.long_name = 'temperature'
    sal.long_name = 'salinity'
    levels.units = "hPa"

    # Fill in lons and lats
    lats = np.linspace(-90, 90, 10)
    lons = np.linspace(-180, 180, 20)
    latitudes[:] = lats
    longitudes[:] = lons
    levels[:] = [100, 75, 50, 25, 0]

    nlats = len(dataset.dimensions['latitude'])
    nlons = len(dataset.dimensions['longitude'])
    temp[:, :, :, :] = uniform(size=(5, 365, nlats, nlons))
    sal[:, :, :, :] = weibull(5, size=(5, 365, nlats, nlons))

    times[:] = np.asarray(range(365)) + 365 * (year - 2000)

    print("New file " + file_name + " created in " + directories.DATA + " directory.")

    dataset.close()


# ------------------------------- CREATE FILES -----------------------------------

# ex(2000)
# ex(2001)
# ex(2002)
# ex(2003)
#
# ex_level(2000, swap=True)
# ex2_level(2000)

# ex(2000, monthly=True)
# ex(2001, monthly=True)
# ex(2002, monthly=True)
# ex(2003, monthly=True)

# ex2(2000)
# ex2(2001)
# ex2(2002)
# ex2(2003)

ex2(2000)
ex2(2001)
ex2(2002)
ex2(2003)