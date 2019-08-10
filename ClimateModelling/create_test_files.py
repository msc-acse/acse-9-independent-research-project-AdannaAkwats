from netCDF4 import Dataset, num2date, date2num
import time
import numpy as np
from numpy.random import uniform
from datetime import datetime, timedelta


def ex():
    # time - 365 daily , 2000 - 2005
    dataset = Dataset('data/ex_ens101_2000.nc', 'w', format='NETCDF4_CLASSIC')

    # Create datasets
    dataset.createDimension('lat', 10)
    dataset.createDimension('lon', 20)
    dataset.createDimension('time', None)

    times = dataset.createVariable('time', np.float64, ('time',))
    latitudes = dataset.createVariable('latitude', np.float32, ('lat',))
    longitudes = dataset.createVariable('longitude', np.float32, ('lon',))

    # Create the 3D variable
    temp = dataset.createVariable('temp', np.float32, ('time', 'lat', 'lon'), fill_value=-1.e+20)

    # Global Attributes
    dataset.description = 'Example NetCDF file'
    dataset.history = 'Created by Adanna Akwataghibe ' + time.ctime(time.time())
    dataset.source = 'netCDF4 python module tutorial'

    # Variable Attributes
    latitudes.units = 'degree_north'
    longitudes.units = 'degree_east'
    temp.units = 'K'
    time.calendar = 'NOLEAP'
    time.calendar_type = 'NOLEAP'
    times.units = 'days since 2000-01-01 00:00:00'
    temp.missing_value = -1.e+20

    # Fill in lons and lats
    lats = np.linspace(-90, 90, 10)
    lons = np.linspace(-180, 180, 20)
    latitudes[:] = lats
    longitudes[:] = lons

    nlats = len(dataset.dimensions['lat'])
    nlons = len(dataset.dimensions['lon'])
    temp[:, :, :] = uniform(size=(365, nlats, nlons))

    # Fill in times
    # dates = []
    # for n in range(temp.shape[1]):
    #     dates.append(datetime(2001, 3, 1))
    # times[:] = date2num(dates, units=times.units, calendar=times.calendar)
    times[:] = list(range(365))

    dataset.close()


def ex_level():
    # time - 365 daily , 2000 - 2005, depth
    dataset = Dataset('data/ex_level_ens101_2000.nc', 'w', format='NETCDF4_CLASSIC')

    dataset.createDimension('level', 5)
    dataset.createDimension('lat', 73)
    dataset.createDimension('lon', 144)
    dataset.createDimension('time', None)

    times = dataset.createVariable('time', np.float64, ('time',))
    levels = dataset.createVariable('level', np.int32, ('level',))
    latitudes = dataset.createVariable('latitude', np.float32, ('lat',))
    longitudes = dataset.createVariable('longitude', np.float32, ('lon',))

    # Create the actual 4D variable
    temp = dataset.createVariable('temp', np.float32, ('level', 'time', 'lat', 'lon'), fill_value=-1.e+20)

    # Global Attributes
    dataset.description = 'Example NetCDF file'
    dataset.history = 'Created by Adanna Akwataghibe ' + time.ctime(time.time())
    dataset.source = 'netCDF4 python module tutorial'

    # Variable Attributes
    latitudes.units = 'degree_north'
    longitudes.units = 'degree_east'
    temp.units = 'K'
    time.calendar = 'NOLEAP'
    time.calendar_type = 'NOLEAP'
    times.units = 'hours since 2000-01-01 00:00:00'
    times.calendar = 'gregorian'
    temp.missing_value = -1.e+20
    temp._FillValue = -1.e+20
    levels.units = "hPa"

    # Fill in lons and lats
    lats = np.arange(-90, 91, 2.5)
    lons = np.arange(-180, 180, 2.5)
    lvs = np.linspace(900, 1000, 5)
    latitudes[:] = lats
    longitudes[:] = lons
    levels[:] = lvs

    nlats = len(dataset.dimensions['lat'])
    nlons = len(dataset.dimensions['lon'])
    temp[0:5, :, :, :] = uniform(size=(5, 365, nlats, nlons))

    # Fill in times
    # dates = []
    # for n in range(temp.shape[1]):
    #     dates.append(datetime(2001, 3, 1))
    # times[:] = date2num(dates, units=times.units, calendar=times.calendar)
    times[:] = list(range(365))

    dataset.close()


def ex_monthly():
    # time - 12 monthly, 2000 - 2005
    dataset = Dataset('data/ex_monthly_ens101_2000.nc', 'w', format='NETCDF4_CLASSIC')

    lat = dataset.createDimension('lat', 73)
    lon = dataset.createDimension('lon', 144)
    time = dataset.createDimension('time', None)

    times = dataset.createVariable('time', np.float64, ('time',))
    latitudes = dataset.createVariable('latitude', np.float32, ('lat',))
    longitudes = dataset.createVariable('longitude', np.float32, ('lon',))

    # Create the 3D variable
    temp = dataset.createVariable('temp', np.float32,('time', 'lat', 'lon'))


ex()
