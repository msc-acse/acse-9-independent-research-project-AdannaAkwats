"""
Created by Adanna Akwataghibe (Github: AdannaAkwats)
"""
import iris

"""
This is an script that given an example of a user-given function
In the command line or input file :
file_name = example_function 
func_name = temp_rh
"""

def temp_rh(cube):
    """
   Perform simple arithmetic with variables
   :param data: a dictionary of variables and their data
   :return: result of equation, name of results, unit
   """

    # Perform equation wirh temperature and salinity
    result = cube['temp'].data - cube['rh'].data
    name = 'diff_temp_rh'
    long_name = 'difference of temperature and relative humidity'
    unit = 'K'

    return result, name, long_name, unit
