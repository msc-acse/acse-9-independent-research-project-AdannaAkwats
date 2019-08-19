import iris

"""
This is an script that given an example of a user-given function
In the command line or input file :
file_name = example_function 
func_name = simple_equation
"""

def simple_equation(data):
    """
    Perform simple arithmetic with variables
    :param data: a dictionary of variables and their data
    :return: result of equation, name of results, unit
    """

    # Perform equation wirh temperature and salinity
    result = data['temp_mean'].data * 2 + data['sal_mean'].data * 100
    name = 'result'
    long_name = 'result calculated using simple equation'
    unit = 'K'

    return result, name, long_name, unit


def difference_temp(data):
    """
    Perform simple arithmetic with variables
    :param data: a dictionary of variables and their data
    :return: result of equation, name of results, unit
    """

    # Perform equation wirh temperature and salinity
    result = data['temp_mean'].data[0] - data['temp_mean'].data[5]
    name = 'diff_temp'
    long_name = 'difference of temperature calculated using difference_temp'
    unit = 'K'

    return result, name, long_name, unit
