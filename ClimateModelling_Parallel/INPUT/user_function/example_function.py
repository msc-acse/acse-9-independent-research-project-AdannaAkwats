"""
Created by Adanna Akwataghibe (Github: AdannaAkwats)
"""
import iris

"""
This is an script that given an example of a user-given function
In the command line or input file :
file_name = example_function 
func_name = max_lon_lat
"""


def max_lon_lat(cube):
    """
    Get maximum longitude and latitude for each time
    :param cube: iris.cube (of ONE variable)
    :return: iris.cube containing time coordinate
    """

    max_cube = cube.collapsed(['longitude', 'latitude'], iris.analysis.MAX)

    # print(max_cube)

    return max_cube
