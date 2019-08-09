import iris

"""
This is an script that given an example of a user-given function
In the command line or input file :
file_name = example_function 
func_name = max_lon_lat
"""

""" INPUT/OUTPUT EXPECTED IN FUNCTION 

INPUT : LIST of DICTIONARIES
-----
- Each element in the LIST represents each ENSEMBLE
- The keys of the DICTIONARY are the VARIABLES
- The values of the DICTIONARY are the iris.cube of the VARIABLE data

EXTRACTING FROM INPUT
-----
E.g. to get variable 'temperature' from ensemble 3
-->  list_ens[3]['temperature']
E.g to get variable 'temperature' from all ensembles
--> temperatures = [le['temperature'] for le in list_ens]

OUTPUT 
-----

"""


def max_lon_lat(cube):
    """
    Get maximum longitude and latitude for each time
    :param cube: iris.cube
    :return: iris.cube containing time coordinate
    """

    max_cube = cube.collapsed(['longitude', 'latitude'], iris.analysis.MAX)

    # print(max_cube)

    return max_cube