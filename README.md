# IRP - A software package for Climate Modelling diagnostics

Written in Python 3 and contains functions that calculate climate model diagnostics. 

## Getting Started
- Download / clone the repository unto your computer
- You must first download and install conda, for example from http://conda.pydata.org/miniconda.html in order to install some dependencies used in the package. 
- Once conda is installed, then you can install **Iris** using conda on any platform with the following command:
    ```
    conda install -c conda-forge iris
    ```
- Install **basemap** for plotting
    ```
    conda install basemap
    ```
- Install the rest of the required python package dependencies with the following command:
```
pip install -r requirements.txt
```

## Usage

### Input and output folders
The `.nc` data files to be analysed should be stored in the folder `DATA`. The `RESULTS` folder will store the analysis, once computed. 
The file names can be changed in the file `directories.py`

There are two ways to call the program: 

### Using an input file

An input file is set up for the user to fill in the required values in ``input.txt``. 
The values to fill in are listed below. The pre-filled values are given as the default. 
```
# REQUIRED ARGUMENTS
# ------------------------------------------------------------------------------
Prefix:
Start date of analysis:
Variables:
Number of ensembles: 1
#
# ------------------------------------------------------------------------------
# OPTIONAL ARGUMENTS
# ------------------------------------------------------------------------------
End date of analysis:
Analysis: mean
Total ensemble stats: True
Plot:
Monthly: False
Grid:
Sample:
Mask file:
Save Output: False
Covary: False
Histogram bin selection:
Longitude centre:
Save extract data: False
User function:
Calculate areas: True
#
```
After filling the values in the file, to run the program simply call: 
```
python main.py
```
#### Example input file
A pre-filled input file is given as an example in `input_example.txt`, and can be run by calling:
```
python main.py -ex
```
This also calculates the mean of `air_temperature` given in `DATA\E1_north_america_ens101_1970.nc`.

### Command line interface (CLI)
```
python main.py -h 
```
shows the list of commands available to use and their descriptions. 
#### Example commands
The file `E1_north_america_ens101_1970.nc` is given as example data in the `DATA` folder. To calculate the mean of variable `air_temperature`, and plot the histogram and map, we would call:
```
python main.py 1970 -pf E1_north_america -v air_temperature -e 1 -a mean -p 1
```

## Assumptions
- Files do not have overlapped data.
- Files given **must** exist for the given time frame
- Daily increments of data, except if the monthly tag is set in the arguments.
- Grids have constant latitude and longitude.
- Files are in the format:
 ```
 "{START OF FILENAME}_ens{NUM}_{YEAR}.nc"
 ```   
  where
  - `{START OF FILENAME}` is the prefix of the file, this can be the algae type.
  - `{NUM}` is the ensemble number
  - `{YEAR}` is the year.

  **OR** if you have multiple years stored in one file then:

  The file names should be in the format:
  ```
   "{START OF FILENAME}_ens{NUM}_{YEAR 1}_{YEAR 2}.nc"
  ```
  where
  - `{START OF FILENAME}` is the prefix of the file, this can be the algae type.
  - `{NUM}` is the ensemble number
  - `{YEAR 1}` and `{YEAR 2}` are the start and end year of the data in the file.
 


## Tests
[Pytests](https://docs.pytest.org/en/latest/index.html) were written to test and support the code. These can be run by:
```
pytest main_tests.py
```
Note: if you do not have pytest installed, then you will need to install it with:
```
pip install -U pytest
```

## Additional info
- The names and path of your `data` and `results` folder can be changed in the python file `directories.py`, if needed.
- Ensure that your files do **not** have spaces between them as commands used in **nco** will **not** work. 


## Built With

Python 3

## Contributing

Please contact a team member directly if you would like to be involved with the development of this software.

## Versioning

Currently released version 1.0.0 

## Author

**Adanna Akwataghibe**

## Acknowledgments 

- Thank you to Dr. Yves Plancherel, my supervisor for his support.
- Thank you to my family and friends for their moral support. 

## License

This project is licensed under the MIT License - see the [LICENSE.md](LICENSE.md) file for details
