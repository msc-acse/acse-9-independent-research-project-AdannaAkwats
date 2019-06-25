# IRP - A software package for Climate Modelling diagnostics

Written in Python 3 and contains functions that enable generate relevant diagnostics. 

Currently when called, the program calculates the mean of each ensemble and saves the data in folder ``analysis\ensemble_means``.
## Getting Started
- Download / clone the repository unto your computer
- Install **netCDF4** ([pip](https://pypi.org/project/netCDF4/) or [conda](https://anaconda.org/anaconda/netcdf4))
    ##### Using PyPI
    ```
    pip install netCDF4
    ```
    ##### Using conda
    ```
    conda install -c anaconda netcdf4 
    ```
- Install **nco** ([pip](https://pypi.org/project/nco/) or [conda](https://anaconda.org/conda-forge/nco))
    ##### Using PyPI
    ```
    pip install nco
    ```
    ##### Using conda
    ```
    conda install -c conda-forge nco
    ```

- Install the required python package dependencies with the following command:
```
pip install -r requirements.txt
```

## Usage

### Input and output folders
The `.nc` data files to be analysed should be stored in the folder `data`. The `results` folder will store the analysis, once computed. 

There are two ways to call the program: 

### Command line interface (CLI)
```
python main.py -h 

usage: CLIMATE_ANALYSIS [-h] -v variables [variables ...] [-p] [-m]
                        [-g lat lon | -s lat lon] [-mk filename] [-o] [-cv] -e
                        number_of_ensembles
                        [-ht [number_of_bins_in_histogram]]
                        prefix start_date [end_date]

The functions will give statistical analysis of the climate data
                                     presented
    FILENAMES FORMAT
    ----------------
    - The filenames should be in the format "{START OF FILENAME}_ens{NUM}_{YEAR}.nc", where {START OF FILENAME} is
    the prefix of the file, this can be the algae type etc, {NUM} is the ensemble number and {YEAR} is the year.
   OR if you have multiple years stored in one file then:
   - The filenames should be in the format "{START OF FILENAME}_ens{NUM}_{YEAR 1}_{YEAR 2}.nc", where
   {START OF FILENAME} is the prefix of the file, this can be the algae type etc, {NUM} is the ensemble number and
   {YEAR 1} and {YEAR 2} are the start and end year of the data in the file.
   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    ASSUMPTIONS
    ------------
    - Files do not have overlapped data.
    - Daily increments of data, except if the monthly tag is set in the arguments.
    - Grids have constant latitude and longitude.
    ------------
    - Some example files are in the data folder.


positional arguments:
  prefix                This is the prefix of the file. This can be the type of algae e.g. dic_deltap, fndet_100, jprod_ndi_100.
  start_date            Start date of analysis
                            Can be in the following formats:
                            ----------------------------------
                            YYYY-MM-DD : e.g. 2020-04-12
                            YYYY-MM    : e.g. 2020-04
                            YYYY       : e.g. 2020
                            - If day is not given, the 1st of the given month will be used i.e 2020-04 => 2020-04-01
                            - If day and month is not given, 1st Jan will be used as the start date i.e 2020 => 2020-01-01
  end_date               <Not required> End date of analysis - format is the same as start_date
                            -----------------------------------end_date not given-------------------------------------
                            - If only start year is given, the end_date is automatically set to the 31 Dec of start year
                            - If start year and month is given, then end_date is set to the end of the start month
                               -----------------------------------end_date given-------------------------------------
                            - If day is not given, the end of the given month will be used i.e 2020-04 => 2020-04-30
                            - If day and month is not given, 31 Dec will be used as the end date i.e 2020 => 2020-12-31

other arguments:
  -h, --help            show this help message and exit
  -v variables [variables ...], --vars variables [variables ...]
                        <Required> Variables of data to analyse
  -p, --plot            Save plots of analysis in results as a .png file.
  -m, --monthly         Data in file is stored in monthly increments.
  -g lat lon, --grid lat lon
                        Uses grid point that latitude and longitude lies in.
  -s lat lon, --sample lat lon
                        Uses sample point given by latitude and longitude using interpolation.
  -mk filename, --mask filename
                        Uses masking grid given as a file (contains boolean array to be imposed on the global grid).
  -o, --output          Save data output of histogram and timeseries analysis in results as a .dat file.
  -cv, --covary         Analysis on how the variables given in -v vary with each other.
  -e number_of_ensembles, --ens number_of_ensembles
                        <Required> The number of ensembles of the data. If not set, the default value = 1
  -ht [number_of_bins_in_histogram], --hist [number_of_bins_in_histogram]
                         Options for bin size selection. If not set, the default value = fd (Freedman Diaconis Estimator). The list of the potential options are listed in:
                        https://docs.scipy.org/doc/numpy/reference/generated/numpy.histogram_bin_edges.html#numpy.histogram_bin_edges
```
#### Example commands
```
python main.py sst 1953 1955 -v sst -e 12
```
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
Plot: False
Monthly: False
Grid:
Sample:
Mask file:
Save Output: False
Covary: False
Histogram bin selection: fd
#
```
After filling the values in the file, to run the program simply call: 
```
python main.py
```
An pre-filled input file is given as an example in `input_example.txt`.

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
