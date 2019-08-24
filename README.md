# IRP - A software package for Climate Modelling diagnostics

The Climate Modelling (Diagnostics) Software is written using Python 3 and enables users to compute climate model diagnostics. This software makes use of many python modules and libraries. These libraries must be installed in order to successfully run the program.

It is **highly recommended** to go through the UserGuide.pdf as this gives detailed instructions on how to setup, run and test the program. 

Below is the *quick setup and run* of the software. 

## Getting Started (Quick Setup)
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
pip install -r INSTALL.txt
```

## Usage

### Input and output folders
The `.nc` data files to be analysed should be stored in the folder `DATA`. The `INPUT` folder stores all input, input masks, input sample points etc. files. The `RESULTS` folder will store the analysis, once computed. 
The file names can be changed in the file `directories.py`

There are two ways to call the program: 

### Using an input file

An input file is set up for the user to fill in the required values in ``INPUT/input.txt``. 
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
Note that the commands must be done within the `ClimateModelling` folder. 
#### Example input file
A pre-filled input file is given as an example in `INPUT/input_example.txt`, and can be run by calling:
```
python main.py -ex
```
This example input file calculates the mean of `air_temperature` given in `DATA\E1_north_america_ens101_1970.nc`.

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

## Tests
[Pytests](https://docs.pytest.org/en/latest/index.html) were written to test and support the code. The unit tests are stored in `ClimateModelling/tests` These can be run by:
```
pytest
```
Note: if you do not have pytest installed, then you will need to install it with:
```
pip install -U pytest
```

### Integration/Scenario tests 
In the folder `ClimateModelling/integration_tests`, a list of scenarios have been run and saved. 
You can run the program using the input and data provided to see if the results obtained match the expected results.

## Brief description of files
* `ClimateModelling`: contains classes and scripts that make up the software. 
    * 

* `ClimateModelling_Parallel`: The parallelised version of the code





## Additional info
- The names and path of your `data` and `results` folder can be changed in the python file `directories.py`, if needed.
- Ensure that your files do **not** have spaces between them as running on Linux or Anaconda Prompt will not work. 


## Built With

Python 3

## Contributing

Please contact the author if you would like to be involved with the development of this software.

## Versioning

Currently released version 1.0.0 

## Author

**Adanna Akwataghibe** (aa14415@ic.ac.uk)

## Acknowledgments 

- Thank you to Dr. Yves Plancherel, my supervisor for his support.
- Thank you to my family and friends for their moral support. 

## License

This project is licensed under the MIT License - see the [LICENSE.md](LICENSE.md) file for details
