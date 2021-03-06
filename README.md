# **EZclim** - *A software package for Climate Modelling diagnostics*

(**IRP Report** : in file [**IRP_REPORT.pdf**](IRP_REPORT.pdf))

The Climate Modelling (Diagnostics) Software is written using Python 3 and enables users to compute climate model diagnostics. This software makes use of many python modules and libraries. These libraries must be installed in order to successfully run the program.

It is **highly recommended** to go through the [**User Guide**](UserGuide.pdf) as this gives detailed instructions on how to setup, run and test the program. 

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
    cd acse-9-independent-research-project-AdannaAkwats

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
# User input for Climate Modelling Diagnostics Program: EZclim
#
# ----------------- PLEASE FILL IN THE ARGUMENTS LISTED BELOW ------------------
#
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
Analysis:
Spatial:
Total ensemble stats:
Plot: 1
Monthly:
Grid:
Sample:
Mask file:
Save Output: True
Covary:
Histogram bin selection:
Longitude centre:
User function:
Calculate areas:
Calculate index:
#
# ------------------------------------------------------------------------------
# HELP : Found in the UserGuide
# ------------------------------------------------------------------------------
```
All ommands must be done within the `EZclim` folder, so:
```
cd EZclim
```
After filling the values in the file, to run the program simply call: 
```
python main.py
```
#### Example input file
A pre-filled input file is given as an example in `INPUT/input_example.txt`, and can be run by calling:
```
python main.py -ex
```
This example input file calculates the mean of `air_temperature` given in `DATA\E1_north_america_ens101_1970.nc`.

### Command line interface (CLI)
```
 python main.py -h 
 less output.log
```
shows the list of commands available (saved in `output.log`) to use and their descriptions. 
#### Example commands
The file `E1_north_america_ens101_1970.nc` is given as example data in the `DATA` folder. To calculate the mean of variable `air_temperature`, and plot the histogram and map, we would call:
```
python main.py 1970 -pf E1_north_america -v air_temperature -e 1 -a mean -p 1
```

## Debugging
When the program is run successfully, the following message prints out:
```
python main.py 
...

PROGRAM FINSISHED.
>> Progress and potential errors are logged in output.log file.
>> To open output.log, type in cmd line: less output.log
```
If no message prints out, it is because an error has been seen in the code, this may be caused by invalid inputs being seen. 
Open the `output.log` file created and this will give you a summary of what the code has been doing, and any errors that may have been thrown. 

## Tests
[Pytests](https://docs.pytest.org/en/latest/index.html) were written to test and support the code. The unit tests are stored in `EZclim/tests` These can be run by:
```
pytest
```
Note: if you do not have pytest installed, then you will need to install it with:
```
pip install -U pytest
```

### Integration/Scenario tests 
In the folder `EZclim/integration_tests`, a list of scenarios have been run and saved. 
You can run the program using the input and data provided to see if the results obtained match the expected results.

## Brief description of files
* `EZclim`: contains classes and scripts that make up the software. 
    * `DATA/` : stores NetCDF model output files 
    * `INPUT/` : stores input.txt, input_example.txt and masks, histogram bins, sample points files. 
    * `RESULTS/` : stores the output of run of program
        * `ensemble_averages/`: stores the NetCDF output of program
    * `tests/` : stores unit tests written with pytest
    * `integration_tests/` : stores tests cases used to check all components are working. 
    * `Analysis.py` : class that performs analysis on data 
    * `Extract.py` : class that extracts useful data from NetCDF files given by user
    * `Months.py` : class that has months of the year (needed by user_entry.py)
    * `ProgressBar` : class that constructs the progress bar used when program is run.
    * `WriteOutput`: class that write final data to NetCDF files. 
    * `calculate_indices.py` : python wrapper that calls bash script calculate_indices.sh given data 
    * `calculate_indices.sh` : bash script that calculates indices (enso, tno, pdo ...)
    * `create_test_data_files.py` : used to create dummy NetCDF files to test with
    * `directories.py` : stores the mappings of folders 
    * `file_entry.py` : script that deals with extracting data from input.txt and input_example.txt
    * `main.py` : the entry point of the program
    * `plots.py` : script that contains functions that plots histograms, time series, boxplots and maps
    * `user_entry` : script that parses command line arguments and calls all functions to extract, analyse, plot and write output. 
    * `utils.py` : script that contains useful functions used across the code base. 
    * (`output.log`) : when the program is run, all progress/error messages is stored here.  

* `EZclim_Parallel`: The parallelised version of the code

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
