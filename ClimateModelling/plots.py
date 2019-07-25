import matplotlib.pyplot as plt
from utils import *
import iris
import iris.plot as iplt
import iris.quickplot as qplt
from iris.time import PartialDateTime
import pandas as pd
import seaborn as sns
from iris.pandas import as_data_frame, as_series
from pandas.plotting import register_matplotlib_converters


def create_histogram(list_ens, start_date, end_date, variables, monthly=False, save_out=True, sel='fd', plot=True, ens_num=1):
    """
    Plot or save histogram data of cube
    :param list_ens: list of ensemble data
    :param start_date: [day, month, year] used for constructing title for plot and file name
    :param end_date: [day, month, year] used for constructing title for plot and file name
    :param variables: If one variable - then have 1D histogram. If list of 2 variables, then 2D histogram
    :param monthly: used for title for plot and file name
    :param save_out: if set, then save histogram frequency and values in text file
    :param sel: type of bin selection for 1D histogram
    :param plot: if set, then plots the histogram in a new window
    :param ens_num: selection of ensemble to use
    """

    # Make sure data structures are not empty
    assert list_ens is not None
    assert variables is not None

    # Daily or monthly
    time_str = "daily"
    if monthly:
        time_str = "monthly"

    # If variables is just one object, cast to list
    if not isinstance(variables, list):
        variables = [variables]
    elif len(variables) > 2:
        print("Error: function create_histogram: cannot construct histogram with more than 2-dimensions.")
        sys.exit()

    cov = len(variables) == 2
    # Get flattened data (without nan values) from cubes
    cubes = []
    full_datum = []
    for var in variables:
        # Get cube from dictionary and flatten data
        cube = list_ens[ens_num][var]
        full_data = cube.data.flatten()
        # Get indices of where value is nan
        indices = np.argwhere(np.isnan(full_data.data))
        # Remove nan values from full data
        full_data = np.delete(full_data, indices)
        cubes.append(cube)
        full_datum.append(full_data)

    # Construct title of histogram
    part_title = " measured " + time_str + " between " + str(start_date[2]) + "-" + str(start_date[1]) + "-" + \
            str(start_date[0]) + " and " + str(end_date[2]) + "-" + str(end_date[1]) + "-" + str(end_date[0]) \
            + " using the E2S2M climate model"
    if not cov:
        title_name = cubes[0].name() + part_title
    else:
        title_name = cubes[0].name() + " and " + cubes[1].name() + part_title

    # Plot histogram
    if plot:
        sns.set()
        fig, ax = plt.subplots()
        if not cov:
            ax.set_title(title_name)
            ax.hist(full_datum[0], bins=sel, color="skyblue")
            ax.set_ylabel("Frequency")
            ax.set_xlabel(cubes[0].name() + ' (' + str(cubes[0].units) + ')')

        else:  # Plot 2D histogram
            plt.title(title_name)
            plt.hist2d(full_datum[0], full_datum[1])
            plt.xlabel(cubes[0].name() + ' (' + str(cubes[0].units) + ')')
            plt.ylabel(cubes[1].name() + ' (' + str(cubes[1].units) + ')')
            plt.colorbar(label='Frequency')

    if save_out:
        # Construct file name
        title_name = " " + time_str + " " + str(start_date[2]) + "-" + str(start_date[1]) + "-" + \
            str(start_date[0]) + "_" + str(end_date[2]) + "-" + str(end_date[1]) + "-" + str(end_date[0])

        if not cov:
            title_name = "hist_" + cubes[0].name() + title_name
            title_name = make_into_file_name(title_name)
            file = open(os.path.join(directories.ANALYSIS, title_name), 'w')
            # Write header of file
            file.write("# " + cubes[0].name() + "," + "Frequency\n")
            # Get histogram data and write to file
            hist, bins = np.histogram(full_datum[0], bins=sel)
            data = list(zip(bins, hist))
            np.savetxt(file, data, delimiter=',')
            file.close()
        else:
            title_name = "hist2D_" + cubes[0].name() + "_" + cubes[1].name() + title_name
            title_name = make_into_file_name(title_name)
            file = open(os.path.join(directories.ANALYSIS, title_name), 'w')
            file.write("# " + cubes[0].name() + "," + cubes[1].name() + "," + "Frequency\n")
            hist, xs, ys = np.histogram2d(full_datum[0].data, full_datum[1].data)
            # Remove unnecessary last element
            xs = xs[:-1]
            ys = ys[:-1]

            # Get each corresponding value and write to file
            for i in range(len(xs)):
                for j in range(len(ys)):
                    file.write(str(xs[i]) + "," + str(ys[j]) + "," + str(hist[i][j]) + "\n")
            file.close()

        print("Histogram data is saved in the " + directories.ANALYSIS + " folder as a txt file.")


def create_timeseries(list_ens, start_date, end_date, variables,  monthly=False, save_out=True, plot=False, ens_num=1):
    """
    Analysis the data given - in this case it computes the timeseries (assumes grid/sample point)
    :param list_ens: the list of ensembles (dicts) containing the data of the climate variables
    :param units: the units matching to each variable
    :param start_date and end_date: extract data within this time frame
    :param variables: If one variable - then have 1D histogram. If list of 2 variables, then 2D histogram
    :param monthly: data is stored in monthly increments (time = 12) else assumed (time = 365)
    :param save_out: if set, then save output of histogram/ rimeseries
    :param cov: if set, then perform covariance analysis
    :return: None
    """

    # Make sure data structures are not empty
    assert list_ens is not None
    assert variables is not None

    # If variables is just one object, cast to list
    if not isinstance(variables, list):
        variables = [variables]

    # Get list of years to use for plots
    years = list(range(start_date[-1], end_date[-1]+1))

    # Daily or monthly
    time_str = "daily"
    if monthly:
        time_str = "monthly"

    # Construct figure and axes
    sns.set(rc={'figure.figsize': (11, 4)})
    # Convert pandas dates to matplotlib date format
    register_matplotlib_converters()

    for variable in variables:
        # Get cube from dictionary
        cube = list_ens[ens_num][variable]

        title_name = cube.name() + " measured " + time_str + " between " + str(start_date[2]) + "-" + str(start_date[1]) + "-" + \
                     str(start_date[0]) + " and " + str(end_date[2]) + "-" + str(end_date[1]) + "-" + str(end_date[0]) \
                     + " using the E2S2M climate model"

        # Convert cube to pandas series dataframe
        pd_data = None
        try:
            pd_data = as_series(cube)
        except Exception:
            print("Error: function create_timeseries: cannot construct timeseries with more than 1-dimension cube.")

        # Get x values in series
        indices = pd_data.index
        # Save dates that are in datetime format
        new_indices = []
        selected_indices = []
        for i in range(len(indices)):
            # Convert cftime.datetime to datetime
            dt = convert_cftime_datetime(indices[i])
            new_indices.append(dt)
            # Get first of the year
            if dt.month == 1 and dt.day == 1:
                selected_indices.append(indices[i])

        indx = pd.DatetimeIndex(new_indices)
        # df = pd.Series(pd_data.to_numpy(), index=indx)
        df = pd.DataFrame(data=pd_data.to_numpy(), index=indx, columns=[variable])

        # Save table to file for each variable
        # Construct file name
        date_str = " " + time_str + " " + str(start_date[2]) + "-" + str(start_date[1]) + "-" + \
                     str(start_date[0]) + "_" + str(end_date[2]) + "-" + str(end_date[1]) + "-" + str(end_date[0])

        file_name = "ts_" + cube.name() + date_str
        file_name = make_into_file_name(file_name)
        if save_out:
            # np.savetxt(file_name, df.values, fmt='%d', comments="dates  " + cube.name())
            df.to_csv(os.path.join(directories.ANALYSIS, file_name), sep=',', index=True, index_label='dates')
            print("Timeseries data is saved in the " + directories.ANALYSIS + " folder as a txt file.")

        # Plot the timeseries
        fig, axs = plt.subplots(2, 1)
        plt.tight_layout()
        # Timeseries
        df.plot(ax=axs[0], title=title_name, grid=True, legend=False, alpha=0.7, color='m', ls='-')
        df.resample('BM').mean().plot(ax=axs[0], style='--')
        rolling_str = 'one-year rolling mean'
        if monthly:
            df.rolling(12, center=True).mean().plot(ax=axs[0], style=':')
        else:
            df.rolling(365, center=True).mean().plot(ax=axs[0], style=':')
        axs[0].legend(['input', 'monthly mean', rolling_str], loc='upper left')
        axs[0].set_xlabel("Time (month/year)")
        axs[0].set_ylabel(cube.name())

        # Boxplot - Yearly seasonality
        df['Month'] = df.index.month
        sns.boxplot(ax=axs[1], data=df, x='Month', y=variable)
        axs[1].set_title("Yearly seasonality of " + cube.name())
        axs[1].set_xlabel("Months")
        axs[1].set_ylabel(cube.name())

        # Boxplot - Weekly seasonality
        # df['Weekday name'] = df.index.weekday_name
        # sns.boxplot(ax=axs[2], data=df, x='Weekday name', y=variable)
        # axs[2].set_title("Weekly seasonality of " + cube.name())
        # axs[2].set_xlabel("Days")
        # axs[2].set_ylabel(cube.name())


def plot_graph(list_ens, current_date, variables, mask=None, ens_num=1, save_out=False, lat=None, lon=None):
    """
    Plot (and save) graphs given date and variable in cube
    :param list_ens: the list of ensembles (dicts) containing the data of the climate variables
    :param current_date: list of [day, month, year]
    :param variables: variables to select in cube
    :param mask: if given, plot with mask applied
    :param ens_num: gives which ensemble to plot, default = 1
    :param save_out: if set, then saves to png file
    :param lat: latitude, float, if given, mark on graph
    :param lon: longitude, float, if given, mark on graph
    """
    # Make sure data structures are not empty
    assert list_ens is not None
    assert variables is not None

    # If variables is just one object, cast to list
    if not isinstance(variables, list):
        variables = [variables]

    # Unpack current date
    cur_day, cur_month, cur_yr = current_date[0], current_date[1], current_date[2]
    # Construct date
    pdt = PartialDateTime(year=cur_yr, month=cur_month, day=cur_day)
    date_range = iris.Constraint(time=pdt)

    for variable in variables:
        # Extract time in cube
        cube = list_ens[ens_num][variable].extract(date_range)

        # Plot title name
        title_name = cube.name() + " at " + str(cur_yr) + "-" + str(cur_month) + "-" + str(
            cur_day) + " with variable " + str(variable)

        sns.set()
        # Plot the graph
        fig, ax = plt.subplots()
        # Get the Purples "Brewer" palette for plot
        brewer_cmap = plt.get_cmap('brewer_Purples_09')

        # Draw the contours, with n-levels set for the map colours
        lats = cube.coord('latitude')
        lons = cube.coord('longitude')
        qplt.contourf(cube, cmap=brewer_cmap)
        # cube.intersection(longitude=(-45, 45))  # get specific range
        plt.xlabel("Longitude (" + str(lons.units) + ")")
        plt.ylabel("Latitude (" + str(lats.units) + ")")

        # x and y axis
        x_btm, x_top = plt.xlim()
        y_btm, y_top = plt.ylim()

        plt.xticks(np.arange(x_btm,  x_top, step=20))
        plt.yticks(np.arange(y_btm,  y_top, step=30))

        # Add a citation to the plot.
        iplt.citation(iris.plot.BREWER_CITE)

        plt.scatter(lon, lat, s=100, c='black', marker='x', label="Point (" + str(lon) + ", " + str(lat) + ")")

        # Overlay mask
        xs, ys = mask
        if xs is not None:
            for i in range(len(xs)):
                plt.fill(xs[i], ys[i], alpha=0.7, label='Mask')
            title_name = "Applied masks to " + cube.name() + " with variable " + str(variable)
        plt.legend(frameon=True)
        plt.title(title_name)

        if save_out:
            title_name = make_into_file_name(title_name)
            plt.savefig(os.path.join(directories.ANALYSIS, title_name))
            print("Map is saved in the " + directories.ANALYSIS + " folder as a png file.")
