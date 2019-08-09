import matplotlib.pyplot as plt
from utils import *
import iris
import iris.plot as iplt
import pandas as pd
import seaborn as sns
from iris.pandas import as_data_frame, as_series
from pandas.plotting import register_matplotlib_converters
import sys


def create_histogram(list_ens, start_date, end_date, variables, monthly=False, save_out=True, sel='fd', ens_num=1,
                     cov=False, mask=None):
    """
    Plot or save histogram data of cube
    :param list_ens: list of ensemble data
    :param start_date: [day, month, year] used for constructing title for plot and file name
    :param end_date: [day, month, year] used for constructing title for plot and file name
    :param variables: If one variable - then have 1D histogram. If list of 2 variables, then 2D histogram
    :param monthly: used for title for plot and file name
    :param save_out: if set, then save histogram frequency and values in text file
    :param sel: type of bin selection for 1D histogram
    :param ens_num: selection of ensemble to use
    :param cov: if covariance between 2 values, then param is set to true
    :param mask: set to true if masked data given
    """

    # Make sure data structures are not empty
    assert list_ens is not None
    assert variables is not None

    # If file is given, select the bin edges in the file
    bins_dict = None
    if '.' in sel:
        bins_dict = get_bins_from_file(sel, variables)

    # Daily or monthly
    time_str = "daily"
    if monthly:
        time_str = "monthly"

    # If variables is just one object, cast to list
    if not isinstance(variables, list):
        variables = [variables]
    elif len(variables) > 2:
        print("ERROR: function create_histogram: cannot construct histogram with more than 2-dimensions.")
        sys.exit()

    # Get flattened data (without nan values) from cubes
    cubes = []
    full_datum = []
    for var in variables:
        # Get cube from dictionary and flatten data
        cube = list_ens[ens_num][var]
        full_data = cube.data.flatten()
        # Get indices of where value is nan
        indices = []
        if mask is None or not mask:
            indices = np.argwhere(np.isnan(full_data.data))
        else:
            x = cube.data.filled()
            indices = np.argwhere(np.isclose(x.flatten().data, cube.data.fill_value))
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
    sns.set()
    for i in range(len(variables)):
        if not cov:
            fig, ax = plt.subplots()
            ax.set_title(title_name)
            if bins_dict is None:
                ax.hist(full_datum[i], bins=sel, color="skyblue")
            else:
                ax.hist(full_datum[i], bins=bins_dict[variables[i]], color="skyblue")
            ax.set_ylabel("Frequency")
            ax.set_xlabel(cubes[i].name() + ' (' + str(cubes[i].units) + ')')

    if cov:  # Plot 2D histogram
        fig, ax = plt.subplots()
        ax.set_title(title_name)
        _, _, _, img = ax.hist2d(full_datum[0], full_datum[1])
        ax.set_xlabel(cubes[0].name() + ' (' + str(cubes[0].units) + ')')
        ax.set_ylabel(cubes[1].name() + ' (' + str(cubes[1].units) + ')')
        plt.colorbar(img, ax=ax, label='Frequency')

    if save_out:
        # Construct file name
        title_name = " " + time_str + " " + str(start_date[2]) + "-" + str(start_date[1]) + "-" + \
            str(start_date[0]) + "_" + str(end_date[2]) + "-" + str(end_date[1]) + "-" + str(end_date[0])

        for i in range(len(variables)):
            if not cov:
                title_name = "hist_" + cubes[i].name() + title_name
                title_name = make_into_file_name(title_name)
                file = open(os.path.join(directories.ANALYSIS, title_name), 'w')
                # Write header of file
                file.write("# " + cubes[i].name() + "," + "Frequency\n")
                # Get histogram data and write to file
                if bins_dict is None:
                    hist, bins = np.histogram(full_datum[i], bins=sel)
                else:
                    hist, bins = np.histogram(full_datum[i], bins=bins_dict[variables[i]])
                data = list(zip(bins, hist))
                np.savetxt(file, data, delimiter=',')
                file.close()
                plt.savefig(os.path.join(directories.ANALYSIS, title_name))
        if cov:
            title_name = "hist2D_" + cubes[0].name() + "_" + cubes[1].name() + title_name
            title_name = make_into_file_name(title_name)
            file = open(os.path.join(directories.ANALYSIS, title_name), 'w')
            file.write("# " + cubes[0].name() + "," + cubes[1].name() + "," + "Frequency\n")
            if bins_dict is None:
                hist, xs, ys = np.histogram2d(full_datum[0].data, full_datum[1].data)
            else:
                hist, xs, ys = np.histogram2d(full_datum[0].data, full_datum[1].data,
                                          bins=[bins_dict[variables[0]], bins_dict[variables[1]]])
            # Remove unnecessary last element
            xs = xs[:-1]
            ys = ys[:-1]

            # Get each corresponding value and write to file
            for i in range(len(xs)):
                for j in range(len(ys)):
                    file.write(str(xs[i]) + "," + str(ys[j]) + "," + str(hist[i][j]) + "\n")
            file.close()

            plt.savefig(os.path.join(directories.ANALYSIS, title_name))

        print("Histogram data is saved in the " + directories.ANALYSIS + " folder as a txt file.")


def create_timeseries(list_ens, start_date, end_date, variables,  monthly=False, save_out=True, ens_num=1,
                      func_name=None):
    """
    Analysis the data given - in this case it computes the timeseries (assumes grid/sample point)
    :param list_ens: the list of ensembles (dicts) containing the data of the climate variables
    :param start_date extract from end date from data, list [d, m, y]
    :param end_date: extract till end date from data, list [d, m, y]
    :param variables: If one variable - then have 1D histogram. If list of 2 variables, then 2D histogram
    :param monthly: data is stored in monthly increments (time = 12) else assumed (time = 365)
    :param save_out: if set, then save output of histogram/ rimeseries
    :param ens_num: selection of ensemble to use
    :param func_name: if user function analysis used, this is the function name
    """

    # Make sure data structures are not empty
    assert list_ens is not None
    assert variables is not None

    # If variables is just one object, cast to list
    if not isinstance(variables, list):
        variables = [variables]

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
        # Construct title name
        title_name = cube.name() + " measured " + time_str + " between " + str(start_date[2]) + "-" + str(start_date[1]) + "-" + \
                     str(start_date[0]) + " and " + str(end_date[2]) + "-" + str(end_date[1]) + "-" + str(end_date[0]) \
                     + " using the E2S2M climate model"
        if func_name is not None:
            title_name = func_name + " of " + title_name
        # Convert cube to pandas series dataframe
        pd_data = None
        try:
            pd_data = as_series(cube)
        except Exception:
            print("WARNING: function create_timeseries: cannot construct timeseries with more than 1-dimension cube.")
            return None

        # Get x values in series
        indices = pd_data.index
        if len(indices) == 1:
            print("WARNING in function create_timeseries: cannot construct timeseries of one value.")
            return None
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
        axs[1].set_title("Yearly seasonality of " +cube.name())
        axs[1].set_xlabel("Months")
        axs[1].set_ylabel(cube.name())

        # Boxplot - Weekly seasonality
        # df['Weekday name'] = df.index.weekday_name
        # sns.boxplot(ax=axs[2], data=df, x='Weekday name', y=variable)
        # axs[2].set_title("Weekly seasonality of " + cube.name())
        # axs[2].set_xlabel("Days")
        # axs[2].set_ylabel(cube.name())


def plot_map(list_ens, variables, analysis_str=None, ens_num=1, save_out=False, lat=None, lon=None, total=False):
    """
    Plot (and save) maps given date and variable in cube
    :param list_ens: the list of ensembles (dicts) containing the data of the climate variables
    :param variables: variables to select in cube
    :param analysis_str: list of analysis names performed on ensemble, list of strings
    :param ens_num: gives which ensemble to plot, default = 1
    :param save_out: if set, then saves to png file
    :param lat: latitude, float, if given, mark on graph
    :param lon: longitude, float, if given, mark on graph
    :param total: set if all ensembles have been calculated together, instead of separately, boolean
    """
    # Make sure data structures are not empty
    assert list_ens is not None
    assert variables is not None
    assert analysis_str is not None

    # If variables is just one object, cast to list
    if not isinstance(variables, list):
        variables = [variables]

    for variable in variables:
        for a in range(len(analysis_str)):
            # Extract time in cube
            cube = None
            if len(analysis_str) > 1:
                if total:
                    cube = list_ens[a][variable]
                else:
                    cube = list_ens[ens_num][a][variable]
            else:
                if total:
                    cube = list_ens[variable]
                else:
                    cube = list_ens[ens_num][variable]

            # Plot title name
            title_name = analysis_str[a] + " of " + cube.name() + " of ensemble " + str(ens_num) + " with variable " + str(variable)
            if total:
                title_name = analysis_str[a] + " of " + cube.name() + " of all ensembles" + " with variable " + str(variable)

            # Plot the graph
            fig, ax = plt.subplots(figsize=(11, 4))
            # Get the Purples "Brewer" palette for plot
            brewer_cmap = plt.get_cmap('brewer_Purples_09')

            # Draw the contours, with n-levels set for the map colours
            lats = cube.coord('latitude')
            lons = cube.coord('longitude')
            masked_array = np.ma.array(cube.data, mask=np.isnan(cube.data))

            try:
                cf = ax.contourf(lons.points, lats.points, masked_array, cmap=brewer_cmap)
            except Exception:
                print("ERROR in function plot_map: cube is not 2D.")
                return None
            # cube.intersection(longitude=(-45, 45))  # get specific range
            fig.colorbar(cf, ax=ax, label=cube.units)
            ax.set_xlabel("Longitude (" + str(lons.units) + ")")
            ax.set_ylabel("Latitude (" + str(lats.units) + ")")

            # Add a citation to the plot.
            iplt.citation(iris.plot.BREWER_CITE)

            # Plot X on map if point is given
            if lon is not None:
                ax.scatter(lon, lat, s=100, c='black', marker='x', label="Point (" + str(lon) + ", " + str(lat) + ")")

            # # Overlay mask
            # xs, ys = mask
            # if xs is not None:
            #     for i in range(len(xs)):
            #         ax.fill([-112, -112, -75, 45, 45], ys[i], alpha=0.7, label='Mask')
            #     title_name = "Applied masks to " + cube.name() + " with variable " + str(variable)

            ax.set_title(title_name)

            # Have legend if mask or point
            if lon is not None:
                plt.legend(frameon=True)

            if save_out:
                title_name = make_into_file_name(title_name)
                plt.savefig(os.path.join(directories.ANALYSIS, title_name))
                print("Map is saved in the " + directories.ANALYSIS + " folder as a png file.")