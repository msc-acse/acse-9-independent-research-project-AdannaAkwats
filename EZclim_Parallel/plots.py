"""
Created by Adanna Akwataghibe (Github: AdannaAkwats)
"""
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap
from utils import *
import iris
import iris.plot as iplt
import pandas as pd
import seaborn as sns
from iris.pandas import as_data_frame, as_series
from pandas.plotting import register_matplotlib_converters
import directories
from collections import defaultdict


def create_histogram(list_ens, ens_stats, start_date, end_date, variables, monthly=False, save_out=True, sel='fd', ens_num=1,
                     cov=False, mask=None, total=False, analysis_str=None, nan_indices=None, plot=None,
                     second_date_given=False, start_date2=None, end_date2=None, spatial=False):
    """
    Plot or save histogram data of cube
    :param list_ens: list of ensemble data
    :param ens_stats: list of analysed ensemble data
    :param start_date: [day, month, year] used for constructing title for plot and file name
    :param end_date: [day, month, year] used for constructing title for plot and file name
    :param variables: If one variable - then have 1D histogram. If list of 2 variables, then 2D histogram
    :param monthly: used for title for plot and file name
    :param save_out: if set, then save histogram frequency and values in text file
    :param sel: type of bin selection for 1D histogram
    :param ens_num: selection of ensemble to use
    :param cov: if covariance between 2 values, then param is set to true
    :param mask: set to true if masked data given
    :param total: set if all ensemble are averaged into one
    :param analysis_str: list of analysis performed on models
    :param nan_indices: indices where values are nan in data (for temporal average)
    :param plot: setting to show plot in a figure popup
    :param second_date_given: if set, multi model averages calculated
    :param start_date2: second model start date
    :param end_date2: second model end date
    :param spatial: spatial averages have been calculated
    """

    # Make sure data structures are not empty
    assert list_ens is not None
    assert variables is not None

    # If file is given, select the bin edges in the file
    bins_dict, x, y = None, None, None
    if '.' in sel[0]:
        hist_file = os.path.join(directories.INPUT, sel[0])
        bins_dict, x, y = get_bins_from_file(hist_file, variables)
    elif len(sel) == 1:
        try:
            sel = int(sel[0])
        except Exception:
            sel = sel[0]
    elif len(sel) == 2:
        try:
            sel[0] = int(sel[0])
            sel[1] = int(sel[1])
        except Exception:
            print("ERROR in function create_histogram: bins size not valid: " + str(sel) + " - both elements should be integers.")
            return None

    # Daily or monthly
    time_str = "daily"
    if monthly:
        time_str = "monthly"

    # If variables is just one object, cast to list
    if not isinstance(variables, list):
        variables = [variables]
    elif len(variables) > 2:
        print("ERROR: function create_histogram: cannot construct histogram with more than 2-dimensions.")
        return None

    # Plot histogram of analysis performed
    if ens_stats is not None and analysis_str:
        create_analysis_histogram(ens_stats, start_date, end_date, variables, time_str=time_str, save_out=save_out, sel=sel,
                         ens_num=ens_num, cov=cov, mask=mask, total=total, analysis_str=analysis_str, bins_dict=bins_dict,
                                  x=x, y=y, nan_indices=nan_indices, plot=plot, second_date_given=second_date_given,
                                  start_date2=start_date2, end_date2=end_date2, spatial=spatial)

    # If total and not analysis
    if total and not analysis_str:
        list_ens, ens_num = ens_stats, 1

    if not analysis_str:
        # Get flattened data (without nan values) from cubes
        cubes = []
        full_datum = []
        # If cov and bin dict set
        var1, var2 = 0, 1
        if cov and bins_dict is not None:
            if x == variables[1] and y == variables[0]:
                var1, var2 = 1, 0

        for var in variables:
            # Get cube from dictionary and flatten data
            cube = list_ens[ens_num-1][var]
            full_data = np.asarray(cube.data.flatten())
            # Get indices of where value is nan
            indices = []
            if mask is None or not mask:
                indices = np.argwhere(np.isnan(full_data))
            else:
                x = cube.data.filled()
                indices = np.argwhere(np.isclose(x.flatten().data, cube.data.fill_value))
            # Remove nan values from full data
            full_data = np.delete(full_data, indices)
            cubes.append(cube)
            full_datum.append(full_data)

        # Construct title of histogram
        part_title = " measured " + time_str + " between " + str(start_date[2]) + "-" + str(start_date[1]) + "-" + \
                str(start_date[0]) + " and " + str(end_date[2]) + "-" + str(end_date[1]) + "-" + str(end_date[0])

        # Construct file name
        part_name = " " + time_str + " " + str(start_date[2]) + "_" + str(end_date[2])

        if plot is not None:
            # Plot histogram
            sns.set(rc={'figure.figsize': (11, 11)})
            for i in range(len(variables)):
                if not cov:
                    fig, ax = plt.subplots()
                    title_name = cubes[i].name() + part_title
                    ax.set_title(title_name)
                    if bins_dict is None:
                        if len(full_datum[i]) == 0:
                            print("ERROR in function create_histogram: No data obtained for histogram.")
                            return None
                        min_bin, max_bin = min(full_datum[i].data), max(full_datum[i].data)
                        try:
                            ax.hist(full_datum[i], bins=sel, range=(min_bin, max_bin), color="purple")
                        except (MemoryError, ValueError):
                            ax.hist(full_datum[i], bins=10, range=(min_bin, max_bin), color="purple")
                    else:
                        ax.hist(full_datum[i], bins=bins_dict[variables[i]], color="purple")
                    ax.set_ylabel("Frequency")
                    ax.set_xlabel(cubes[i].name() + ' (' + str(cubes[i].units) + ')')

                    f = "hist_" + cubes[i].name() + part_name
                    f = make_into_file_name(f)
                    plt.savefig(os.path.join(directories.ANALYSIS, f))

            if cov:  # Plot 2D histogram
                if len(full_datum[0]) == 0 or len(full_datum[1]) == 0:
                    print("ERROR in function create_histogram: No data obtained for histogram.")
                    return None
                fig, ax = plt.subplots()
                title_name = cubes[0].name() + " and " + cubes[1].name() + part_title
                ax.set_title(title_name)
                if bins_dict is None:
                    _, _, _, img = ax.hist2d(full_datum[0], full_datum[1], bins=sel)
                else:
                    _, _, _, img = ax.hist2d(full_datum[var1], full_datum[var2], bins=[bins_dict[x], bins_dict[y]])
                ax.set_xlabel(cubes[var1].name() + ' (' + str(cubes[var1].units) + ')')
                ax.set_ylabel(cubes[var2].name() + ' (' + str(cubes[var1].units) + ')')
                plt.colorbar(img, ax=ax, label='Frequency')

                f = "hist2D_" + cubes[0].name() + "_" + cubes[1].name() + part_name
                f = make_into_file_name(f)
                plt.savefig(os.path.join(directories.ANALYSIS, f))

        if save_out:
            for i in range(len(variables)):
                if not cov:
                    title_name = "hist_" + cubes[i].name() + part_name
                    title_name = make_into_file_name(title_name) + '.txt'
                    file = open(os.path.join(directories.ANALYSIS, title_name), 'w')
                    # Write header of file
                    file.write("# " + cubes[i].name() + "," + "Frequency\n")
                    # Get histogram data and write to file
                    if bins_dict is None:
                        if len(full_datum[i]) == 0:
                            print("ERROR in function create_histogram: No data obtained for histogram.")
                            return None
                        min_bin, max_bin = min(full_datum[i].data), max(full_datum[i].data)
                        if min_bin == max_bin:
                            print("ERROR in function create_histogram: min and max bin in histogram are equal.")
                            return None
                        hist, bins = None, None
                        try:
                            hist, bins = np.histogram(full_datum[i], range=(min_bin, max_bin), bins=sel)
                        except (MemoryError, ValueError):
                            hist, bins = np.histogram(full_datum[i], range=(min_bin, max_bin), bins=10)
                    else:
                        hist, bins = np.histogram(full_datum[i], bins=bins_dict[variables[i]])
                    data = list(zip(bins, hist))
                    np.savetxt(file, data, delimiter=',')
                    file.close()
            if cov:
                if len(full_datum[0]) == 0 or len(full_datum[1]) == 0:
                    print("ERROR in function create_histogram: No data obtained for histogram.")
                    return None
                title_name = "hist2D_" + cubes[0].name() + "_" + cubes[1].name() + part_name
                title_name = make_into_file_name(title_name) + '.txt'
                file = open(os.path.join(directories.ANALYSIS, title_name), 'w')
                file.write("# " + cubes[var1].name() + "," + cubes[var2].name() + "," + "Frequency\n")
                if bins_dict is None:
                    hist, xs, ys = np.histogram2d(full_datum[0], full_datum[1], bins=sel)
                else:
                    hist, xs, ys = np.histogram2d(full_datum[var1], full_datum[var2],
                                                  bins=[bins_dict[x], bins_dict[y]])

                # Get each corresponding value and write to file
                for i in range(len(xs)):
                    for j in range(len(ys)):
                        file.write(str(xs[i]) + "," + str(ys[j]) + "," + str(hist[i][j]) + "\n")
                file.close()

            print("Histogram data is saved in the " + directories.ANALYSIS + " folder as a txt file.")


def create_analysis_histogram(list_ens, start_date, end_date, variables, time_str=False, save_out=True, sel='fd',
                     ens_num=1, cov=False, mask=None, total=False, analysis_str=None, bins_dict=None,
                                  x=None, y=None, nan_indices=None, plot=None, second_date_given=False,
                              start_date2=None, end_date2=None, spatial=False):
    """
    Plot or save histogram data of cube
    :param list_ens: list of analysed ensemble data
    :param start_date: [day, month, year] used for constructing title for plot and file name
    :param end_date: [day, month, year] used for constructing title for plot and file name
    :param variables: If one variable - then have 1D histogram. If list of 2 variables, then 2D histogram
    :param time_str: used for title for plot and file name, either monthly or daily
    :param save_out: if set, then save histogram frequency and values in text file
    :param sel: type of bin selection for 1D histogram
    :param ens_num: selection of ensemble to use
    :param cov: if covariance between 2 values, then param is set to true
    :param mask: set to true if masked data given
    :param total: set if all ensemble are averaged into one
    :param analysis_str: list of analysis performed on models
    :param bins_dict: from histogrsm_bins file
    :param x: from histogram_bins file
    :param y: from histogram_bins file
    :param nan_indices: indices where values are nan in data (for temporal average)
    :param plot: setting to show plot in a figure popup
    :param second_date_given: if set, multi model averages calculated
    :param start_date2: second model start date
    :param end_date2: second model end date
    :param spatial: spatial averages have been calculated
    """

    # Make sure data structures are not empty
    assert list_ens is not None
    assert variables is not None

    # Get flattened data (without nan values) from cubes
    cubes, full_datum = defaultdict(list), defaultdict(list)
    # If cov and bin dict set
    var1, var2 = 0, 1
    if cov and bins_dict is not None:
        if x == variables[1] and y == variables[0]:
            var1, var2 = 1, 0

    for var in variables:
        for a in range(len(analysis_str)):
            # Get cube from dictionary and flatten data
            cube = None
            if total:
                cube = list_ens[a][var]
            else:
                cube = list_ens[ens_num - 1][a][var]
            full_data = np.asarray(cube.data.flatten())
            # Remove nan values from full data
            if not spatial:
                full_data = np.delete(full_data, nan_indices[var])
            full_data[full_data >= 0.5e+20] = np.nan
            full_data[full_data <= -1e+20] = np.nan
            full_data = full_data[~np.isnan(full_data)]
            cubes[a].append(cube)
            full_datum[a].append(full_data)

    # Construct title of histogram
    part_title = " measured " + time_str + " between " + str(start_date[2]) + "-" + str(start_date[1]) + "-" + \
                 str(start_date[0]) + " and " + str(end_date[2]) + "-" + str(end_date[1]) + "-" + str(
        end_date[0])
    # Construct file name
    part_name = " " + time_str + " " + str(start_date[2]) + "_" + str(end_date[2])
    if second_date_given:
        part_title = " measured " + time_str + " between " + str(start_date[2]) + " - " + str(end_date[2]) + " and " + \
                     str(start_date2[2]) + " - " + str(end_date2[2])
        part_name = "_multi_model"


    # Plot histograms of each analysis
    for a in range(len(analysis_str)):
        file_names = defaultdict(list)
        if plot is not None:
            # Plot histogram
            sns.set(rc={'figure.figsize': (11, 11)})
            for i in range(len(variables)):
                if not cov:
                    fig, ax = plt.subplots()
                    title_name = analysis_str[a] + " of " + cubes[a][i].name() + part_title
                    if second_date_given:
                        title_name = "multi_model " + title_name
                    ax.set_title(title_name)
                    if bins_dict is None:
                        if len(full_datum[a][i]) == 0:
                            print("ERROR in function create_histogram: No data obtained for histogram.")
                            return None
                        min_bin, max_bin = min(full_datum[a][i]), max(full_datum[a][i])
                        if min_bin == max_bin:
                            print("ERROR in function create_histogram: min and max bin in histogram are equal.")
                            return None
                        try:
                            ax.hist(full_datum[a][i], bins=sel, range=(min_bin, max_bin), color="purple")
                        except (MemoryError, ValueError):
                            ax.hist(full_datum[a][i], bins=10, range=(min_bin, max_bin), color="purple")
                    else:
                        ax.hist(full_datum[a][i], bins=bins_dict[variables[i]], color="purple")
                    ax.set_ylabel("Frequency")
                    ax.set_xlabel(cubes[a][i].name() + ' (' + str(cubes[a][i].units) + ')')

                    # Get file name
                    file_name = "hist_" + analysis_str[a] + "_" + cubes[a][i].name() + part_name
                    file_name = make_into_file_name(file_name)
                    plt.savefig(os.path.join(directories.ANALYSIS, file_name))
                    file_name = file_name + '.txt'
                    file_names[a].append(file_name)

            if cov:  # Plot 2D histogram
                if len(full_datum[a][0]) == 0 or len(full_datum[a][1]) == 0:
                    print("ERROR in function create_histogram: No data obtained for histogram.")
                    return None
                fig, ax = plt.subplots()
                title_name = analysis_str[a] + " of " + cubes[a][0].name() + " and " + cubes[a][1].name() + part_title
                if second_date_given:
                    title_name = "multi_model " + title_name
                ax.set_title(title_name)
                if bins_dict is None:
                    _, _, _, img = ax.hist2d(full_datum[a][0], full_datum[a][1], bins=sel)
                else:
                    _, _, _, img = ax.hist2d(full_datum[a][var1], full_datum[a][var2], bins=[bins_dict[x], bins_dict[y]])
                ax.set_xlabel(cubes[a][var1].name() + ' (' + str(cubes[a][var1].units) + ')')
                ax.set_ylabel(cubes[a][var2].name() + ' (' + str(cubes[a][var2].units) + ')')
                plt.colorbar(img, ax=ax, label='Frequency')

        if save_out:
            for i in range(len(variables)):
                if not cov:
                    if plot is None:
                        file_name = "hist_" + analysis_str[a] + "_" + cubes[a][i].name() + part_name
                        file_name = make_into_file_name(file_name) + '.txt'
                        file_names[a].append(file_name)
                    file = open(os.path.join(directories.ANALYSIS, file_names[a][i]), 'w')
                    # Write header of file
                    file.write("# " + cubes[a][i].name() + "," + "Frequency\n")
                    # Get histogram data and write to file
                    if bins_dict is None:
                        if len(full_datum[a][i]) == 0:
                            print("ERROR in function create_histogram: No data obtained for histogram.")
                            return None
                        min_bin, max_bin = min(full_datum[a][i]), max(full_datum[a][i])
                        if min_bin == max_bin:
                            print("ERROR in function create_histogram: min and max bin in histogram are equal.")
                            return None
                        hist, bins = None, None
                        try:
                            hist, bins = np.histogram(full_datum[a][i], bins=sel, range=(min_bin, max_bin))
                        except (MemoryError, ValueError):
                            hist, bins = np.histogram(full_datum[a][i], bins=10, range=(min_bin, max_bin))
                    else:
                        hist, bins = np.histogram(full_datum[a][i], bins=bins_dict[variables[i]])
                    data = list(zip(bins, hist))
                    np.savetxt(file, data, delimiter=',')
                    file.close()
            if cov:
                if len(full_datum[a][0]) == 0 or len(full_datum[a][1]) == 0:
                    print("ERROR in function create_histogram: No data obtained for histogram.")
                    return None
                file_name = "hist2D_" + analysis_str[a] + "_" + cubes[a][0].name() + "_" + cubes[a][1].name() + part_name
                file_name = make_into_file_name(file_name)
                plt.savefig(os.path.join(directories.ANALYSIS, file_name))
                file_name = file_name + '.txt'
                file = open(os.path.join(directories.ANALYSIS, file_name), 'w')
                file.write("# " + cubes[a][var1].name() + "," + cubes[a][var2].name() + "," + "Frequency\n")
                if bins_dict is None:
                    hist, xs, ys = np.histogram2d(full_datum[a][0], full_datum[a][1], bins=sel)
                else:
                    hist, xs, ys = np.histogram2d(full_datum[a][var1], full_datum[a][var2],
                                                  bins=[bins_dict[x], bins_dict[y]])

                # Get each corresponding value and write to file
                for i in range(len(xs) - 1):
                    for j in range(len(ys) - 1):
                        file.write(str(xs[i]) + "," + str(ys[j]) + "," + str(hist[i][j]) + "\n")
                file.close()

        print("Histogram data is saved in the " + directories.ANALYSIS + " folder as a txt file.")


def create_timeseries(list_ens, start_date, end_date, variables,  monthly=False, save_out=True, ens_num=1,
                      func_name=None, second_date_given=False, plot=None):
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
    :param second_date_given: if set, multi model averages calculated
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
        cube = list_ens[ens_num-1][variable]
        # Construct title name
        title_name = cube.name() + " measured " + time_str + " between " + str(start_date[2])  + " and " + str(end_date[2])

        if func_name is not None:
            title_name = func_name + " of " + title_name
        if second_date_given and func_name is not None:
            title_name = "Multi model " + title_name
        # Convert cube to pandas series dataframe
        pd_data = None
        try:
            pd_data = as_series(cube)
        except Exception:
            print("WARNING in function create_timeseries: cannot construct timeseries with more than 1-dimension cube.")
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
        data = np.asarray(cube.data)
        df = pd.DataFrame(data=data, index=indx, columns=[variable])

        # Save table to file for each variable
        # Construct file name
        date_str = " " + time_str + " " + str(start_date[2]) + "_" + str(end_date[2])
        file_name = "ts_" + cube.name() + date_str
        file_name = make_into_file_name(file_name)
        if save_out:
            # np.savetxt(file_name, df.values, fmt='%d', comments="dates  " + cube.name())
            df.to_csv(os.path.join(directories.ANALYSIS, file_name  + '.txt'), sep=',', index=True, index_label='dates')
            print("Timeseries data is saved in the " + directories.ANALYSIS + " folder as a txt file.")

        if plot:
            num_years = len(np.unique(indx.year))
            # Plot the timeseries
            fig, axs = plt.subplots(1, 1)
            # Timeseries
            df.plot(ax=axs, title=title_name, grid=True, legend=False, alpha=0.7, color='m', ls='-')
            if not monthly:
                df.resample('BM').mean().plot(ax=axs, style='-')
            rolling_str = 'one-year rolling mean'
            if num_years > 3:
                if monthly:
                    df.rolling(12).mean().plot(ax=axs, style='-')
                else:
                    df.rolling(365).mean().plot(ax=axs, style='-')
            if not monthly and num_years > 3:
                axs.legend(['input', 'monthly mean', rolling_str], loc='upper left')
            elif monthly and num_years > 3:
                axs.legend(['input', rolling_str], loc='upper left')
            elif num_years <= 2 and monthly:
                axs.legend(['input'], loc='upper left')
            elif num_years <= 2:
                axs.legend(['input', 'monthly mean'], loc='upper left')
            axs.set_xlabel("Time (month/year)")
            axs.set_ylabel(cube.name())

            if save_out:
                plt.savefig(os.path.join(directories.ANALYSIS, file_name))
                print("Timeseries plot is saved in the " + directories.ANALYSIS + " folder as a png file.")

            # Boxplot - Yearly seasonality
            fig, axs1 = plt.subplots(1, 1)
            df['Month'] = df.index.month
            sns.boxplot(ax=axs1, data=df, x='Month', y=variable)
            axs1.set_title("Yearly seasonality of " + cube.name())
            axs1.set_xlabel("Months")
            axs1.set_ylabel(cube.name())


def create_timeseries_helper(cube, variable, start_date, end_date, time_str, monthly, title_name, save_out,
                             spatial, second_date_given, plot, analysis, start2, end2):
    """
    Create timeseries helper function
    """

    # Convert cube to pandas series dataframe
    pd_data = None
    try:
        pd_data = as_series(cube)
    except Exception:
        print("WARNING in function create_timeseries: cannot construct timeseries with more than 1-dimension cube.")
        return None

    # Get x values in series
    indices = pd_data.index
    if len(indices) == 1:
        print("WARNING in function create_timeseries: cannot construct timeseries of one value.")
        return None
    # Save dates that are in datetime format
    new_indices = []
    for i in range(len(indices)):
        # Convert cftime.datetime to datetime
        dt = convert_cftime_datetime(indices[i])
        new_indices.append(dt)

    indx = pd.DatetimeIndex(new_indices)
    data = np.asarray(cube.data)
    # df = pd.Series(pd_data.to_numpy(), index=indx)
    df = pd.DataFrame(data=data, index=indx, columns=[variable])

    # Get names of indexes for which column is > 1e+20
    indexNames = df[df[variable] >= 1e+20].index
    # Delete these row indexes from dataFrame
    df.drop(indexNames, inplace=True)
    df.dropna(inplace=True)

    # Save table to file for each variable
    # Construct file name
    date_str = " " + time_str + "_multimodel"
    if not second_date_given:
        date_str = " " + time_str + " "

    file_name = "ts_" + analysis + '_' + cube.name() + date_str + "_" + str(start_date[2]) + "-" + str(end_date[2])
    if second_date_given:
        file_name = file_name + "_" + str(start2[2]) + '-' + str(end2[2]) + '_multi_model'
    file_name = make_into_file_name(file_name)
    if save_out:
        # np.savetxt(file_name, df.values, fmt='%d', comments="dates  " + cube.name())
        df.to_csv(os.path.join(directories.ANALYSIS, file_name + '.txt'), sep=',', index=True, index_label='dates')
        print("Timeseries data is saved in the " + directories.ANALYSIS + " folder as a txt file.")

    if plot:
        # Construct figure and axes
        sns.set(rc={'figure.figsize': (11, 4)})
        num_years = len(np.unique(indx.year))
        # Plot the timeseries
        fig, axs = plt.subplots(1, 1)
        # Timeseries
        df.plot(ax=axs, title=title_name, grid=True, legend=False, alpha=0.7, color='m', ls='-')
        if not monthly:
            df.resample('BM').mean().plot(ax=axs, style='-')
        rolling_str = 'one-year rolling mean'

        if num_years > 3:
            if monthly:
                df.rolling(12).mean().plot(ax=axs, style='-')
            else:
                df.rolling(365).mean().plot(ax=axs, style='-')
        if not monthly and num_years > 3:
            axs.legend(['input', 'monthly mean', rolling_str], loc='upper left')
        elif monthly and num_years > 3:
            axs.legend(['input', rolling_str], loc='upper left')
        elif num_years <= 2 and monthly:
            axs.legend(['input'], loc='upper left')
        elif num_years <= 2:
            axs.legend(['input', 'monthly mean'], loc='upper left')

        if second_date_given:
            axs.set_xlabel("Time (the control year)")
        else:
            axs.set_xlabel("Time (month/year)")
        axs.set_ylabel(cube.name())

        if save_out:
            plt.savefig(os.path.join(directories.ANALYSIS, file_name))
            print("Timeseries plot is saved in the " + directories.ANALYSIS + " folder as a png file.")

        # Boxplot - Yearly seasonality
        fig, axs1 = plt.subplots(1, 1)
        df['Month'] = df.index.month
        sns.boxplot(ax=axs1, data=df, x='Month', y=variable)
        axs1.set_title("Yearly seasonality of " + cube.name())
        axs1.set_xlabel("Months")
        axs1.set_ylabel(cube.name())

        if save_out:
            file_name = "bp_" +  analysis + '_' + cube.name() + date_str
            file_name = make_into_file_name(file_name)
            plt.savefig(os.path.join(directories.ANALYSIS, file_name))
            print("Box plot is saved in the " + directories.ANALYSIS + " folder as a png file.")


def create_ts_both_graphs(calc, calc2, variable, analysis, save_out):
    """
    Used to plot the control and future run in one graph
    :param calc: control cube
    :param calc2: future cube
    :param variable: variable
    :param analysis: analysis string
    :return: Plot
    """
    # Construct figure and axes
    sns.set(rc={'figure.figsize': (11, 4)})

    # Convert cube to pandas series dataframe
    pd_data, pd_data2 = None, None
    try:
        pd_data = as_series(calc)
        pd_data2 = as_series(calc2)
    except Exception:
        print("WARNING in function create_timeseries: cannot construct timeseries with more than 1-dimension cube.")
        return None

    # Get x values in series
    indices, indices2 = pd_data.index, pd_data2.index
    if len(indices) == 1:
        print("WARNING in function create_timeseries: cannot construct timeseries of one value.")
        return None
    # Save dates that are in datetime format
    new_indices, new_indices2 = [], []
    for i in range(len(indices)):
        # Convert cftime.datetime to datetime
        dt = convert_cftime_datetime(indices[i])
        dt2 = convert_cftime_datetime(indices2[i])
        new_indices.append(dt)
        new_indices2.append(dt2)

    indx, indx2 = pd.DatetimeIndex(new_indices), pd.DatetimeIndex(new_indices2)
    # df = pd.Series(pd_data.to_numpy(), index=indx)
    df = pd.DataFrame(data=pd_data.to_numpy(), index=indx, columns=[variable])
    df2 = pd.DataFrame(data=pd_data2.to_numpy(), index=indx2, columns=[variable])

    # Get names of indexes for which column is > 1e+20
    indexNames = df[df[variable] >= 1e+20].index
    # Delete these row indexes from dataFrame
    df.drop(indexNames, inplace=True)
    df.dropna(inplace=True)
    # Get names of indexes for which column is > 1e+20
    indexNames = df2[df2[variable] >= 1e+20].index
    # Delete these row indexes from dataFrame
    df2.drop(indexNames, inplace=True)
    df2.dropna(inplace=True)

    # Plot the timeseries
    fig, axs = plt.subplots(1, 1)
    plt.tight_layout()
    title_name = "Control and future run of spatial " + analysis + " of " + calc.name()

    # Timeseries
    axs = df.plot(title=title_name, grid=True, alpha=0.7, color='b', ls='-')
    df2.plot(ax=axs, color='g', ls='-')
    axs.legend(["Control", "Future"])
    axs.set_xlabel("Time (month/year)")
    axs.set_ylabel(calc.name())

    if save_out:
        file_name = make_into_file_name(title_name)
        plt.savefig(os.path.join(directories.ANALYSIS, file_name))
        print("Timeseries plot is saved in the " + directories.ANALYSIS + " folder as a png file.")


def create_timeseries_analysis(list_ens, start_date, end_date, variables, analysis_str, monthly=False, save_out=True, ens_num=1,
                       second_date_given=False, total=False, spatial=False,
                               calcs=None, calcs2=None, plot=None, start2=None, end2=None):
    """
    Analysis the data given - in this case it computes the timeseries (assumes grid/sample point)
    :param list_ens: the list of ensembles (dicts) containing the data of the climate variables
    :param start_date extract from end date from data, list [d, m, y]
    :param end_date: extract till end date from data, list [d, m, y]
    :param variables: If one variable - then have 1D histogram. If list of 2 variables, then 2D histogram
    :param analysis_str: list of analysis performed on models
    :param monthly: data is stored in monthly increments (time = 12) else assumed (time = 365)
    :param total: set if all ensemble are averaged into one
    :param second_date_given: if set, multi model averages calculated
    :param spatial: spatial averages have been calculated
    :param calc: control cube
    :param calc2: future cube
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

    # Convert pandas dates to matplotlib date format
    register_matplotlib_converters()

    for variable in variables:
        for a in range(len(analysis_str)):
            # Extract time in cube
            cube, calc, calc2 = None, None, None
            if total:
                cube = list_ens[a][variable]
                if spatial and second_date_given:
                    calc = calcs[a][variable]
                    calc2 = calcs2[a][variable]
            else:
                cube = list_ens[ens_num-1][a][variable]
                if spatial and second_date_given:
                    calc = calcs[ens_num-1][a][variable]
                    calc2 = calcs2[ens_num-1][a][variable]

            # Construct title name
            title_name =  "spatial " + str(analysis_str[a]) + " of " + cube.name() + " measured " + time_str + " between " + str(start_date[2]) + " and " + str(
                end_date[2])

            if second_date_given:
                title_name = "Multi model " + "spatial " + str(analysis_str[a]) + " of " + cube.name() + " measured " + \
                             time_str + " between " + str(start_date[2]) + "-" + str(end_date[2]) + " and " + \
                             str(start2[2]) + "-" + str(end2[2])

            if len(cube.dim_coords) == 0:
                iris.util.promote_aux_coord_to_dim_coord(cube, 'time')
            if second_date_given and spatial:
                if len(calc.dim_coords) == 0:
                    iris.util.promote_aux_coord_to_dim_coord(calc, 'time')
                if len(calc2.dim_coords) == 0:
                    iris.util.promote_aux_coord_to_dim_coord(calc2, 'time')

            create_timeseries_helper(cube, variable, start_date, end_date, time_str, monthly, title_name, save_out,
                                     spatial, second_date_given, plot, str(analysis_str[a]), start2, end2)

            if second_date_given and plot is not None:  # Plot both graphs in one
                create_ts_both_graphs(calc, calc2, variable, str(analysis_str[a]), save_out)


def plot_map_helper(cube, title_name, save_out=False):
    """
    Plot (and save) maps given date and variable in cube
    :param cube: data
    :param title_name: title of graph
    :param save_out: if set, then saves to png file
    """

    # Plot the graph
    fig, ax = plt.subplots(figsize=(11, 4))
    # Get the Purples "Brewer" palette for plot
    brewer_cmap = plt.get_cmap('brewer_Purples_09')

    # Draw the contours, with n-levels set for the map colours
    lats = cube.coord('latitude')
    lons = cube.coord('longitude')
    masked_array = np.ma.array(cube.data, mask=np.isnan(cube.data))

    # Get max and min of lons and lats
    lon_min, lon_max, lat_min, lat_max = min(lons.points), max(lons.points), min(lats.points), max(lats.points)

    # Only set if one lat/lon
    frame, cf = False, None

    try:
        # Overlay world map
        m = Basemap(projection='mill', lat_ts=10, resolution='c', ax=ax)
        # Plot actual data
        scat = False
        if len(lons.points) < 30 or len(lats.points) < 30:
            scat = True
            if len(lons.points) == 1 and len(lats.points) == 1:
                frame = True
                label = "Point (" + str(lons.points[0]) + ", " + str(lons.points[0]) + ")"
                m.scatter(lons.points, lats.points, latlon=True, s=20, c='blue', marker='o', label=label)
            else:
                m.scatter(lons.points, lats.points, latlon=True, s=20, c='blue', marker='o')
        else:
            # Overlay world map
            m = Basemap(projection='mill', lat_ts=10, llcrnrlon=lon_min, \
                        urcrnrlon=lon_max, llcrnrlat=lat_min, urcrnrlat=lat_max, \
                        resolution='c', ax=ax)
            x, y = m(*np.meshgrid(lons.points, lats.points))
            cf = m.contourf(x, y, masked_array, shading='flat', cmap=brewer_cmap)
        m.drawcoastlines()
        m.drawparallels(np.arange(-90., 90., 30.), labels=[1, 0, 0, 0], ax=ax)
        m.drawmeridians(np.arange(-180., 180., 60.), labels=[0, 0, 0, 1], ax=ax)
    except Exception as err:
        print("ERROR in function plot_map: " + str(err))
        plt.close(fig)
        return None
    if not scat:
        fig.colorbar(cf, ax=ax, label=cube.units)

    # Add a citation to the plot.
    iplt.citation(iris.plot.BREWER_CITE)

    ax.set_title(title_name)

    # Have legend if point
    if frame:
        plt.legend(frameon=True)

    if save_out:
        title_name = make_into_file_name(title_name)
        plt.savefig(os.path.join(directories.ANALYSIS, title_name))
        print("Map is saved in the " + directories.ANALYSIS + " folder as a png file.")




def plot_map_analysis(list_ens, variables, analysis_str=None, ens_num=1, save_out=False, total=False,
             second_date_given=False):
    """
    Plot (and save) maps given date and variable in cube
    :param list_ens: the list of ensembles (dicts) containing the data of the climate variables
    :param variables: variables to select in cube
    :param analysis_str: list of analysis names performed on ensemble, list of strings
    :param ens_num: gives which ensemble to plot, default = 1
    :param save_out: if set, then saves to png file
    :param total: set if all ensembles have been calculated together, instead of separately, boolean
    :param second_date_given: set if multi model stats are calculated
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
            if total:
                cube = list_ens[a][variable]
            else:
                cube = list_ens[ens_num-1][a][variable]

            if cube.ndim != 2:
                print("WARNING in function plot_map: Data must be 2D to plot map.")
                return None

            # Plot title name
            title_name = analysis_str[a] + " of " + cube.name() + " of ensemble " + str(ens_num) + " with variable " + str(variable)
            if total:
                title_name = analysis_str[a] + " of " + cube.name() + " of all ensembles" + " with variable " + str(variable)
            if second_date_given:
                title_name = "Multi model " + title_name

            plot_map_helper(cube, title_name, save_out=save_out)


def plot_map(list_ens, variables, ens_num=1, save_out=False, total=False, second_date_given=False):
    """
    Plot (and save) maps given date and variable in cube
    :param list_ens: the list of ensembles (dicts) containing the data of the climate variables
    :param variables: variables to select in cube
    :param ens_num: gives which ensemble to plot, default = 1
    :param save_out: if set, then saves to png file
    :param total: set if all ensembles have been calculated together, instead of separately, boolean
    :param second_date_given: set if multi model stats are calculated
    """
    # Make sure data structures are not empty
    assert list_ens is not None
    assert variables is not None

    # If variables is just one object, cast to list
    if not isinstance(variables, list):
        variables = [variables]

    for variable in variables:
        # Extract time in cube
        cube = None
        if total:
            cube = list_ens[0][variable]
        else:
            cube = list_ens[ens_num-1][variable]

            if cube.ndim != 2:
                print("WARNING in function plot_map: Data must be 2D to plot map.")
                return None

        # Plot title name
        title_name = cube.name() + " of ensemble " + str(ens_num) + " with variable " + str(variable)
        if total:
            title_name = cube.name() + " of all ensembles" + " with variable " + str(variable)
        if second_date_given:
            title_name = "Multi model " + title_name

        plot_map_helper(cube, title_name, save_out=save_out)