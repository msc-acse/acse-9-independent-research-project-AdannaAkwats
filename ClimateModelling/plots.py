import matplotlib.pyplot as plt
from utils import *


def create_histogram(list_ens, units, start_date, end_date, nan_values, monthly=False, save_out=None,
                     cov=None, sel=None, plot=False):
    """
    Analysis the data given - in this case it computes the histogram (assumes grid/sample point)
    :param list_ens: the list of ensembles (dicts) containing the data of the climate variables
    :param units: the units matching to each variable
    :param start_date: start date given to user
    :param end_date: end date given to user
    :param nan_values: the missing values in the variables data array
    :param monthly: data is stored in monthly increments (time = 12) else assumed (time = 365)
    :param save_out: if set, then save output of histogram/ rimeseries
    :param cov: if set, then perform covariance analysis
    :param sel: selection option for bin siz, default is fd - Freedman Diaconis Estimator
    :param plot: if plot is true, then shows plot of histogram
    :return: None
    """

    if not sel:
        sel = 'fd'
    fig, axs = plt.subplots(len(list_ens), len(list_ens[0]), squeeze=False)
    time_str = "daily"
    if monthly:
        time_str = "monthly"
    fig.suptitle("Variables " + str(list(list_ens[0])) + " measured " + time_str + " between " + str(start_date[0]) +
                 "-" + str(start_date[1]) + "-" + str(start_date[2]) + " and " + str(end_date[0]) + "-" +
                 str(end_date[1]) + "-" + str(end_date[2]) + " using the E2S2M climate model")
    a, e = 0, 0
    for dict_ in list_ens:
        axs[e, a].set_title("Ensemble " + str(e))
        for d in dict_:
            ens = dict_[d].flatten()
            indices = np.argwhere(np.isclose(ens, nan_values[d]))
            ens = np.delete(ens, indices)
            hist, bin_edges = np.histogram(ens, bins=sel)
            print(ens)
            if plot and not cov:
                axs[e, a].hist(ens, bins=sel)
                axs[e, a].set_ylabel("Frequency")
                axs[e, a].set_xlabel(d + ' (' + units[d] + ')')
                # a += 1

            if plot and cov:  # If covariance between 2 variables, plot a 2d histogram
                axs[a].hist2d(ens, bins=sel)
        e += 1

    plt.show()

    return None


def create_timeseries(list_ens, units, start_date, end_date, monthly=False, save_out=None, cov=None):
    """
    Analysis the data given - in this case it computes the timeseries (assumes grid/sample point)
    :param list_ens: the list of ensembles (dicts) containing the data of the climate variables
    :param units: the units matching to each variable
    :param start_date and end_date: extract data within this time frame
    :param monthly: data is stored in monthly increments (time = 12) else assumed (time = 365)
    :param save_out: if set, then save output of histogram/ rimeseries
    :param cov: if set, then perform covariance analysis
    :return: None
    """
    return None

#
#
# def plot_graph(file):
#     """
#     Plot the data given in file
#     :param file: file from analysis
#     :return: graph plot in output + saves as png file
#     """
#
#     # Open file
#     dataset = Dataset(file, 'r')
#
#     d = np.array(dataset.variables['mean_air_temperature'])
#
#     fig, ax = plt.subplots()
#     ax.imshow(d)
#
#     png_file = file.rstrip('nc') + 'png'
#     fig.savefig(png_file)
#
#     print("Image is saved in the " + directories.ANALYSIS + " folder as a png file.")
#
#     return None
