import iris
from utils import *


def compute_stats_analysis(list_ens, time_name, analysis='mean'):
    """
       Analysis the data given - in this case it computes the mean, std, median and rms
       :param list_ens: the list of ensembles (dicts) containing the data of the climate variables
       :param time_name: name of time variable in original file
       :param analysis: type of computation
       :return:
           ens_means: list of averages of the different ensembles
       """
    # Holds the means for each ensemble
    ens_calcs = []
    for dict_ in list_ens:
        # Save the mean of each variable in a dict of list
        calcs = {}
        # Calculate the mean of each variable in the dictionary given
        for d in dict_:
            calc = None
            if analysis == 'mean':
                calc = dict_[d].collapsed(time_name, iris.analysis.MEAN).data
            elif analysis == 'std':
                calc = dict_[d].collapsed(time_name, iris.analysis.STD_DEV).data
            elif analysis == 'median':
                calc = dict_[d].collapsed(time_name, iris.analysis.MEDIAN).data
            elif analysis.lower() == 'rms':
                calc = dict_[d].collapsed(time_name, iris.analysis.RMS).data
            calcs[d] = calc
        ens_calcs.append(calcs)

    return ens_calcs


def compute_average(list_ens, nan_values):
    """
    Analysis the data given - in this case it computes the mean
    :param list_ens: the list of ensembles (dicts) containing the data of the climate variables
    :param nan_values: missing values in data set
    :return:
        ens_means: list of averages of the different ensembles
    """

    # Holds the means for each ensemble
    ens_means = []
    for dict_ in list_ens:
        # Save the mean of each variable in a dict of list
        means = {}
        # Calculate the mean of each variable in the dictionary given
        for d in dict_:
            # Select the parts of the data within timeframe
            mean = np.mean(dict_[d], axis=0)
            # Replace values close to nan values to actual nan values
            if mean.shape:
                mean[np.isclose(mean, nan_values[d], rtol=1)] = nan_values[d]
            # Save mean for variable
            means[d] = mean
        ens_means.append(means)

    return ens_means

