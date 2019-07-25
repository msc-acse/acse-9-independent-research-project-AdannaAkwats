import iris
from cdo import Cdo


def compute_stats_analysis(list_ens, analysis='mean'):
    """
       Analysis the data given - in this case it computes the mean, std, median and rms
       :param list_ens: the list of ensembles (dicts) containing the data of the climate variables
       :param analysis: type of computation
       :return:
           ens_means: list of averages of the different ensembles
           string : type of analysis computed
       """
    # Assertions
    assert list_ens is not None
    assert analysis.lower() in ['mean', 'rms', 'std', 'median', 'all']

    time_name = 'time'
    # Holds the means for each ensemble
    ens_calcs = []
    for dict_ in list_ens:
        # Save the mean of each variable in a dict of list
        calcs = {}
        # Used only if we analyse all at the same time
        mean_calcs, std_calcs, median_calcs, rms_calcs = {}, {}, {}, {}
        # Calculate the mean of each variable in the dictionary given
        for d in dict_:
            calc, mean_calc, std_calc, median_calc, rms_calc = None, None, None, None, None
            if analysis == 'all':  # Calculate all stats
                mean_calc = dict_[d].collapsed(time_name, iris.analysis.MEAN).data
                std_calc = dict_[d].collapsed(time_name, iris.analysis.STD_DEV).data
                median_calc = dict_[d].collapsed(time_name, iris.analysis.MEDIAN).data
                rms_calc = dict_[d].collapsed(time_name, iris.analysis.RMS).data
            elif analysis.lower() == 'mean':
                calc = dict_[d].collapsed(time_name, iris.analysis.MEAN).data
            elif analysis.lower() == 'std':
                calc = dict_[d].collapsed(time_name, iris.analysis.STD_DEV).data
            elif analysis.lower() == 'median':
                calc = dict_[d].collapsed(time_name, iris.analysis.MEDIAN).data
            elif analysis.lower() == 'rms':
                calc = dict_[d].collapsed(time_name, iris.analysis.RMS).data
            # Save values to dicts
            if analysis.lower() == 'all':
                mean_calcs[d] = mean_calc
                std_calcs[d] = std_calc
                median_calcs[d] = median_calc
                rms_calcs[d] = rms_calc
            else:
                calcs[d] = calc
        if analysis.lower() == 'all':
            ens_calcs.append([mean_calcs, std_calcs, median_calcs, rms_calcs])
        else:
            ens_calcs.append(calcs)
    print("function compute_stats_analysis: Averages of data successfully computed.")

    return ens_calcs, analysis.lower()


def compute_enso_indices():
    cdo = Cdo()
    print(cdo.operators)


