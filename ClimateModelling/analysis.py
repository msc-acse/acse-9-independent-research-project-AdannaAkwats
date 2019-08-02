import iris
from cdo import Cdo
import importlib
import directories


def compute_stats_analysis(list_ens, analysis='mean'):
    """
       Analyse the data given - in this case it computes the mean, std, median and rms
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
                mean_calc = dict_[d].collapsed(time_name, iris.analysis.MEAN)
                std_calc = dict_[d].collapsed(time_name, iris.analysis.STD_DEV)
                median_calc = dict_[d].collapsed(time_name, iris.analysis.MEDIAN)
                rms_calc = dict_[d].collapsed(time_name, iris.analysis.RMS)
            elif analysis.lower() == 'mean':
                calc = dict_[d].collapsed(time_name, iris.analysis.MEAN)
            elif analysis.lower() == 'std':
                calc = dict_[d].collapsed(time_name, iris.analysis.STD_DEV)
            elif analysis.lower() == 'median':
                calc = dict_[d].collapsed(time_name, iris.analysis.MEDIAN)
            elif analysis.lower() == 'rms':
                calc = dict_[d].collapsed(time_name, iris.analysis.RMS)
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
    # TODO
    cdo = Cdo()
    print(cdo.operators)
    return None


def compute_user_analysis(list_ens, file_name, func_name):
    """
    Use function given by user for analysis
    :param list_ens: the list of ensembles (dicts) containing the data of the climate variables
    :param file_name: name of python script where function is
    :param func_name: name of function to call
    :param args: other arguments (excluding cube)
    :return: analysed data of each variable and ensemble
    """
    # Assertions
    assert list_ens is not None
    assert file_name is not None
    assert func_name is not None

    # Get folder name
    pkg = directories.USER_FUNCTION_PACKAGE

    # If .py is at the end of file name, remove it
    if '.py' in file_name:
        file_name = file_name.replace('.py', '')

    # Construct module to import
    module = pkg + '.' + file_name

    # Call user function
    user_script = importlib.import_module(module)
    user_func = getattr(user_script, "max_lon_lat")

    # Holds the user analysis for each ensemble
    ens_calcs = []
    for dict_ in list_ens:
        # Save the analysis of each variable in a dict of list
        calcs = {}
        # Calculate the mean of each variable in the dictionary given
        for d in dict_:
            # Save values to dicts
            calcs[d] = user_func(dict_[d])
        # Save for each ensemble
        ens_calcs.append(calcs)
    print("function compute_user_analysis: User analysis of data successfully computed.")

    return ens_calcs



