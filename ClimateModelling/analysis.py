import iris
from cdo import Cdo
import importlib
import directories


def compute_stats_analysis(list_ens, analysis):
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

    analysis = [a.lower() for a in analysis]

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
            for a in analysis:
                if a == 'mean':
                    mean_calc = dict_[d].collapsed(time_name, iris.analysis.MEAN)
                    mean_calcs[d] = mean_calc
                if a == 'std':
                    std_calc = dict_[d].collapsed(time_name, iris.analysis.STD_DEV)
                    std_calcs[d] = std_calc
                if a == 'median':
                    median_calc = dict_[d].collapsed(time_name, iris.analysis.MEDIAN)
                    median_calcs[d] = median_calc
                if a == 'rms':
                    rms_calc = dict_[d].collapsed(time_name, iris.analysis.RMS)
                    rms_calcs[d] = rms_calc
        calcs = []
        if 'mean' in analysis:
            calcs.append(mean_calcs)
        if 'std' in analysis:
            calcs.append(std_calcs)
        if 'median' in analysis:
            calcs.append(median_calcs)
        if 'rms' in analysis:
            calcs.append(rms_calcs)
        ens_calcs.append(calcs)
    print("function compute_stats_analysis: Averages of data successfully computed.")

    return ens_calcs, analysis


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



