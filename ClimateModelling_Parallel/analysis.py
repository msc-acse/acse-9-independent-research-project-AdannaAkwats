import iris
from cdo import Cdo
import importlib
import directories
import numpy as np
from multiprocessing import Pool
from functools import partial
import parallel_settings


def compute_stats_parallel(total, time_name, indices_filled, analysis, ens):

    nan_indices, variables = None, []
    # Used only if we analyse all at the same time
    mean_calcs, std_calcs, median_calcs, rms_calcs = {}, {}, {}, {}
    # Calculate the mean of each variable in the dictionary given
    for d in ens:
        variables.append(d)
        for a in analysis:
            if a == 'mean':
                mean_calc = ens[d].collapsed(time_name, iris.analysis.MEAN)
                if total:
                    if not indices_filled:
                        for s in ens[d].slices_over('time'):
                            if np.ma.is_masked(s.data):
                                x = s.data.filled()
                                nan_indices = np.argwhere(np.isclose(x.flatten().data, s.data.fill_value))
                            else:
                                nan_indices = np.argwhere(np.isnan(s.data.flatten().data))
                            break
                        indices_filled = True
                mean_calcs[d] = mean_calc
            if a == 'std':
                std_calc = ens[d].collapsed(time_name, iris.analysis.STD_DEV)
                std_calcs[d] = std_calc
            if a == 'median':
                median_calc = ens[d].collapsed(time_name, iris.analysis.MEDIAN)
                median_calcs[d] = median_calc
            if a == 'rms':
                rms_calc = ens[d].collapsed(time_name, iris.analysis.RMS)
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

    return calcs, nan_indices, variables


def compute_stats_analysis(list_ens, analysis, total=False):
    """
    Analyse the data given - in this case it computes the mean, std, median and rms
    :param list_ens: the list of ensembles (dicts) containing the data of the climate variables
    :param analysis: type of computation
    :param total: if total set, then analyse all ensembles together, else separately, boolean
    :return:
       ens_means:
            - averages of ensembles together (if total)
            - list of averages of the different ensembles (if not total)
       string : type of analysis computed
       """
    # Assertions
    assert list_ens is not None

    analysis = [a.lower() for a in analysis]

    indices_filled = False

    time_name = 'time'

    pool = Pool(processes=parallel_settings.NUM_PROCESSORS)
    # Go through variables in each ensemble

    func = partial(compute_stats_parallel, total, time_name, indices_filled, analysis)

    r = pool.map(func, list_ens)

    pool.close()
    pool.join()

    # Get results from processes
    ens_calcs = [tup[0] for tup in r]
    nan_indices = [tup[1] for tup in r]
    variables = r[0][2]

    # Find ensembles average
    if total:
        mean_calcs, rms_calcs, median_calcs, std_calcs = {}, {}, {}, {}
        for j in range(len(analysis)):
            for v in variables:
                each_a = [t[j][v] for t in ens_calcs]
                av_cube = each_a[0]
                each_a_data = [t[j][v].data for t in ens_calcs]
                av = np.mean(each_a_data, axis=0)
                shape = av.shape
                av_flat = av.flatten()
                np.put(av_flat, nan_indices, np.nan)
                av_flat = av_flat.reshape(shape)
                av = av_flat
                # Get cube and change the data to total averages
                av_cube.data = av

                if analysis[j] == 'mean':
                    mean_calcs[v] = av_cube
                if analysis[j] == 'std':
                    std_calcs[v] = av_cube
                if analysis[j] == 'median':
                    median_calcs[v] = av_cube
                if analysis[j] == 'rms':
                    rms_calcs[v] = av_cube
        # define new ens_calcs to store totals
        ens_calcs = []
        if 'mean' in analysis:
            ens_calcs.append(mean_calcs)
        if 'std' in analysis:
            ens_calcs.append(std_calcs)
        if 'median' in analysis:
            ens_calcs.append(median_calcs)
        if 'rms' in analysis:
            ens_calcs.append(rms_calcs)

    print("function compute_stats_analysis: Averages of data successfully computed.")

    return ens_calcs, analysis


def compute_enso_indices():
    # TODO
    cdo = Cdo()
    print(cdo.operators)
    return None

#
# def compute_user_analysis(list_ens, file_name, func_name):
#     """
#     Use function given by user for analysis
#     :param list_ens: the list of ensembles (dicts) containing the data of the climate variables
#     :param file_name: name of python script where function is
#     :param func_name: name of function to call
#     :param args: other arguments (excluding cube)
#     :return: analysed data of each variable and ensemble
#     """
#     # Assertions
#     assert list_ens is not None
#     assert file_name is not None
#     assert func_name is not None
#
#     # Get folder name
#     pkg = directories.USER_FUNCTION_PACKAGE
#
#     # If .py is at the end of file name, remove it
#     if '.py' in file_name:
#         file_name = file_name.replace('.py', '')
#
#     # Construct module to import
#     module = pkg + '.' + file_name
#
#     # Call user function
#     user_script = importlib.import_module(module)
#     user_func = getattr(user_script, func_name)
#
#     # Call and return user function
#     res = user_func(list_ens)
#
#     print("function compute_user_analysis: User analysis of data successfully computed.")
#
#     return res
#

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
    user_func = getattr(user_script, func_name)

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



