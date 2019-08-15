import iris
import importlib
import directories
import numpy as np
from cdo import Cdo
import xarray as xr


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

    nan_indices, indices_filled = None, False

    if analysis:
        analysis = [a.lower() for a in analysis]
    elif total and not analysis:
        # Get variables from one ensemble
        variables = list(list_ens[0])
        # Save averages
        ens_calcs, dict_vars = [], {}
        for var in variables:
            unit = list_ens[0][var].units
            long_name = list_ens[0][var].long_name
            # convert cubes to xarray
            xr_var = [xr.DataArray.from_iris(ens[var]) for ens in list_ens]
            # Combine the dataarrays
            mg = xr.concat(xr_var, 'ensemble_member')
            # Find average across ensembles
            av_mg = mg.mean(dim='ensemble_member', skipna=True)
            # Convert back to cube
            cube = xr.DataArray.to_iris(av_mg)
            cube.units = unit
            cube.long_name = long_name

            dict_vars[var] = cube
            # Get nan indices
            if not indices_filled:
                if np.ma.is_masked(cube.data):
                    x = cube.data.filled()
                    nan_indices = np.argwhere(np.isclose(x.flatten().data, cube.data.fill_value))
                else:
                    nan_indices = np.argwhere(np.isnan(cube.data.flatten().data))
                break
                indices_filled = True
        ens_calcs.append(dict_vars)
        return ens_calcs, False, nan_indices
    elif not total and not analysis:
        return None, False, None

    # Rearrange analysis in this order
    a_order = ['mean', 'std', 'median', 'rms']
    analysis_new_order = []
    for o in a_order:
        if o in analysis:
            analysis_new_order.append(o)

    analysis = analysis_new_order

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
                # Get nan indices
                if not indices_filled:
                    for s in dict_[d].slices_over('time'):
                        if np.ma.is_masked(s.data):
                            x = s.data.filled()
                            nan_indices = np.argwhere(np.isclose(x.flatten().data, s.data.fill_value))
                        else:
                            nan_indices = np.argwhere(np.isnan(s.data.flatten().data))
                        break
                    indices_filled = True
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

    # Find ensembles average
    if total:
        mean_calcs, rms_calcs, median_calcs, std_calcs = {}, {}, {}, {}
        # Get variables from one ensemble
        variables = list(list_ens[0])
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

    return ens_calcs, analysis, nan_indices


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
    pkg = directories.INPUT + '.' + directories.USER_FUNCTION_PACKAGE

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



