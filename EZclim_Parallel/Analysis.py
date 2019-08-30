"""
Created by Adanna Akwataghibe (Github: AdannaAkwats)
"""
import iris
import importlib
import directories
import numpy as np
import xarray as xr
import iris.analysis.maths as iam
import sys
from multiprocessing import Pool
from functools import partial
import parallel_settings


class Analysis:
    def __init__(self, list_ens):
        self.list_ens = list_ens

    def calculate_avg_ensembles(self, analysis):
        """
        Averages all ensembles into one ensembles if no analysis given
        :param list_ens: the list of ensembles (dicts) containing the data of the climate variables
        :param analysis: should be False
        :return:
           ens_calcs:
                - averages of ensembles together
           string : type of analysis computed
           nan indices: indices where the cube is nan
           """
        # Assertions
        assert self.list_ens is not None
        assert analysis is False

        nan_indices = {}

        # Get variables from one ensemble
        variables = list(self.list_ens[0])
        # Save averages
        ens_calcs, dict_vars = [], {}
        for var in variables:
            unit = self.list_ens[0][var].units
            long_name = self.list_ens[0][var].long_name
            var_name = self.list_ens[0][var].var_name
            s_name = self.list_ens[0][var].standard_name
            atrr = self.list_ens[0][var].attributes
            # convert cubes to xarray
            xr_var = [xr.DataArray.from_iris(ens[var]) for ens in self.list_ens]
            # Combine the dataarrays
            mg = xr.concat(xr_var, 'ensemble_member')
            # Find average across ensembles
            av_mg = mg.mean(dim='ensemble_member', skipna=True)
            # Convert back to cube
            cube = xr.DataArray.to_iris(av_mg)
            cube.units = unit
            cube.long_name = long_name
            cube.var_name = var_name
            cube.attributes = attr
            cube.standard_name = s_name

            dict_vars[var] = cube
            # Get nan indices
            if np.ma.is_masked(cube.data):
                x = cube.data.filled()
                nan_indices[var] = np.argwhere(np.isclose(x.flatten().data, cube.data.fill_value))
            else:
                nan_indices[var] = np.argwhere(np.isnan(cube.data.flatten().data))
            break
        ens_calcs.append(dict_vars)

        print("function calculate_avg_ensembles: Averages of data successfully computed.")

        return ens_calcs, False, nan_indices


    def get_nan_indices(self):
        """
        Calculate the nan indices of one time step of list_ens cube
        :param list_ens: list of ensembles
        :return: list of indices where cube is nan
        """
        nan_indices = {}

        # Get variables from one ensemble
        variables = list(self.list_ens[0])

        # Get nan indices
        for var in variables:
            cube = self.list_ens[0][var]
            for s in cube.slices_over('time'):
                if np.ma.is_masked(s.data):
                    x = s.data.filled()
                    nan_indices[var] = np.argwhere(np.isclose(x.flatten().data, s.data.fill_value))
                else:
                    nan_indices[var] = np.argwhere(np.isnan(s.data.flatten().data))
                break
        return nan_indices


    def compute_total_stats_analysis(self, analysis, nan_indices, spatial, ens_calcs):
        """
        Analyse the data given - in this case it computes the mean, std, median and rms
        :param ens_calcs: the list of averaged ensembles (dicts) containing the data of the climate variables
        :param analysis: type of computation
        :param total: if total set, then analyse all ensembles together, else separately, boolean
        :return:
           ens_means:
                - averages of ensembles together (if total)
           string : type of analysis computed
           list: nan indices
           """
        # Assertions
        assert self.list_ens is not None

        # Get variables from one ensemble
        variables = list(self.list_ens[0])

        # Find ensembles average
        mean_calcs, rms_calcs, median_calcs, std_calcs = {}, {}, {}, {}
        for j in range(len(analysis)):
            for v in variables:
                each_a = [t[j][v] for t in ens_calcs]
                av_cube = each_a[0]
                each_a_data = [t[j][v].data for t in ens_calcs]
                av = np.mean(each_a_data, axis=0)
                if not spatial:
                    shape = av.shape
                    av_flat = av.flatten()
                    np.put(av_flat, nan_indices[v], np.nan)
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

        print("function compute_total_stats_analysis: Averages of data successfully computed.")

        return ens_calcs, analysis, nan_indices


    def reorder_analysis_str(self, analysis):
        """
        Rearrange analysis list to  ['mean', 'std', 'median', 'rms']
        :param analysis: analysis list
        :return: reordered list
        """
        # Rearrange analysis in this order
        a_order = ['mean', 'std', 'median', 'rms']
        analysis_new_order = []
        for o in a_order:
            if o in analysis:
                analysis_new_order.append(o)
        return analysis_new_order


    def compute_stats_parallel(self, time_name, analysis, spatial, spatial_names, dict_):
        """
        Helper function for parallelisation
        """
        mean_calcs, std_calcs, median_calcs, rms_calcs = {}, {}, {}, {}
        # Calculate the mean of each variable in the dictionary given
        for d in dict_:
            for a in analysis:
                if a == 'mean':
                    mean_calc = None
                    if spatial:
                        mean_calc = dict_[d].collapsed(spatial_names, iris.analysis.MEAN, mdtol=0)
                    else:
                        mean_calc = dict_[d].collapsed(time_name, iris.analysis.MEAN, mdtol=0)
                    mean_calcs[d] = mean_calc
                if a == 'std':
                    std_calc = None
                    if spatial:
                        std_calc = dict_[d].collapsed(spatial_names, iris.analysis.STD_DEV, mdtol=0)
                    else:
                        std_calc = dict_[d].collapsed(time_name, iris.analysis.STD_DEV, mdtol=0)
                    std_calcs[d] = std_calc
                if a == 'median':
                    median_calc = None
                    if spatial:
                        median_calc = dict_[d].collapsed(spatial_names, iris.analysis.MEDIAN, mdtol=0)
                    else:
                        median_calc = dict_[d].collapsed(time_name, iris.analysis.MEDIAN, mdtol=0)
                    median_calcs[d] = median_calc
                if a == 'rms':
                    rms_calc = None
                    if spatial:
                        rms_calc = dict_[d].collapsed(spatial_names, iris.analysis.RMS, mdtol=0)
                    else:
                        rms_calc = dict_[d].collapsed(time_name, iris.analysis.RMS, mdtol=0)
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

        return calcs


    def compute_stats_analysis(self, analysis, total=False, spatial=False, dim_coords=None):
        """
        Analyse the data given - in this case it computes the mean, std, median and rms
        :param list_ens: the list of ensembles (dicts) containing the data of the climate variables
        :param analysis: type of computation
        :param total: if total set, then analyse all ensembles together, else separately, boolean
        :param spatial: set if spatial avergaes will be calculated
        :param dim_coords: dictionary containing original (long) names of dimensions
        :return: average, nan_indices, analysis
        """
        # Assertions
        assert self.list_ens is not None

        if analysis:
            analysis = [a.lower() for a in analysis]
        elif total and not analysis:
            return calculate_avg_ensembles(analysis)
        elif not total and not analysis:
            return None, False, None

        # Get nan indices
        nan_indices = self.get_nan_indices()

        # Rearrange analysis
        analysis = self.reorder_analysis_str(analysis)

        time_name = 'time'
        # If spatial, then get dim names
        spatial_names = None
        if spatial:
            if len(dim_coords) == 3:
                spatial_names = ['latitude', 'longitude']
            else:
                spatial_names = ['latitude', 'longitude', dim_coords['depth']]

        # Assign process to each ensemble
        pool = Pool(processes=parallel_settings.NUM_PROCESSORS)
        func = partial(self.compute_stats_parallel, time_name, analysis, spatial, spatial_names)
        ens_calcs = pool.map(func, self.list_ens)
        pool.close()
        pool.join()

        # Find ensembles average
        if total:
            return self.compute_total_stats_analysis(analysis, nan_indices, spatial, ens_calcs)

        print("function compute_stats_analysis: Averages of data successfully computed.")

        return ens_calcs, analysis, nan_indices


    def compute_user_analysis(self, file_name, func_name):
        """
        Use function given by user for analysis
        :param list_ens: the list of ensembles (dicts) containing the data of the climate variables
        :param file_name: name of python script where function is
        :param func_name: name of function to call
        :param args: other arguments (excluding cube)
        :return: analysed data of each variable and ensemble
        """
        # Assertions
        assert self.list_ens is not None
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
        for dict_ in self.list_ens:
            # Pass in dictionary of vsriable for user function , get back changed
            calcs, name, long_name, unit = user_func(dict_)
            # Convert to cube
            try:
                calcs.var_name
            except Exception:  # Not a cube
                calcs = iris.cube.Cube(calcs, var_name=name, long_name=long_name, units=unit)
            # Save for each ensemble
            ens_calcs.append(calcs)
        print("function compute_user_analysis: User analysis of data successfully computed.")

        return ens_calcs


    def calc_stats_difference(self, list_ens2, analysis, total=False, spatial=False, dim_coords=None):
        """
        Calculate multi model analysis between models list_ens and list_ens2
        :param list_ens: the list of ensembles of 1st model(dicts) containing the data of the climate variables
        :param list_ens: the list of ensembles od 2nd model (dicts) containing the data of the climate variables
        :param analysis: type of computation
        :param total: if total set, then analyse all ensembles together, else separately, boolean
        :param spatial: set if spatial avergaes will be calculated
        :param dim_coords: dictionary containing original (long) names of dimensions
        :return: average difference, individual average (if spatial set), nan_indices, analysis
        """

        # Assertions
        assert self.list_ens is not None
        assert list_ens2 is not None

        if analysis:
            analysis = [a.lower() for a in analysis]
        elif not analysis:
            print("ERROR in function calc_stats_difference: No analysis argument given.")
            sys.exit()

        # Get nan indices
        nan_indices = self.get_nan_indices()

        # Rearrange analysis
        analysis = self.reorder_analysis_str(analysis)

        time_name = 'time'

        # If spatial, then get dim names
        spatial_names = None
        if spatial:
            if len(dim_coords) == 3:
                spatial_names = ['latitude', 'longitude']
            else:
                spatial_names = ['latitude', 'longitude', dim_coords['depth']]

        # Holds the means for each ensemble
        ens_calcs, spat_calcs, spat_calcs2 = [], [], []
        for i in range(len(self.list_ens)):
            dict_ = self.list_ens[i]
            dict2_ = list_ens2[i]
            # Used only if we analyse all at the same time
            mean_calcs, std_calcs, median_calcs, rms_calcs = {}, {}, {}, {}
            mean_sp, std_sp, median_sp, rms_sp = {}, {}, {}, {}
            mean_sp2, std_sp2, median_sp2, rms_sp2 = {}, {}, {}, {}
            # Calculate the mean of each variable in the dictionary given
            for d in dict_:
                for a in analysis:
                    calc, calc2, abs_diff_cube = None, None, None
                    if a == 'mean':
                        if spatial:
                            calc = dict_[d].collapsed(spatial_names, iris.analysis.MEAN)
                            calc2 = dict2_[d].collapsed(spatial_names, iris.analysis.MEAN)
                            c = calc.copy()
                            mean_sp[d] = c
                            mean_sp2[d] = calc2
                        else:
                            calc = dict_[d].collapsed(time_name, iris.analysis.MEAN)
                            calc2 = dict2_[d].collapsed(time_name, iris.analysis.MEAN)
                        sub_cube = np.subtract(calc.data, calc2.data)
                        abs_diff_cube = np.abs(sub_cube)
                        calc.data = abs_diff_cube
                        mean_calcs[d] = calc
                    if a == 'std':
                        if spatial:
                            calc = dict_[d].collapsed(spatial_names, iris.analysis.STD_DEV)
                            calc2 = dict2_[d].collapsed(spatial_names, iris.analysis.STD_DEV)
                            c = calc.copy()
                            std_sp[d] = c
                            std_sp2[d] = calc2
                        else:
                            calc = dict_[d].collapsed(time_name, iris.analysis.STD_DEV)
                            calc2 = dict2_[d].collapsed(time_name, iris.analysis.STD_DEV)
                        sub_cube = np.subtract(calc.data, calc2.data)
                        abs_diff_cube = np.abs(sub_cube)
                        calc.data = abs_diff_cube
                        std_calcs[d] = calc
                    if a == 'median':
                        if spatial:
                            calc = dict_[d].collapsed(spatial_names, iris.analysis.MEDIAN)
                            calc2 = dict2_[d].collapsed(spatial_names, iris.analysis.MEDIAN)
                            c = calc.copy()
                            median_sp[d] = c
                            median_sp2[d] = calc2
                        else:
                            calc = dict_[d].collapsed(time_name, iris.analysis.MEDIAN)
                            calc2 = dict2_[d].collapsed(time_name, iris.analysis.MEDIAN)
                        sub_cube = np.subtract(calc.data, calc2.data)
                        abs_diff_cube = np.abs(sub_cube)
                        calc.data = abs_diff_cube
                        median_calcs[d] = calc
                    if a == 'rms':
                        if spatial:
                            calc = dict_[d].collapsed(spatial_names, iris.analysis.RMS)
                            calc2 = dict2_[d].collapsed(spatial_names, iris.analysis.RMS)
                            c = calc.copy()
                            rms_sp[d] = c
                            rms_sp2[d] = calc2
                        else:
                            calc = dict_[d].collapsed(time_name, iris.analysis.RMS)
                            calc2 = dict2_[d].collapsed(time_name, iris.analysis.RMS)
                        sub_cube = np.subtract(calc.data, calc2.data)
                        abs_diff_cube = np.abs(sub_cube)
                        calc.data = abs_diff_cube
                        rms_calcs[d] = calc
            calcs, c_sp, c_sp2 = [], [], []
            if 'mean' in analysis:
                calcs.append(mean_calcs)
                if spatial:
                    c_sp.append(mean_sp)
                    c_sp2.append(mean_sp2)
            if 'std' in analysis:
                calcs.append(std_calcs)
                if spatial:
                    c_sp.append(std_sp)
                    c_sp2.append(std_sp2)
            if 'median' in analysis:
                calcs.append(median_calcs)
                if spatial:
                    c_sp.append(median_sp)
                    c_sp2.append(median_sp2)
            if 'rms' in analysis:
                calcs.append(rms_calcs)
                if spatial:
                    c_sp.append(rms_sp)
                    c_sp2.append(rms_sp2)
            ens_calcs.append(calcs)
            if spatial:
                spat_calcs.append(c_sp)
                spat_calcs2.append(c_sp2)

        # Find ensembles average
        if total:
            return self.compute_total_stats_analysis(analysis, nan_indices, spatial, ens_calcs)

        print("function compute_stats_analysis: Averages of data successfully computed.")

        return ens_calcs, spat_calcs, spat_calcs2, analysis, nan_indices



