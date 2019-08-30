"""
Created by Adanna Akwataghibe (Github: AdannaAkwats)
"""
import sys

sys.path.append(".")

from utils import *
import random
import pytest


def test_valid_order():
    start = [21, 5, 1988]
    end = [19, 6, 1990]
    assert check_valid_order(start, end)


def test_invalid_order():
    start = [21, 5, 2020]
    end = [19, 6, 1990]
    assert not check_valid_order(start, end)


def test_invalid_order_2():
    start = [20, 6, 1990]
    end = [19, 6, 1990]
    assert not check_valid_order(start, end)


def test_check_analysis_valid():
    as_ = ['mean', 'rms', 'std', 'median']
    a = [random.choice(as_)]
    res = check_analysis(a)
    assert res == True


def test_check_analysis_not_valid():
    a = ['not_valid']
    with pytest.raises(SystemExit):
        check_analysis(a)


def test_check_analysis_not_list():
    a = 'not-a-list'
    with pytest.raises(AssertionError):
        check_analysis(a)


def test_check_variables_true():
    varbs = ['var1', 'var2']
    res = check_variables_covary(varbs)
    assert res == True


def test_check_variables_system_exit():
    varbs = ['var1']
    with pytest.raises(SystemExit):
        check_variables_covary(varbs)


def test_check_variables_assertion_error():
    varbs = 'not-a-list'
    with pytest.raises(AssertionError):
        check_variables_covary(varbs)


def test_no_overlap():
    x1, x2 = 2005, 2008
    y1, y2 = 2003, 2004
    assert overlaps(x1, x2, y1, y2) == False


def test_overlap():
    x1, x2 = 2005, 2008
    y1, y2 = 2003, 2007
    assert overlaps(x1, x2, y1, y2) == True


def test_get_ens_num_match():
    f = 'ex_ens101_2000.nc'
    match = get_ens_num(f)
    assert match == 101


def test_get_ens_num_no_ens():
    f = 'ex_101_2000.nc'
    match = get_ens_num(f)
    assert match == 101


def test_get_file_two_years_match():
    f = 'ex_ens101_2000_2005.nc'
    match1, match2 = get_file_two_years(f)
    assert match1 == 2000 and match2 == 2005


def test_get_file_two_years_not_match():
    f = 'ex_ens101_2000.nc'
    match1 = get_file_two_years(f)
    assert match1 == False


def test_ens_to_indx_valid():
    indx = ens_to_indx(1001, 2)
    assert indx == 0


def test_ens_to_indx_larger_than_ensembles():
    indx = ens_to_indx(1003, 2)
    assert indx == -1


def test_ens_to_indx_not_valid():
    indx = ens_to_indx(10, 2)
    assert indx == -1


def test_ens_to_indx_system_exit():
    with pytest.raises(SystemExit):
        ens_to_indx(1000, 2, max_start=100)


def test_find_middle():
    arr = [1,2,3]
    mid, mid_indx = find_middle(arr)
    assert mid == 2 and mid_indx == 1


def test_find_middle_even_list():
    arr = [1,2]
    mid, mid_indx = find_middle(arr)
    assert mid == 1 and mid_indx == 0


def test_find_middle_one_elem_list():
    arr = [1]
    mid, mid_indx = find_middle(arr)
    assert mid == 1 and mid_indx == 0



def test_nested_list():
    l = [[1]]
    assert is_nested_list(l)


def test_not_nested_list():
    l = [1]
    assert not is_nested_list(l)


def test_check_list_date_valid():
    date = [2000, 3, 4]
    assert check_list_date(date)


def test_check_list_date_valid():
    date = [2000]
    assert not check_list_date(date)


def test_make_into_file_name():
    str = 'file 1'
    res = make_into_file_name(str)
    assert res == 'file_1'


def test_make_into_file_name_no_change():
    str = 'file'
    res = make_into_file_name(str)
    assert res == 'file'
