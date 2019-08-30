"""
Created by Adanna Akwataghibe (Github: AdannaAkwats)
"""
import sys

sys.path.append(".")

from user_entry import get_date
from Months import Month

"""
Pyrest script for testing functions used for user input
"""


def test_get_date_start_given_year():
    test_start = "2000"
    day, month, year = get_date(test_start)
    assert day == 1
    assert month == Month.January
    assert year == 2000


def test_get_date_end_given_year():
    test_end = "2000"
    day, month, year = get_date(test_end, start=False)
    assert day == 31
    assert month == Month.December
    assert year == 2000


def test_get_date_start_given_year_month():
    test_date = "2000-11"
    day, month, year = get_date(test_date)
    assert day == 1
    assert month == Month.November
    assert year == 2000


def test_get_date_end_given_year_month():
    test_date = "2000-11"
    day, month, year = get_date(test_date, start=False)
    assert day == 30
    assert month == Month.November
    assert year == 2000


def test_get_date_given_year_month_day():
    test_date = "2000-11-05"
    day, month, year = get_date(test_date)
    assert day == 5
    assert month == Month.November
    assert year == 2000


def test_date_not_valid():
    test_date = "2010-10-43"
    assert not get_date(test_date)

