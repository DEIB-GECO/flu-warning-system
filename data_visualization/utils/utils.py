import pandas as pd
import numpy as np
from datetime import datetime, timedelta
import streamlit as st
import toml

class DateRange:
    def __init__(self):
        self.start = None
        self.end = None

    def __str__(self):
        return f"{self.start}_{self.end}"

    @staticmethod
    def from_dates(start: str = '2000-01-01', end: str = 'today'):
        """
        The date range is defined as [start,end).
        :param str start: the first date included in the range
        :param str end: the closing date of the range, that should not be included
        """
        obj = DateRange()
        obj.start = pd.to_datetime(start)
        obj.end = pd.to_datetime(end)
        return obj

    @staticmethod
    def from_iso_weeks(start: str, end: str):
        start += "-1" if not start.endswith("-1") else ""
        end += "-1" if not end.endswith("-1") else ""
        obj = DateRange()
        obj.start = pd.to_datetime(start, format="%G-%V-%u") #2023-46-1: Monday of 46th week of 2023 according to ISO8601
        obj.end = pd.to_datetime(end, format="%G-%V-%u")     
        return obj
    
    @staticmethod
    def from_iso_weeks_right_included(start: str, end: str):
        start += "-1" if not start.endswith("-1") else ""
        end += "-7" if not end.endswith("-1") else ""
        obj = DateRange()
        obj.start = pd.to_datetime(start, format="%G-%V-%u") #2023-46-1: Monday of 46th week of 2023 according to ISO8601
        obj.end = pd.to_datetime(end, format="%G-%V-%u") + pd.Timedelta(days=1)
        return obj

    @staticmethod
    def week_step_date_ranges(last_datetime: pd.Timestamp = pd.to_datetime('today'), start_datetime: pd.Timestamp = pd.to_datetime('2000-01-01'), weeks=16):
        """
        Returns pairs of <start,end> pd.Timestamp spanning the given number of weeks. Edges of the range should be considered as left-included, right-excluded. 
        The argument last_datetime is always included in the last range. The argument start_datetime is always included in the first range and may fall well after the 
        left-edge of the first range. 
        :param pd.Timestamp last_datetime: last date (most recent) that must be included in the range. This date is included by the last date range returned and the 
        start date of the range is positioned the argment number of weeks backward.
        :param pd.Timestamp start_datetime: date included in the oldest date range. 
        :param int weeks: the number of weeks defining the date range width.
        """
        end_interval = last_datetime.round(freq="D") + pd.Timedelta(days=1) # is excluded (round(freq="D") sets to 00:00, +1d includes the given date)
        start_interval = end_interval
        output = []
        while start_interval > start_datetime:
            start_interval = end_interval - pd.Timedelta(weeks=weeks)
            output.append((start_interval, end_interval))
            end_interval = start_interval
        return output[::-1]
    
    @staticmethod
    def week_step_date_ranges_split_train_test_moving_window(last_date: pd.Timestamp = pd.to_datetime('today'), start_date: pd.Timestamp = pd.to_datetime('2000-01-01'), train_weeks=16, test_weeks=4, window_weeks=4):
        """
        Returns paired time intervals. The first set of intervals is the training set. The second set is the test set. 
        The i-th interval in the test set always follows the i-th interval in the training set. 
        An interval is expressed as a tuple of <start,end> pd.Timestamp objects.
        The i-th+1 interval starts after the number of weeks specified in window_weeks. If this number is equal to training_weeks+test_weeks, 
        the i-th+1 interval does not overlap with the i-th interval. 
        :param pd.Timestamp last_datetime: last date (most recent). The produced ranges will stop before this date
        :param pd.Timestamp start_datetime: date included in the oldest date range. 
        :param int train_weeks: the number of weeks defining the date range width.
        :param int test_weeks: the number of weeks defining the date range width.
        :param int window_weeks: the number of weeks defining the date range width.
        """
        output_train, output_test = [], []
        window_start = start_date
        while(window_start + pd.Timedelta(weeks=train_weeks+test_weeks)) < last_date:
            train_end = window_start + pd.Timedelta(weeks=train_weeks)
            output_train.append((window_start, train_end))
            test_start = train_end
            test_end = test_start + pd.Timedelta(weeks=test_weeks)
            output_test.append((test_start, test_end))
            window_start = window_start + pd.Timedelta(weeks=window_weeks)
        return output_train, output_test
    
    @staticmethod
    def week_step_date_ranges_split_train_test_moving_window_v2(date_range, train_weeks=16, test_weeks=4, window_shift=4):
        """
        As week_step_date_ranges_split_train_test_moving_window, but last_date and start_date are referred to the test window.
        Receives a date range of type DateRange
        The first window is the week containing the start date, the last window the last possible window that do not iclude the end date of the range
        """
        output_train, output_test = [], []
        start_date, end_date = date_range.start, date_range.end

        # ensure start_date begins with a monday and last_date ends with monday (excluded)
        if start_date.isoweekday() > 1:
            start_date -= pd.Timedelta(days=start_date.isoweekday() - 1)
        if end_date.isoweekday() > 1:
            end_date -= pd.Timedelta(days=end_date.isoweekday() - 1)

        delta_train = pd.Timedelta(weeks=train_weeks)
        delta_test = pd.Timedelta(weeks=test_weeks)
        delta_shift = pd.Timedelta(weeks=window_shift)

        window_start = start_date
        while window_start + delta_test < end_date:
            # training windows
            train_start = window_start - delta_train
            train_end = window_start
            # test windows
            test_start = window_start
            test_end = test_start + delta_test
            # update rule
            output_train.append((train_start, train_end))
            output_test.append((test_start, test_end))
            window_start += delta_shift
        return output_train, output_test
    
    @staticmethod
    def week_step_date_ranges_split_train_test_moving_window_v3(date_range, train_weeks=16, test_weeks=4, window_shift=4):
        """
        As week_step_date_ranges_split_train_test_moving_window, but last_date and start_date are referred to the test window.
        Receives a date range of type DateRange
        The first window is the week containing the start date, the last window the last possible window that do not iclude the end date of the range
        """
        output_train, output_test = [], []
        start_date, end_date = date_range.start, date_range.end

        # ensure start_date begins with a monday and last_date ends with monday (excluded)
        if start_date.isoweekday() > 1:
            start_date -= pd.Timedelta(days=start_date.isoweekday() - 1)
        if end_date.isoweekday() > 1:
            end_date -= pd.Timedelta(days=end_date.isoweekday() - 1)

        delta_train = pd.Timedelta(weeks=train_weeks)
        delta_test = pd.Timedelta(weeks=test_weeks)
        delta_shift = pd.Timedelta(weeks=window_shift)

        window_start = start_date
        while window_start + delta_test <= end_date:
            # training windows
            train_start = window_start - delta_train
            train_end = window_start
            # test windows
            test_start = window_start
            test_end = test_start + delta_test
            # update rule
            output_train.append((train_start, train_end))
            output_test.append((test_start, test_end))
            window_start += delta_shift
        return output_train, output_test
    

class ReleaseNotes:

    file_path = "./release_notes.md"

    @staticmethod
    @st.dialog("Release Notes")
    def print_release_notes():
        with open(ReleaseNotes.file_path, "r") as f:
            st.markdown(f.read())

    @staticmethod
    @st.cache_resource
    def current_version():
        with open(ReleaseNotes.file_path, "r") as f:
            rnotes = f.read()
            return rnotes.lstrip("# v").split(" ", maxsplit=1)[0]
        

class StreamlitConfig:

    @staticmethod
    def hide_streamlit_top_bar():
        st.markdown("""
        <style>
        MainMenu {visibility: hidden;}
        footer {visibility: hidden;}
        .stAppHeader {visibility: hidden; hight: 0px; background: rgba(0,0,0,0); !important;} 
        .stAppToolbar {visibility: hidden; hight: 0px; background: rgba(0,0,0,0); !important;} 
        root > div:nth-child(1) > div > div > div > div > section > div {padding-top: 0rem;}
        .stMainBlockContainer {padding-top: 1rem; !important;}
        </style>
        """, unsafe_allow_html=True)

    @staticmethod
    def make_streamlit_top_bar_translucent():
        st.markdown("""
        <style>
        .stAppHeader {background: rgba(0,0,0,0.05); !important;} 
        .stAppToolbar {background: rgba(0,0,0,0.05); !important;} 
        </style>
        """, unsafe_allow_html=True)

    @staticmethod
    @st.cache_resource
    def get_config():
        with open(".streamlit/config.toml", "r") as f:
            return toml.load(f)

    @staticmethod
    def autohide_streamlit_top_bar():
        streamlit_config = StreamlitConfig.get_config()
        if streamlit_config.get('ui', False) and streamlit_config['ui'].get('hideTopBar', False):
            StreamlitConfig.hide_streamlit_top_bar()
        else:
            StreamlitConfig.make_streamlit_top_bar_translucent()

@st.cache_resource
def read_parquet_cached(path):
    return pd.read_parquet(path)

@st.cache_resource
def table_host_types(table: pd.DataFrame):
    return table["Host_Type"].unique().tolist()

@st.cache_resource
def table_countries_and_continents(table: pd.DataFrame):
    location_cont_and_countries = pd.Series(table.Location.unique()).str.strip().str.split(" / ", n=2).str[:2]
    unique_continents = location_cont_and_countries.str[0].str.strip().str.capitalize().unique().tolist()
    unique_countries = location_cont_and_countries.str[1].str.strip().str.capitalize().unique().tolist()
    return unique_continents + unique_countries

@st.cache_resource
def table_countries_and_states(table: pd.DataFrame):
    return table.country_or_state.sort_values().unique().tolist()

@st.cache_resource
def table_genotypes(table: pd.DataFrame):
    return pd.Series(table["Genotype"].unique()).sort_values().tolist()

class Table:

    @staticmethod
    def filter(table: pd.DataFrame, sel_host_type: list = [], sel_country_or_state: list = [], sel_locations: list = [], date_range: DateRange = None, sel_genotypes: list = []):
        mask = np.ones(table.shape[0], dtype=bool)  # incrementally add elements to the mask (AND between different filter categories, OR within values of the same filter)
        if sel_host_type:
            mask = table.Host_Type.isin(sel_host_type)  # stands for mask = mask & ...  but the ".isin" does the OR of the argument values 
        if sel_country_or_state:
            mask = mask & (table.country_or_state.isin(sel_country_or_state))
        if sel_locations:
            location_mask = np.zeros(table.shape[0], dtype=bool)     # OR mask
            for sel_location in sel_locations:
                location_mask = location_mask | (table.Location.str.contains(sel_location, regex=False, case=False))
            mask = mask & location_mask
        if date_range:        
            mask = mask & (table.Collection_Date >= date_range.start.to_numpy()) & (table.Collection_Date < date_range.end.to_numpy())
        if sel_genotypes: 
            mask = mask & (table['Genotype'].isin(sel_genotypes))
        return table[mask]


def read_stray_configurations(file_path):
    try:
        with open(file_path) as f:
            configs = [line.replace(',', ' ').split() for line in f if not line.startswith("#")]
            configs = [tuple(map(int,pair)) for pair in configs]
            assert all([len(pairs) == 2 for pairs in configs])
            assert len(configs) > 0
    except:
        raise ValueError(f"Malformed file {file_path}. Make sure the file contains at least one configuration and that each line is formatted as 'number,number'. Lines starting with # are ignored.")
    return configs