import streamlit as st
from streamlit_extras.stylable_container import stylable_container
import numpy as np
import pandas as pd
from streamlit.components.v1 import html
import utils.db_queries as dbq
from utils.db_queries import CachedDB
from utils.utils import DateRange, read_parquet_cached, table_host_types, table_countries_and_states, Table
import seaborn as sns
import plotly.express as px
import os
from datetime import datetime, timedelta
from utils.aggregated_view import *
from utils.individual_sequences_view import tabs_individual_sequences
from utils import filters_view
from utils.map_view import map_warnings_selector, time_periods, mapbox_input_plus_warnings_in_daterange

# Database engine
db_engine = st.connection(url="sqlite:///data/h1n1_08-09.sqlite?mode=ro", type="sql", name="h1n1_08-09")  

# Page structure
window_sizes = (5, 10, 50, 100)
ks = (1, 3, 5, 10, 15)
combinations = [(w, k) for w in window_sizes for k in ks if k < w]
tab_titles = ["Input sequences"] + [f"Window {w} k {k}" for w, k in combinations]
tab2table_name_dict = dict(
    zip(tab_titles, ["input_data"] + [f"window{w}_k{k}" for w, k in combinations])
)

input_sequences_table = CachedDB.get_input_data_table(db_engine)
table_shown_columns = [
    "Lineage",
    "Clade",
    "bmc_cluster_label",
    "Host_Type",
    "Collection_Date",
    "Week",
    "country_or_state",
    "IsolateID_Name"
]



##### FILTERS
st.write("Browse and filter the sequences")
unique_host_types = table_host_types(input_sequences_table)
sel_host_type = filters_view.host_type_multiselect(unique_host_types)

unique_country_or_state = table_countries_and_states(input_sequences_table)
sel_country_or_state = filters_view.locations_multiselect(unique_country_or_state)


##### AGGREGATED VIEW DF
st.markdown("### Aggregated view")
aggregated_view = compute_agg_view_df(
    date_range=DateRange.from_iso_weeks_right_included("2008-31", "2010-01"), 
    input_sequences_table=input_sequences_table, 
    db_engine=db_engine,
    st_tab_name2db_table_name={k: tab2table_name_dict[k] for k in tab_titles if k.startswith("Window")}, 
    combinations=combinations,
    sel_host_type=sel_host_type, sel_country_or_state=sel_country_or_state)

##### AGGREGATED VIEW STYLE
#### cell background color mapping
agg_view_style = aggregated_view.style.set_properties(**{'background-color': 'white'})
aggview_background_controller = BackgroundColor(aggregated_view, agg_view_style, combinations)
aggview_background_controller.enable_aggview_background_color_switch()

#### table style
agg_view_style = (agg_view_style.set_sticky(axis="columns").hide(axis='index'))
htmltable_helper = XMLTableEditor(agg_view_style.to_html())
XMLTableEditor.add_attribute(htmltable_helper.table_root, "style", "width:100%")

st.html(htmltable_helper.format_with_style())



##### TABS
if not "individual_sequences_date_range_h1n1_08-09" in st.session_state:
    st.session_state["individual_sequences_date_range_h1n1_08-09"] = DateRange.from_dates("2008-07-28", "2010-01-10")
tabs_individual_sequences(tab_titles, tab2table_name_dict, db_engine, table_shown_columns, st.session_state["individual_sequences_date_range_h1n1_08-09"], 
                          sel_host_type=sel_host_type, sel_country_or_state=sel_country_or_state)



##### MAP
st.markdown("### Map of warnings")
sel_map_warning_config = map_warnings_selector(window_sizes, ks, combinations)
display_warnings = all(sel_map_warning_config)

if display_warnings:
    warning_table = CachedDB.get_table_named_as(db_engine, f"window{sel_map_warning_config[0]}_k{sel_map_warning_config[1]}")
else:
    warning_table = None

time_bins, period_labels, test_ranges = time_periods(DateRange.from_iso_weeks_right_included("2008-31", "2010-01"))
period_label2range = dict(zip(period_labels, test_ranges))
col0, col1 = st.columns([1,2])
with col0:
    st.markdown("**Time Period**")
with col1:
    default_value = st.session_state.get("map_view_time_period_slider", period_labels[0])
    st.select_slider("whatever", options=period_labels, value=default_value, key="map_view_time_period_slider", 
                     help="Move the slider to the desired time period", label_visibility="hidden")

selected_time_period = period_label2range[st.session_state['map_view_time_period_slider']]
selected_time_period = DateRange.from_dates(*period_label2range[st.session_state['map_view_time_period_slider']])

map_data, fig = mapbox_input_plus_warnings_in_daterange(
    CachedDB.get_table_named_as(db_engine, "input_data"),
    selected_time_period,
    warning_table, 
    sel_host_type, sel_country_or_state)
number_all_sequences = map_data[map_data.type=="all"].quantity.sum()
number_warnings = map_data[map_data.type=="warning"].quantity.sum()
number_unique_locations = len(map_data.country_or_state.unique())

col0, col1 = st.columns([1, 2])
with col0:
    st.markdown("Currently displaying:")
with col1:
    st.markdown(f"{number_unique_locations} unique locations")   
    if display_warnings: 
        st.markdown(f"{number_warnings} warnings in {number_all_sequences} sequences")
    else:
        st.markdown(f"{number_all_sequences} sequences")

st.plotly_chart(fig, use_container_width=True, on_select="ignore")
st.markdown("""
            <style>            
            .stPlotlyChart { margin-top: -60px;    }
            </style>
            """, unsafe_allow_html=True)
st.markdown("**Summary**")
resume = map_data.rename(columns={'country_or_state':'locations'}).pivot(index='locations', columns='type', values='quantity').fillna(0).rename(columns={'all':'sequences', 'warning': 'warnings'})
if 'warnings' not in resume.columns:
    resume['warnings'] = 0
st.dataframe(resume)
