import streamlit as st
import numpy as np
import pandas as pd
import utils.db_queries as dbq
from utils.db_queries import CachedDB
from utils.utils import DateRange, table_host_types, table_genotypes, table_countries_and_states, Table
from datetime import datetime
from utils.aggregated_view import *
from utils.individual_sequences_view import tabs_individual_sequences
from utils import filters_view
from utils.map_view import map_warnings_selector, time_periods, mapbox_input_plus_warnings, mapbox_fig

# Database engine
db_engine = st.connection(url="sqlite:///data/h5n1.sqlite?mode=ro", type="sql", name="h1n1")  

# Page structure
# window_sizes = (5, 10, 50, 100)
# ks = (1, 3, 5, 10, 15)
window_sizes = (50,)
ks = (10,)
combinations = [(w, k) for w in window_sizes for k in ks if k < w]
tab_titles = ["Input sequences"] + [f"Window {w} k {k}" for w, k in combinations]
tab2table_name_dict = dict(
    zip(tab_titles, ["input_data"] + [f"window{w}_k{k}" for w, k in combinations])
)

input_sequences_table = CachedDB.get_input_data_table(db_engine)
table_shown_columns = [
    "Lineage",
    "Clade",
    "Genotype",
    "Host_Type",
    "Collection_Date",
    "Week",
    "country_or_state",
    "Location",
    "IsolateID_Name"
]



##### FILTERS
st.write("Browse and filter the sequences")
unique_host_types = table_host_types(input_sequences_table)
sel_host_type = filters_view.host_type_multiselect(unique_host_types)

unique_genotypes = table_genotypes(input_sequences_table)
sel_genotypes = filters_view.genotype_multiselect(unique_genotypes)

unique_country_or_state = table_countries_and_states(input_sequences_table)
sel_country_or_state = st.multiselect(
        label="North American States or USA States",
        options=unique_country_or_state,
        help="Select one or more values to enable the filter",
        placeholder="Select one or more values to enable the filter",
    )

intial_range = datetime(2019, 1, 7), datetime(2025, 5, 5, 23, 59, 59, 999999) 
date_range = filters_view.collection_date_slider(intial_range)

filtered_input_sequences_table = Table.filter(input_sequences_table, sel_host_type, sel_country_or_state, date_range=date_range, sel_genotypes=sel_genotypes)
st.write(f"Selected isolates (inputs) {filtered_input_sequences_table.shape[0]}")

##### AGGREGATED VIEW DF
st.markdown("### Aggregated view")
aggregated_view = compute_agg_view_df(
    date_range=date_range, 
    input_sequences_table=input_sequences_table, 
    db_engine=db_engine,
    st_tab_name2db_table_name={k: tab2table_name_dict[k] for k in tab_titles if k.startswith("Window")}, 
    combinations=combinations,
    sel_host_type=sel_host_type, sel_country_or_state=sel_country_or_state, sel_genotypes=sel_genotypes)

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
tabs_individual_sequences(tab_titles, tab2table_name_dict, db_engine, table_shown_columns, date_range, 
                          sel_host_type=sel_host_type, sel_country_or_state=sel_country_or_state, sel_genotypes=sel_genotypes)

##### MAP
@st.cache_resource
def coordinates_table():
    return pd.read_csv("data/coordinates.csv").replace("", np.nan).replace(" ", np.nan).astype({"lat": pd.Float64Dtype(), "lon": pd.Float64Dtype()})

st.markdown("### Map of warnings")
sel_map_warning_config = map_warnings_selector(window_sizes, ks, combinations)
display_warnings = all(sel_map_warning_config)

if display_warnings:
    warning_table = CachedDB.get_table_named_as(db_engine, f"window{sel_map_warning_config[0]}_k{sel_map_warning_config[1]}")
else:
    warning_table = None

time_bins, period_labels, test_ranges = time_periods(date_range)    # <- reflects user choice filter 
max_period_label = max(period_labels)
min_period_label = min(period_labels)
period_label2range = dict(zip(period_labels, test_ranges))
col0, col1 = st.columns([1,2])
with col0:
    st.markdown("**Time Period**")
with col1:
    default_value = st.session_state.get("map_view_time_period_slider", period_labels[0])
    if default_value > max_period_label:
        default_value = max_period_label
    elif default_value < min_period_label:
        default_value = min_period_label
    st.select_slider("whatever", options=period_labels, value=default_value, key="map_view_time_period_slider", 
                     help="Move the slider to the desired time period", label_visibility="hidden")

selected_time_period = DateRange.from_dates(*period_label2range[st.session_state['map_view_time_period_slider']])

# attach coordinates
map_filtered_input_sequences_table = pd.merge(
    Table.filter(filtered_input_sequences_table, date_range=selected_time_period), 
    coordinates_table(), how="left", on="country_or_state")
if warning_table is None:
    map_filtered_warning_table = None
else:
    map_filtered_warning_table = pd.merge(
        Table.filter(warning_table, sel_host_type, sel_country_or_state, date_range=selected_time_period, sel_genotypes=sel_genotypes), 
        coordinates_table(), how="left", on="country_or_state")

map_data = mapbox_input_plus_warnings(map_filtered_input_sequences_table, map_filtered_warning_table)

# before plotting, remove rows where country_or_state is "USA" 
figure_map_data = map_data[(~map_data.country_or_state.isna()) & (map_data.country_or_state != "USA")]

## Currently displaying:
if figure_map_data.empty:
    st.markdown("<center>No sequences in the selected time frame.\n\nMove the slider to select another time frame or relax the filters at the top of the page.<center>", unsafe_allow_html=True)
else:
    number_all_sequences = figure_map_data[figure_map_data.type=="all"].quantity.sum()
    number_warnings = figure_map_data[figure_map_data.type=="warning"].quantity.sum()
    number_unique_locations = len(figure_map_data.country_or_state.unique())
    col0, col1 = st.columns([1, 2])
    with col0:
        st.markdown("Currently displaying:")
    with col1:
        st.markdown(f"{number_unique_locations} unique locations")   
        if display_warnings: 
            st.markdown(f"{number_warnings} warnings in {number_all_sequences} sequences")
        else:
            st.markdown(f"{number_all_sequences} sequences")

    # append info for isolates in generic "USA"
    not_displayed = (map_data[(map_data.country_or_state.isna()) | (map_data.country_or_state == "USA")])
    if not not_displayed.empty:
        number_all_sequences = not_displayed[not_displayed.type=="all"].quantity.sum()
        number_warnings = not_displayed[not_displayed.type=="warning"].quantity.sum()
        number_unique_locations = len(not_displayed.country_or_state.unique())
        col0, col1 = st.columns([1, 2])
        with col0:
            st.markdown("Sequences not geolocated that could not be displayed:")
        with col1:
            if display_warnings: 
                st.markdown(f"{number_warnings} warnings in {number_all_sequences} sequences")
            else:
                st.markdown(f"{number_all_sequences} sequences")


    fig = mapbox_fig(figure_map_data)
    st.plotly_chart(fig, use_container_width=True, on_select="ignore")
    st.markdown("""
                <style>            
                .stPlotlyChart { margin-top: -60px;    }
                </style>
                """, unsafe_allow_html=True)
    st.markdown(f"**Number of sequences/anomalies by location in the period {st.session_state['map_view_time_period_slider']}**")
    resume = map_data.rename(columns={'country_or_state':'locations'}).pivot(index='locations', columns='type', values='quantity').fillna(0).rename(columns={'all':'sequences', 'warning': 'warnings'})
    if 'warnings' not in resume.columns:
        resume['warnings'] = 0
    st.dataframe(resume)
