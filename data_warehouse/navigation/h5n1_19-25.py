import streamlit as st
from streamlit_extras.stylable_container import stylable_container
import numpy as np
import pandas as pd
from streamlit.components.v1 import html
import utils.db_queries as dbq
from utils.db_queries import CachedDB
from utils.utils import DateRange, read_parquet_cached, table_host_types, table_countries_and_continents
import seaborn as sns
from datetime import datetime, timedelta
import os
from utils.aggregated_view import *
from utils.individual_sequences_view import tabs_individual_sequences
from utils import filters_view

# Database engine
db_engine = st.connection(url="sqlite:///data/h5n1_2019-2025-01-27.sqlite?mode=ro", type="sql", name="h5n1_19-25")  

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
    "Pathogenicity",
    "Host_Type",
    "Collection_Date",
    "Week",
    "Location",
    "IsolateID_Name"
    # "country_or_state",
]



##### FILTERS
st.write("Browse and filter the sequences")
unique_host_types = table_host_types(input_sequences_table)
sel_host_type = filters_view.host_type_multiselect(unique_host_types)

unique_countries_and_continents = table_countries_and_continents(input_sequences_table)
sel_locations = filters_view.locations_multiselect(unique_countries_and_continents)

intial_range = datetime(2019, 1, 7), datetime(2025, 1, 27, 23, 59, 59, 999999) 
date_range = filters_view.collection_date_slider(intial_range)



##### AGGREGATED VIEW DF
st.markdown("### Aggregated view")
aggregated_view = compute_agg_view_df(
    date_range=date_range, 
    input_sequences_table=input_sequences_table, 
    db_engine=db_engine,
    st_tab_name2db_table_name={k: tab2table_name_dict[k] for k in tab_titles if k.startswith("Window")}, 
    combinations=combinations,
    sel_host_type=sel_host_type, sel_locations=sel_locations)

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
                          sel_host_type=sel_host_type, sel_locations=sel_locations)
