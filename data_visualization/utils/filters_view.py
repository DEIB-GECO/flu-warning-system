import streamlit as st
from utils.utils import DateRange
import pandas as pd
from datetime import datetime, timedelta

def host_type_multiselect(unique_host_types) -> list:
    return st.multiselect(
        label="Host Type",
        options=unique_host_types,
        help="Select one or more values to enable the filter",
        placeholder="Select one or more values to enable the filter",
    )

def locations_multiselect(unique_countries_and_continents) -> list:
    return st.multiselect(
        label="Location",
        options=unique_countries_and_continents,
        help="Select one or more values to enable the filter",
        placeholder="Select one or more values to enable the filter",
    )

def collection_date_slider(intial_range: tuple[datetime], default_range: list=None):
    default_range = default_range if default_range is not None else intial_range
    sel_collection_date_range = st.slider(
        label="Collection Date", 
        min_value=intial_range[0], max_value=intial_range[1], value=default_range, 
        step=timedelta(days=14), format="YYYY MMM D", 
        key="coldate-selector", help="Move the right/left slider to enable the filter", 
        on_change=None, args=None, kwargs=None
    )
    date_range = DateRange()
    date_range.start = pd.Timestamp(sel_collection_date_range[0])
    date_range.end = pd.Timestamp(sel_collection_date_range[1])
    return date_range

def genotype_multiselect(unique_genotypes) -> list:
    return st.multiselect(
        label="Genotypes",
        options=unique_genotypes,
        help="Select one or more values to enable the filter",
        placeholder="Select one or more values to enable the filter",
    )