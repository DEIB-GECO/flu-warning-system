from utils.utils import Table, DateRange
import streamlit as st
import plotly.express as px
import pandas as pd
from datetime import timedelta

def map_warnings_selector(window_sizes, ks, combinations):
    col0, col1, col2 = st.columns([1, 1,1])
    with col0:
        st.markdown("**Configuration of the displayed warnings**")
        st.markdown("<small>Select one value of window size and K to display the related warnings.</small>", unsafe_allow_html=True)
    with col1:
        w_selection = st.pills(
            "window size", 
            window_sizes, 
            selection_mode="single", 
            key="w_selection",
            default=st.session_state.get("w_selection", None))
    with col2:
        previous_value = st.session_state.get("k_selection", None)
        previous_value = previous_value if previous_value in [k for w,k in combinations if w == w_selection] else None    # if the previous_value is not valid, set it to None
        k_selection = st.pills(
            "K", 
            ks if w_selection is None else [k for w,k in combinations if w == w_selection], 
            selection_mode="single", 
            key="k_selection", 
            default=previous_value)
    return w_selection, k_selection


def mapbox_input_plus_warnings_in_daterange(input_table, date_range, warning_table = None, sel_host_type = [], sel_country_or_state=[], sel_location=[]):
    """
    Returns a DataFrame with columns "country_or_state", "lat", "lon", "quantity", "type"
    where "type" can be "all"/"warning"
    """
    map_input_seq_table = (Table.filter(input_table, sel_host_type, sel_country_or_state, sel_location, date_range)[["country_or_state", "lat", "lon"]]
                        .groupby(["country_or_state", "lat", "lon"]).size()
                        .reset_index().rename(columns={0:'quantity'}))
    map_input_seq_table['type'] = 'all'

    if warning_table is not None:
        map_warning_table = (Table.filter(warning_table, sel_host_type, sel_country_or_state, sel_location, date_range)
                             [["country_or_state", "lat", "lon"]]
                             .groupby(["country_or_state", "lat", "lon"]).size()
                             .reset_index().rename(columns={0:'quantity'}))
        map_warning_table['type'] = 'warning'
        map_data = pd.concat([map_input_seq_table, map_warning_table], ignore_index=True).astype({"quantity": int}).sort_values(by="type")
    else:
        map_data = map_input_seq_table

    fig = px.scatter_mapbox(map_data, lat='lat', lon='lon', size ='quantity', 
                            color="type", color_discrete_map={'warning': 'rgba(255, 20, 0, 1)', 'all': 'rgba(0, 0, 0, 1)'}, 
                            hover_name="country_or_state", hover_data = {'type': True, 'quantity': True, 'lat': False, 'lon': False},
                            center={'lat':map_data.lat.mean(), 'lon': map_data.lon.mean()}, 
                            zoom=2.5, mapbox_style='carto-positron', 
                            height=900, width=1200, template='simple_white')
    fig.update(layout_coloraxis_showscale=False)
    fig.update_layout({
        'plot_bgcolor': 'rgba(0, 0, 0, 0)',
        'paper_bgcolor': 'rgba(0, 0, 0, 0)',
        'modebar': {
                   'orientation': 'h',
                   'bgcolor': '#ffffff',
                   'color':'rgb(250, 125, 127)',
                   'activecolor':'rgb(250, 125, 127)'}
               })
    return map_data, fig
    
    
def mapbox_input_plus_warnings_in_daterange_by_week(input_table, date_range, warning_table = None, sel_host_type = [], sel_country_or_state=[], sel_location=[], animation_frame=None):
    """
    Returns a DataFrame with columns "country_or_state", "lat", "lon", "Week", "quantity", "type"
    where "type" can be "all"/"warning"
    """
    if warning_table is not None:
        map_warning_table = (Table.filter(warning_table, sel_host_type, sel_country_or_state, sel_location, date_range)
                             [["country_or_state", "lat", "lon", "Week"]]
                             .groupby(["country_or_state", "lat", "lon", "Week"]).size()
                             .reset_index().rename(columns={0:'quantity'}))
    else:
        map_warning_table = pd.DataFrame([], columns=["country_or_state", "lat", "lon", "Week", 'quantity'])
    map_warning_table['type'] = 'warning'

    map_input_seq_table = (Table.filter(input_table, sel_host_type, sel_country_or_state, sel_location, date_range)[["country_or_state", "lat", "lon", "Week"]]
                        .groupby(["country_or_state", "lat", "lon", "Week"]).size()
                        .reset_index().rename(columns={0:'quantity'}))
    map_input_seq_table['type'] = 'all'

    map_data = pd.concat([map_input_seq_table, map_warning_table], ignore_index=True).astype({"quantity": int}).sort_values(by=["Week","type"])

    fig = px.scatter_mapbox(map_data, lat='lat', lon='lon', hover_name="country_or_state", 
                            color="type", size ='quantity', color_discrete_map={'warning': 'red', 'all': 'black'},
                            hover_data = {'type': True, 'quantity': True, 'lat': False, 'lon': False, 'Week': False},
                            center={'lat':map_data.lat.mean(), 'lon': map_data.lon.mean()}, 
                            zoom=2.5, mapbox_style='carto-positron', 
                            height=900, width=1200, template='simple_white',
                            animation_frame=animation_frame)
    fig.update(layout_coloraxis_showscale=False)
    return map_data, fig


def mapbox_input_plus_warnings_v2(input_sequences_table, date_range, warning_table = None, sel_host_type = [], sel_country_or_state=[], sel_locations=[]):
    # define period for binning sequences
    _, test_ranges = DateRange.week_step_date_ranges_split_train_test_moving_window_v3(
        date_range,
        train_weeks=1,
        test_weeks=2,
        window_shift=2,
    )
    bins = [test_ranges[0][0]] + [d[1] for d in test_ranges]
    period_labels = [f'{b1.strftime("%G-%V")}  {b2.strftime("%G-%V")}' for b1,b2 in zip(bins[:-1], bins[1:])]
    
    # filter warnings and assign period
    if warning_table is not None:
        map_warning_table = (Table.filter(warning_table, sel_host_type, sel_country_or_state, [], date_range)
                             [["country_or_state", "lat", "lon", "Collection_Date"]])
        map_warning_table["Time Period"] = pd.cut(map_warning_table["Collection_Date"], bins, labels=period_labels, right=False, include_lowest=True)
        map_warning_table = map_warning_table.drop(columns=["Collection_Date"])
        # group and count warnings
        map_warning_table = (map_warning_table
                             .groupby(by=["country_or_state", "lat", "lon", "Time Period"])
                             .size()
                             .reset_index().rename(columns={0:'quantity'}))
        map_warning_table = map_warning_table[map_warning_table.quantity > 0]
        map_warning_table['type'] = 'warning'
    else:
        map_warning_table = pd.DataFrame([], columns=["country_or_state", "lat", "lon", "Time Period", 'quantity', "type"], dtype={"country_or_state": str, "lat": float, "lon": float, "Time Period": str, 'quantity': int, "type": str})
    
    # filter input sequences and assign period
    map_input_seq_table = (Table.filter(input_sequences_table, sel_host_type, sel_country_or_state, [], date_range)
                           [["country_or_state", "lat", "lon", "Collection_Date"]])
    map_input_seq_table["Time Period"] = pd.cut(map_input_seq_table["Collection_Date"], bins, labels=period_labels, right=False, include_lowest=True)
    map_input_seq_table = map_input_seq_table.drop(columns=["Collection_Date"])
    # group and count input sequences
    map_input_seq_table = (map_input_seq_table
                            .groupby(by=["country_or_state", "lat", "lon", "Time Period"])
                            .size()
                            .reset_index().rename(columns={0:'quantity'}))
    map_input_seq_table = map_input_seq_table[map_input_seq_table.quantity > 0]
    map_input_seq_table['type'] = 'all'

    # join tables
    map_data = pd.concat([map_input_seq_table, map_warning_table], ignore_index=True).astype({"quantity": int}).sort_values(by=["Time Period","type"])
    st.dataframe(map_data)
    st.write(map_data.shape[0])
    fig = px.scatter_mapbox(map_data, lat='lat', lon='lon', hover_name="country_or_state", 
                            color="type", size ='quantity', color_discrete_map={'warning': 'red', 'all': 'black'},
                            hover_data = {'type': True, 'quantity': True, 'lat': False, 'lon': False, 'Time Period': False},
                            center={'lat':map_data.lat.mean(), 'lon': map_data.lon.mean()}, 
                            zoom=2.5, mapbox_style='carto-positron', 
                            height=900, width=1200, template='simple_white',
                            animation_frame="Time Period")
    fig.update(layout_coloraxis_showscale=False)
    return fig


def time_periods(date_range: DateRange):
    _, test_ranges = DateRange.week_step_date_ranges_split_train_test_moving_window_v3(
        date_range,
        train_weeks=1,
        test_weeks=2,
        window_shift=2,
    )
    bins = [test_ranges[0][0]] + [d[1] for d in test_ranges]
    period_labels = [f"{t1.strftime('%G-%V')}  {(t2-timedelta(days=1)).strftime('%G-%V')}" for t1,t2 in test_ranges]
    return bins, period_labels, test_ranges


def coordinates_and_count_of_input_and_warnings(input_sequences_table, date_range, warning_table = None, sel_host_type = [], sel_country_or_state=[], sel_locations=[]):
    """
    Returns a DataFrame with columns "country_or_state", "lat", "lon", "Time Period", "quantity", "type"
    where "type" can be "all"/"warning"
    """
     # define period for binning sequences
    _, test_ranges = DateRange.week_step_date_ranges_split_train_test_moving_window_v3(
        date_range,
        train_weeks=1,
        test_weeks=2,
        window_shift=2,
    )
    bins = [test_ranges[0][0]] + [d[1] for d in test_ranges]
    period_labels = [f'{b1.strftime("%G-%V")}  {b2.strftime("%G-%V")}' for b1,b2 in zip(bins[:-1], bins[1:])]
    
    # filter warnings and assign period
    if warning_table is not None:
        map_warning_table = (Table.filter(warning_table, sel_host_type, sel_country_or_state, [], date_range)
                             [["country_or_state", "lat", "lon", "Collection_Date"]])
        map_warning_table["Time Period"] = pd.cut(map_warning_table["Collection_Date"], bins, labels=period_labels, right=False, include_lowest=True)
        map_warning_table = map_warning_table.drop(columns=["Collection_Date"])
        # group and count warnings
        map_warning_table = (map_warning_table
                             .groupby(by=["country_or_state", "lat", "lon", "Time Period"])
                             .size()
                             .reset_index().rename(columns={0:'quantity'}))
        map_warning_table = map_warning_table[map_warning_table.quantity > 0]
        map_warning_table['type'] = 'warning'
    else:
        map_warning_table = pd.DataFrame([], columns=["country_or_state", "lat", "lon", "Time Period", 'quantity', "type"], dtype={"country_or_state": str, "lat": float, "lon": float, "Time Period": str, 'quantity': int, "type": str})
    
    # filter input sequences and assign period
    map_input_seq_table = (Table.filter(input_sequences_table, sel_host_type, sel_country_or_state, [], date_range)
                           [["country_or_state", "lat", "lon", "Collection_Date"]])
    map_input_seq_table["Time Period"] = pd.cut(map_input_seq_table["Collection_Date"], bins, labels=period_labels, right=False, include_lowest=True)
    map_input_seq_table = map_input_seq_table.drop(columns=["Collection_Date"])
    # group and count input sequences
    map_input_seq_table = (map_input_seq_table
                            .groupby(by=["country_or_state", "lat", "lon", "Time Period"])
                            .size()
                            .reset_index().rename(columns={0:'quantity'}))
    map_input_seq_table = map_input_seq_table[map_input_seq_table.quantity > 0]
    map_input_seq_table['type'] = 'all'

    # join tables
    return pd.concat([map_input_seq_table, map_warning_table], ignore_index=True).astype({"quantity": int}).sort_values(by=["Time Period","type"])


