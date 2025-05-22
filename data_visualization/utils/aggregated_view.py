from lxml import etree
from lxml import html
import re
from utils.utils import DateRange, Table
import numpy as np
import pandas as pd
from utils.db_queries import CachedDB
from streamlit_extras.stylable_container import stylable_container
import streamlit as st
import seaborn as sns
from decimal import Decimal
import os


class XMLTableEditor:
    def __init__(self, table_html: str, style_html=""):
        # separate style part containing CSS only from table part containing HTML only
        if "<style" in table_html:
            self.style_root = table_html[: table_html.index("<table")]
            self.table_root = table_html[table_html.index("<table") :]
        else:
            self.style_root = """<style> 
            </style>"""
            self.table_root = table_html
        self.table_root = html.fragment_fromstring(self.table_root)

    class CSS:
        top_border_class = "my-top-border"
        right_border_class = "my-right-border"
        bottom_border_class = "my-bottom-border"
        left_border_class = "my-left-border"

        additional_css = f"""
        table {{
                font-family: arial; 
            }}
        /*table tr td,th {{ 
            width:100% 
        }}*/
        td {{
            padding: 5px; 
        }}
        th {{
            background: #3f577c !important;
            color: white; 
            border:3px solid white; 
            text-align:left !important;
        }}
        .{top_border_class} {{
            border-top: 4px solid black;
            border-collapse: collapse;
        }}
        .{right_border_class} {{
            border-right: 4px solid black;
            border-collapse: collapse;
        }}
        .{bottom_border_class} {{
            border-bottom: 4px solid  black;
            border-collapse: collapse;
        }}
        .{left_border_class} {{
            border-left: 4px solid  black;
            border-collapse: collapse;
        }}     
        """
        re.sub("\n\b+", "\n", additional_css)

    def prettyprint_with_style(self, **kwargs):
        style = self.style_root[: self.style_root.index("</style>")]
        style += XMLTableEditor.CSS.additional_css + "</style>"
        print(style)
        xml = etree.tostring(self.table_root, pretty_print=True, **kwargs)
        print(xml.decode(), end="")

    def format_with_style(self, **kwargs):
        style = self.style_root[: self.style_root.index("</style>")]
        style += XMLTableEditor.CSS.additional_css + "</style>"
        xml = etree.tostring(self.table_root, pretty_print=True, **kwargs)
        return style + "\n" + xml.decode()

    @staticmethod
    def prettyprint(element, **kwargs):
        xml = etree.tostring(element, pretty_print=True, **kwargs)
        print(xml.decode(), end="")

    def cells_at_row_between(self, row, left, right):
        n = self.table_root.xpath(f"tbody/tr[{row + 1}]")[0]
        return [n.xpath(f"td[{i + 1}]")[0] for i in range(left, right + 1)]

    def cells_at_column_between(self, col, top, bottom):
        return [
            self.table_root.xpath(f"tbody/tr[{i + 1}]/td[{col + 1}]")[0]
            for i in range(top, bottom + 1)
        ]

    @staticmethod
    def add_attribute(e, att_name, att_value):
        e.attrib[att_name] = e.attrib.get(att_name, "") + " " + att_value

    def border_to_group_of_cells(self, top, right, bottom, left):
        for c in self.cells_at_row_between(top, left, right):
            XMLTableEditor.add_attribute(
                c, "class", XMLTableEditor.CSS.top_border_class
            )
        for c in self.cells_at_column_between(right, top, bottom):
            XMLTableEditor.add_attribute(
                c, "class", XMLTableEditor.CSS.right_border_class
            )
        for c in self.cells_at_row_between(bottom, left, right):
            XMLTableEditor.add_attribute(
                c, "class", XMLTableEditor.CSS.bottom_border_class
            )
        for c in self.cells_at_column_between(left, top, bottom):
            XMLTableEditor.add_attribute(
                c, "class", XMLTableEditor.CSS.left_border_class
            )
    
    def side_border_to_group_of_cells(self, top, right, bottom, left):
        for c in self.cells_at_column_between(right, top, bottom):
            XMLTableEditor.add_attribute(
                c, "class", XMLTableEditor.CSS.right_border_class
            )
        for c in self.cells_at_column_between(left, top, bottom):
            XMLTableEditor.add_attribute(
                c, "class", XMLTableEditor.CSS.left_border_class
            )


def compute_agg_view_df(date_range, input_sequences_table, 
                        db_engine, st_tab_name2db_table_name, combinations,
                        sel_host_type=[], sel_locations=[], sel_country_or_state=[], sel_genotypes=[]):
    _, test_ranges = DateRange.week_step_date_ranges_split_train_test_moving_window_v3(
        date_range,
        train_weeks=1,
        test_weeks=2,
        window_shift=2,
    )
    bins = [test_ranges[0][0]] + [d[1] for d in test_ranges]

    # count input sequences by time
    filtered_input_sequences = Table.filter(input_sequences_table, sel_host_type, sel_country_or_state, sel_locations, sel_genotypes=sel_genotypes)
    input_count_by_time = filtered_input_sequences.groupby(
        pd.cut(
            filtered_input_sequences.Collection_Date, bins, right=False, include_lowest=True
        ),
        observed=False,
    ).size()

    # count warnings by time for all algorithm variations
    warnings_by_alg = dict()  # alg_name: [] <- warnings_by_time
    for st_tab_name, table_name in st_tab_name2db_table_name.items():
        if not st_tab_name.startswith("Window"):
            continue
        table = CachedDB.get_table_named_as(db_engine, table_name)
        filtered_table = Table.filter(table, sel_host_type, sel_country_or_state, sel_locations, sel_genotypes=sel_genotypes)
        warnings_count_by_time = filtered_table.groupby(
            pd.cut(filtered_table.Collection_Date, bins, right=False, include_lowest=True),
            observed=False,
        ).size()
        warnings_by_alg[st_tab_name] = warnings_count_by_time


    # populate aggregated view
    aggregated_view = pd.DataFrame(
        [],
        index=pd.Series([d[0] for d in test_ranges]),
        columns=pd.MultiIndex.from_tuples(
            combinations,
            names=("window", "k"),
        ),
    )
    aggregated_view["Week"] = aggregated_view.index.strftime("%G-%V").astype(str)
    aggregated_view["N. Sequences"] = input_count_by_time
    for tabname, warnins_by_time in warnings_by_alg.items():
        _, w, _, k = tabname.split(" ")
        w, k = int(w), int(k)
        aggregated_view.loc[:, (w, k)] = warnins_by_time


    # display options 
    aggregated_view = aggregated_view.reset_index().rename(columns={"index": "1st day of the Week"})
    aggregated_view = aggregated_view.loc[:,[("Week", ""), ("1st day of the Week", ""), ("N. Sequences", "")] + combinations]
    aggregated_view["1st day of the Week"] = aggregated_view["1st day of the Week"].dt.strftime("%Y %b %d")

    return aggregated_view


def compute_warnings_ratio_df(aggregated_view: pd.DataFrame, stray_combinations: list[tuple]):
    warn_prop = aggregated_view.copy()
    for c in stray_combinations:
        warn_prop[c] = warn_prop[c].astype("Float64").div(warn_prop[("N. Sequences", "")].astype("Float64")).replace(np.inf, 0)
    return warn_prop

# @st.cache_resource
# def aggview_warning_column_names(aggregated_view: pd.DataFrame):
#     return [(c1,c2) for c1,c2 in aggregated_view.columns if c2 != '']

@st.cache_resource
def aggview_warning_column_indexes(aggregated_view: pd.DataFrame, stray_combinations: list[tuple]):  
    return [aggregated_view.columns.get_loc(c) for c in stray_combinations]

@st.cache_resource
def aggview_warning_window_column_index_range(aggregated_view: pd.DataFrame, window_size, stray_combinations: list[tuple]):
    all_columns = [aggregated_view.columns.get_loc(c) for c in stray_combinations if c[0] == window_size]
    return all_columns[0], all_columns[-1]



class BackgroundColor:

    def __init__(self, aggregated_view, agg_view_style, combinations):
        self.aggregated_view = aggregated_view
        self.warn_prop = compute_warnings_ratio_df(aggregated_view, combinations)
        self.agg_view_style = agg_view_style
        self.combinations = combinations

    @staticmethod
    @st.cache_resource
    def background_color_palette():
        return sns.light_palette("#0080ff", as_cmap=True)

    def color_by_proportion(self):
        for idx in self.aggregated_view.index:
            self.agg_view_style = self.agg_view_style.background_gradient(cmap=BackgroundColor.background_color_palette(), gmap=self.warn_prop.loc[idx, self.combinations], subset=(idx, self.combinations), vmin=0, vmax=0.5, axis=1)

    def color_by_absolute(self):
        self.agg_view_style = self.agg_view_style.background_gradient(cmap=BackgroundColor.background_color_palette(), subset=self.combinations)

    @staticmethod
    def st_aggview_colormap_switch(disabled_range, enabled_range, disabled_callback, enabled_callback, callback_args):
        with stylable_container(key="toggle_gmap_div", css_styles="div {overflow: hidden;}"):
            with st.container(height=50, border=False, key="toggle_gmap_container"):
                col0, col1, col2, col3 = st.columns([6, 6, 1, 10])
                with col0:
                    st.write("Color cells by:")
                with col2:
                    toggle_value = st.toggle(label="whatever", label_visibility="hidden",  value=False, key="toggle_gmap_setting", help="Choose to color the aggregated view by absolute warning count or by percentage with respect to the number of sequences in the period of time.")
                if toggle_value:
                    with col1:
                        st.html('<p style="text-align: right; margin-right: 12px">Absolute amount of warnings</p>')
                    with col3:
                        st.html('<p style="text-align: left"><strong>Percentage of warnings in the period of time</strong>'
                                f'</br><small>Range: {enabled_range[0]} (gray) -  {enabled_range[1]} (bright blue)</small></p>')
                    enabled_callback(**callback_args)
                else:
                    with col1:
                        st.html('<p style="text-align: right; margin-right: 12px"><strong>Absolute amount of warnings</strong>'
                                f'</br><small>Range: {disabled_range[0]} (gray) - {disabled_range[1]} (bright blue)</small></p>')
                    with col3:
                        st.html('<p style="text-align: left">Percentage of warnings in the period of time</p>')
                    disabled_callback(**callback_args)

    def enable_aggview_background_color_switch(self):
        BackgroundColor.st_aggview_colormap_switch(
            disabled_range=[0, self.aggregated_view[self.combinations].max(axis=None)], 
            enabled_range=["0 %", "50 % or higher"], 
            disabled_callback=self.color_by_absolute, 
            enabled_callback=self.color_by_proportion, 
            callback_args={})



# WINDOW RECOMMENDATION
class StrayRecommendation:
    
    @staticmethod
    def st_recommender_selector(stray_window_advice_options, default_option):
        col0, col1 = st.columns([1, 2])
        with col0:
            st.markdown("<span>Recommended window for Stray:</br>(computed on unfiltered data)</span>", unsafe_allow_html=True)
        with col1:
            st.radio(label="align-this", label_visibility="collapsed", 
                     options=stray_window_advice_options, index=stray_window_advice_options.index(default_option), 
                     key="stray-window-recom-selector", 
                     help="Choose the method for highlighting the best window size for Stray",
                     horizontal=True)
            
    @staticmethod
    def insert_decision_metric_in_aggview(aggregated_view, avg_df: pd.DataFrame):
        aggregated_view.insert(loc=aggregated_view.columns.get_loc(('N. Sequences', ''))+1, 
                               column=("N. Sequences", "Unfiltered 4 Weeks Avg"), 
                               value=avg_df.loc[aggregated_view[('Week', '')],"mov.avg."].reset_index(drop=True))


class StrayWindowRecommendationAvg6w(StrayRecommendation):
    """
    Considering an aggregated view without filters applied, for any period of time of 2 weeks, computes the average number of sequences 
    in the sorrounding 6 weeks (2 past + 2 current + 2 future). 
    """

    @staticmethod
    def avg_6w(aggregated_view, combinations):
        """
        returns: a DataFrame with the recommended stray window size for each period of time. 
        In addition, week name of the period and two columns containing the index of the leftmost 
        and rightmost columns in the aggregated view corresponding to the recommended window.
        """
        num_sequences_moving_average = aggregated_view[('N. Sequences', '')].rolling(window=3, center=True, min_periods=1).mean().round(1).astype(str).apply(Decimal)
        suggested_stray_window = StrayWindowRecommendationAvg6w.recommended_stray_window_size(num_sequences_moving_average)
        recommended_column_index_range = [aggview_warning_window_column_index_range(aggregated_view, w, combinations) for w in suggested_stray_window]
        result = pd.concat([
            pd.DataFrame([aggregated_view[('Week', '')]]).reset_index(drop=True).rename({0: "week"}, axis=0).T,
            pd.DataFrame([*recommended_column_index_range], columns=["left_cell_idx", "right_idx"]),
            pd.DataFrame([num_sequences_moving_average, suggested_stray_window], index=["mov.avg.", "stray_window"]).T
        ], axis=1)
        return result
    
    @staticmethod
    def save_or_get_cached_avg_6w(aggregated_view, combinations, cache_dir_path):
        """
        Wrapper around StrayWindowRecommendationAvg6w.avg_6w that provides read/store from a cached file
        """
        os.makedirs(cache_dir_path, exist_ok=True)
        cache_path = os.path.join(cache_dir_path, "recommended_stray_window_size_mov_avg_6w.parquet")
        if not os.path.exists(cache_path):
            stray_recom_cache = StrayWindowRecommendationAvg6w.avg_6w(aggregated_view, combinations)
            stray_recom_cache.to_parquet(cache_path)
        return pd.read_parquet(cache_path).set_index("week")
            
    @staticmethod
    def insert_avg_6w_in_aggview(aggregated_view, avg_6w_df: pd.DataFrame):
        super().insert_decision_metric_in_aggview(aggregated_view, avg_6w_df)

    @staticmethod
    def recommended_stray_window_size(number_of_sequences: pd.Series):
        return pd.cut(number_of_sequences, [0, 5, 10, 50, 100], right=False, include_lowest=True, 
                    labels=[5, 10, 50, 100]).fillna(100)
    

class StrayWindowRecommendationAvg4wBefore(StrayRecommendation):
    """
    Considering an aggregated view without filters applied, for any period of time of 2 weeks, computes the average number of sequences 
    in the sorrounding current and past 2 weeks. 
    """

    @staticmethod
    def avg_4w_before(aggregated_view, combinations):
        """
        returns: a DataFrame with the recommended stray window size for each period of time. 
        In addition, week name of the period and two columns containing the index of the leftmost 
        and rightmost columns in the aggregated view corresponding to the recommended window.
        """
        num_sequences_moving_average = aggregated_view[('N. Sequences', '')].rolling(window=2, center=False, min_periods=1).mean().round(1).astype(str).apply(Decimal)
        suggested_stray_window = StrayWindowRecommendationAvg4wBefore.recommended_stray_window_size(num_sequences_moving_average)
        recommended_column_index_range = [aggview_warning_window_column_index_range(aggregated_view, w, combinations) for w in suggested_stray_window]
        result = pd.concat([
            pd.DataFrame([aggregated_view[('Week', '')]]).reset_index(drop=True).rename({0: "week"}, axis=0).T,
            pd.DataFrame([*recommended_column_index_range], columns=["left_cell_idx", "right_idx"]),
            pd.DataFrame([num_sequences_moving_average, suggested_stray_window], index=["mov.avg.", "stray_window"]).T
        ], axis=1)
        return result
    
    @staticmethod
    def save_or_get_cached_avg_4w_before(aggregated_view, combinations, cache_dir_path):
        """
        Wrapper around StrayWindowRecommendationAvg4wBefore.avg_4w_before that provides read/store from a cached file
        """
        os.makedirs(cache_dir_path, exist_ok=True)
        cache_path = os.path.join(cache_dir_path, "recommended_stray_window_size_mov_avg_4w_before.parquet")
        if not os.path.exists(cache_path):
            stray_recom_cache = StrayWindowRecommendationAvg4wBefore.avg_4w_before(aggregated_view, combinations)
            stray_recom_cache.to_parquet(cache_path)
        return pd.read_parquet(cache_path).set_index("week")

    @staticmethod
    def recommended_stray_window_size(number_of_sequences: pd.Series):
        return pd.cut(number_of_sequences, [0, 5, 10, 50, 100], right=False, include_lowest=True, 
                    labels=[5, 10, 50, 100]).fillna(100)
    

