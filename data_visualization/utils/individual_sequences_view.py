import streamlit as st
from utils.utils import DateRange, Table
from utils.db_queries import CachedDB

def tabs_individual_sequences(tab_titles: list, tab2table_name_dict: dict, db_engine, 
                              table_shown_columns: list, date_range: DateRange, 
                              sel_host_type = [], sel_country_or_state = [], 
                              sel_locations = [], sel_genotypes = []):
    st.markdown("### Individual sequences")
    tabs = st.tabs(tab_titles)
    tabs = dict(zip(tab_titles, tabs))

    for tname in tab_titles:
        with tabs[tname]:
            table_name = tab2table_name_dict[tname]
            table = CachedDB.get_table_named_as(db_engine, table_name).reset_index(drop=False).rename(columns={"index":"IsolateID_Name"})
            # anticipate date_range filter to make sequences start from index 0 regardless of input filters
            table = table[(table.Collection_Date >= date_range.start.to_numpy()) & (table.Collection_Date < date_range.end.to_numpy())].reset_index(drop=True)
            row_selection = Table.filter(table, sel_host_type, sel_country_or_state, sel_locations, date_range, sel_genotypes)
            st.dataframe(
                row_selection[table_shown_columns],
                column_config={
                    "Collection_Date": st.column_config.DateColumn(
                        "Collection_Date", format="YYYY MMM DD"
                    )
                },
                use_container_width=True,
            )
            st.write(f"Rows: {row_selection.shape[0]}")

