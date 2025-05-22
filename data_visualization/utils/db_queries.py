import pandas as pd
import streamlit as st
import sqlalchemy as sa

# class Tables:
#     input_data = "input_data"

# def query_input_sequences(host_type: list = [], location: list =[]):
#     base_query = f"SELECT * FROM {Tables.input_data} "
#     if not host_type or location:
#         return base_query
#     where_stmt = "WHERE "
#     if host_type:
#         "host_type IN ("+",".join(",")

#     where_stmt = "" if not kwagrs else "WHERE "
#     for k:v in kwagrs.items():
#         where_stmt += f"{k}==v" 

class CachedDB:
    
    @staticmethod
    def host_type_values(db_engine: sa.engine.Engine):
        with db_engine.connect() as _conn:
            return pd.read_sql_table(table_name="input_data", con=_conn).Host_Type.unique().tolist()

    @staticmethod
    def location_values(db_engine: sa.engine.Engine):
        with db_engine.connect() as _conn:
            return pd.read_sql_table(table_name="input_data", con=_conn).state_or_country.unique().tolist()

    @staticmethod
    def get_input_data_table(db_engine: sa.engine.Engine):
        with db_engine.connect() as _conn:
            return pd.read_sql_table(table_name="input_data", con=_conn)

    @staticmethod
    def get_warning_detail_table(db_engine: sa.engine.Engine, window_size, k):
        with db_engine.connect() as _conn:
            return pd.read_sql_table(table_name=f"window{window_size}_k{k}", con=_conn)

    @staticmethod
    def get_table_named_as(db_engine: sa.engine.Engine, name):
        with db_engine.connect() as _conn:
            return pd.read_sql_table(table_name=name, con=_conn)


class DB:
    
    @staticmethod
    def host_type_values(_db_engine: sa.engine.Engine):
        with _db_engine.connect() as _conn:
            return pd.read_sql_table(table_name="input_data", con=_conn).Host_Type.unique().tolist()

    @staticmethod
    def location_values(_db_engine: sa.engine.Engine):
        with _db_engine.connect() as _conn:
            return pd.read_sql_table(table_name="input_data", con=_conn).state_or_country.unique().tolist()

    @staticmethod
    def get_input_data_table(_db_engine: sa.engine.Engine):
        with _db_engine.connect() as _conn:
            return pd.read_sql_table(table_name="input_data", con=_conn)

    @staticmethod
    def get_warning_detail_table(_db_engine: sa.engine.Engine, window_size, k):
        with _db_engine.connect() as _conn:
            return pd.read_sql_table(table_name=f"window{window_size}_k{k}", con=_conn)

    @staticmethod
    def get_table_named_as(_db_engine: sa.engine.Engine, name):
        with _db_engine.connect() as _conn:
            return pd.read_sql_table(table_name=name, con=_conn)