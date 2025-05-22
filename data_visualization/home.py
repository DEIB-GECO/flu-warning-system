import streamlit as st
from utils.utils import StreamlitConfig
import os
from datetime import datetime

# Config Page (must be the first function called)
st.set_page_config(layout="wide", initial_sidebar_state="expanded")

# hide streamlit top bar if it is set in the config file
StreamlitConfig.autohide_streamlit_top_bar()

# Header of streamlit
st.header("Data Warehouse for the Flu Warning System")


st.navigation([st.Page("navigation/h1n1.py", title="H1N1 08-09"),
               st.Page("navigation/h5n1.py", title="H5N1 19-25")
                ]).run()
