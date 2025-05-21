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


st.navigation([st.Page("navigation/h1n1_08-09.py", title="H1N1 08-09"),
               st.Page("navigation/h5n1_19-25.py", title="H5N1 19-25"),
               st.Page("navigation/ss_v10.1.py", title="Sensitivity & Specificity")
                ]).run()

os.makedirs("generated/logs", exist_ok=True)
with open("generated/logs/sessions.log", "a") as f:
    f.write("User interaction at {}\n".format(datetime.now().strftime('%Y %b %d (%a) - %H:%M\' %S\"')))