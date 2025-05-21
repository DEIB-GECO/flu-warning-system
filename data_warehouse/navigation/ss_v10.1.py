import streamlit as st
import pandas as pd
from utils.aggregated_view import *

st.markdown("**Input data**")
st.markdown("Data source: H1N1 sequences as clustered in [A codon usage-based approach for the stratification of Influenza A; Alfonsi T, Chiara M., Bernasconi A.; https://doi.org/10.21203/rs.3.rs-5737660/v1](https://doi.org/10.21203/rs.3.rs-5737660/v1).")
st.markdown("- G1: 100 human sequences within the cluster 1 and collection date period 2009/01/26 - 2009/02/20 (25 days)")
st.markdown("- G2: 100 human sequences within the cluster 0, Lineage assigned and collection date period 2009/12/14 - 2010/01/10 (27 days)")
st.markdown("**Noise**")
st.markdown("- noise proportions: 10/20/30%")
st.markdown("- noise source: 30 human sequences within the cluster 1 and collection date period 2009/02/20 - 2009/03/01 (9 days)")
st.markdown("- noise insertion method: a portion of the sequences og G2 gets replaced with noise sequences at random positions")


st.markdown("### Senstivity and Specificity")
st.latex(r"sensitivity = \frac{K}{PLW} \,\,\,\,\,\,\,\,\,  specificity=\frac{N+1-PLW}{N+1-K}")

# table columns ["N", "K", "Noise_perc", "Noise", "PLW", "N.Warnings", "PLW_global", "N.Warnings_global",  "LNP", "LNP_global", "idx_G1>=G2", "NP"]
table = (pd.read_parquet("data/sensitivity_specificity10.parquet")
         .drop(columns=['PLW_global', 'N.Warnings_global', "LNP", "LNP_global", "idx_G1>=G2", "NP"])
         .rename(columns={'Noise_perc': '% Noise', 'N.Warnings': 'Warnings'})
         )
table[['% Noise']] = table[['% Noise']].map(lambda x: f"{int(x*100)} %")

# table = table.astype({x: int for x in ['N', 'K', 'Noise_perc', 'Noise', 'PLW', 'N.Warnings']})

combinations = table[['N', 'K', '% Noise']].sort_values(by=['N', 'K', '% Noise']).value_counts().index


# sensitivity / specificity
table['Sensitivity'] = table.K / table.PLW
table['Specificity'] = (table.N + 1 - table.PLW) / (table.N +1 -table.K)

# average values
average_values = table.groupby(by=["N", "K", "% Noise", "Noise"]).mean().drop(index=[(5, 3, "20 %", 1)])
average_values_cleaned_index = average_values.index.drop([(5, 3,), (10, 5)])
average_values = average_values.loc[average_values_cleaned_index]


# cosemtics
def pretty_html_table(t: pd.DataFrame):    
    t_style = t.style.set_properties(**{'background-color': 'white'})
    t_style = (t_style.set_sticky(axis="columns"))
    t_style = t_style.format('{:.2f}', na_rep='MISS', subset=['PLW', 'Warnings', 'Sensitivity', 'Specificity'])
    html_table = t_style.to_html()
    htmltable_helper = XMLTableEditor(html_table)
    XMLTableEditor.add_attribute(htmltable_helper.table_root, "style", "width:100%")
    return htmltable_helper.format_with_style()


st.markdown("### Average Sensitivity / Specificity (50 runs)")
st.html(pretty_html_table(average_values))


