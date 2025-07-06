import streamlit as st
import pandas as pd
import re
import json
import os

st.set_page_config(page_title="PLAbDab Plus", page_icon="🧬")
st.title('🧬 PLAbDab Plus Database')
st.info('The Extended Patent and Literature Antibody Database')

@st.cache_data
def load_data():
    df = pd.read_csv('./data.csv')
    return df

df = load_data()

search_term = st.text_input("Search by num or ID", "")

if search_term:
    # Search in 'num' (as string) and 'ID' columns.
    search_mask = (df['num'].astype(str).str.contains(search_term, case=False, na=False) | 
                   df['ID'].str.contains(search_term, case=False, na=False))
    results = df[search_mask]

    if results.empty:
        st.warning("No results found.")
    else:
        st.success(f"Found {len(results)} result(s).")
        for index, row in results.iterrows():
            with st.expander(f"Result: {row['ID']} (Num: {row['num']})", expanded=True):
                st.write("### Details")
                st.dataframe(row.to_frame().T)

                url = row.get('url')
                if not isinstance(url, str):
                    st.write("No associated URL or URL is not a string.")
                    continue

                st.write("### Associated Data")

                try:
                    if 'patents.google.com' in url:
                        match = re.search(r'patent/([^/]+)', url)
                        if match:
                            patent_id = match.group(1)
                            file_path = f'./data_google_patents/patents_json/{patent_id}.json'
                            if os.path.exists(file_path):
                                with open(file_path, 'r') as f:
                                    st.json(json.load(f))
                            else:
                                st.error(f"File not found: {file_path}")
                        else:
                            st.warning("Could not extract patent ID from URL.")

                    elif 'www.ncbi.nlm.nih.gov/protein/' in url:
                        protein_id = url.split('/')[-1]
                        file_path = f'./data_google_patents/data_ncbi/{protein_id}.json'
                        if os.path.exists(file_path):
                            with open(file_path, 'r') as f:
                                st.json(json.load(f))
                        else:
                            st.error(f"File not found: {file_path}")

                    elif 'sabdab/structureviewer' in url:
                        match = re.search(r'pdb=([^&]+)', url)
                        if match:
                            pdb_id = match.group(1)
                            file_path = f'./data_google_patents/data_sabdab/{pdb_id}.tsv'
                            if os.path.exists(file_path):
                                sabdab_df = pd.read_csv(file_path, sep='\\t', engine='python')
                                st.dataframe(sabdab_df)
                            else:
                                st.error(f"File not found: {file_path}")
                        else:
                            st.warning("Could not extract PDB ID from URL.")
                    else:
                        st.info("URL does not match known patterns for additional data.")

                except Exception as e:
                    st.error(f"An error occurred while trying to display associated data: {e}")
