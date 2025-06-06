import streamlit as st
import streamlit_authenticator as stauth
import yaml
from yaml.loader import SafeLoader
import mariadb
import pandas as pd
from datetime import date

st.set_page_config(page_title="PLAbDab", page_icon="🔒", layout="wide")

# --- DATABASE CONNECTION ---
@st.cache_resource
def get_db_connection():
    try:
        conn = mariadb.connect(
            user="root",
            password="12341234",
            host="127.0.0.1",
            port=3306,
            database="plabdab"
        )
        return conn
    except mariadb.Error as e:
        st.error(f"Error connecting to MariaDB Platform: {e}")
        return None

# --- USER AUTHENTICATION ---
with open('config.yaml') as file:
    config = yaml.load(file, Loader=SafeLoader)

authenticator = stauth.Authenticate(
    config['credentials'],
    config['cookie']['name'],
    config['cookie']['key'],
    config['cookie']['expiry_days']
)

authenticator.login()

if st.session_state["authentication_status"] is False:
    st.error('Username/password is incorrect')
elif st.session_state["authentication_status"] is None:
    st.warning('Please enter your username and password')
elif st.session_state["authentication_status"]:
    st.sidebar.title(f'Welcome *{st.session_state["name"]}*')
    st.sidebar.write(f'Username: `{st.session_state["username"]}`') 
    authenticator.logout('Logout', 'sidebar')
    st.sidebar.header('PLAbDab')
    st.sidebar.info('The Patent and Literature Antibody Database')
    st.sidebar.markdown('---')
    st.sidebar.subheader('Project Info')
    st.sidebar.markdown('**Author:** Alireza Habibzadeh')
    st.sidebar.markdown('**Supervisor:** Dr. Ali Sharifi-Zarchi')
    st.sidebar.markdown('**Sharif University of Technology**')
    

    st.sidebar.markdown('---')
    action = st.sidebar.selectbox('User Actions', ['None', 'Update user details', 'Reset password'])

    if action == 'Update user details':
        try:
            if authenticator.update_user_details(st.session_state['username']):
                st.success('Entries updated successfully')
                # Update config file
                with open('config.yaml', 'w') as file:
                    yaml.dump(config, file, default_flow_style=False)
        except Exception as e:
            st.error(e)

    elif action == 'Reset password':
        try:
            if authenticator.reset_password(st.session_state['username']):
                st.success('Password modified successfully')
                # Update config file
                with open('config.yaml', 'w') as file:
                    yaml.dump(config, file, default_flow_style=False)
        except Exception as e:
            st.error(e)

    st.title('PLAbDab Database Search')

    conn = get_db_connection()
    if conn:
        cursor = conn.cursor(dictionary=True)

        # --- SEARCH AND FILTER UI ---
        st.header('Search and Filter')
        search_id = st.text_input('Search by ID')
        
        today = date.today()
        default_start_date = date(2000, 1, 1)
        
        col1, col2 = st.columns(2)
        with col1:
            start_date = st.date_input('Start date', value=default_start_date, min_value=date(1900, 1, 1), max_value=today)
        with col2:
            end_date = st.date_input('End date', value=today, min_value=start_date, max_value=today)

        if st.button('Search'):
            query = "SELECT * FROM paired_sequences WHERE 1=1"
            params = []
            if search_id:
                query += " AND ID LIKE %s"
                params.append(f"%{search_id}%")
            if start_date and end_date:
                query += " AND update_date BETWEEN %s AND %s"
                params.append(start_date)
                params.append(end_date)
            
            cursor.execute(query, tuple(params))
            results = cursor.fetchall()
            
            if results:
                df = pd.DataFrame(results)
                st.session_state['search_results'] = df
            else:
                st.warning('No results found.')
                st.session_state['search_results'] = pd.DataFrame() # Empty dataframe
    else:
        st.error("Database connection is not available.")


    # --- DISPLAY RESULTS ---
    if 'search_results' in st.session_state and not st.session_state.search_results.empty:
        st.header('Search Results')
        df_results = st.session_state.search_results
        st.dataframe(df_results[['ID', 'organism', 'update_date', 'pairing', 'targets_mentioned']])

        # --- EXPANDED VIEW ---
        st.header('Detailed Information')
        selected_id = st.selectbox('Select an ID to see full details', options=df_results['ID'].tolist())

        if selected_id:
            selected_item = df_results[df_results['ID'] == selected_id].iloc[0]
            with st.expander(f"Details for {selected_id}", expanded=True):
                for col, value in selected_item.items():
                    st.write(f"**{col.replace('_', ' ').title()}:** {value}")