import streamlit_authenticator as stauth
import yaml
from yaml.loader import SafeLoader

with open('config.yaml') as file:
    config = yaml.load(file, Loader=SafeLoader)

# This will hash all plain text passwords in the config file
stauth.Hasher.hash_passwords(config['credentials'])

with open('config.yaml', 'w') as file:
    yaml.dump(config, file, default_flow_style=False)

print("Password for 'user' has been hashed and config.yaml is updated.") 