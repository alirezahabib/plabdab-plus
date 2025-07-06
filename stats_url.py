import pandas as pd
import os
import re

# Load the CSV file
df = pd.read_csv('data.csv')

# Directory containing downloaded patent HTML files
patents_dir = './data_google_patents/patents/'

# Get set of already downloaded patent IDs (without .html extension)
if os.path.isdir(patents_dir):
    downloaded_patent_ids = {os.path.splitext(fname)[0] for fname in os.listdir(patents_dir) if fname.endswith('.html')}
else:
    downloaded_patent_ids = set()

# Regex to extract patent id from Google Patents URL
patent_url_re = re.compile(r'^https://patents\.google\.com/patent/([^/]+)/')

unique_urls = set()

for url in df['url'].dropna().unique():
    match = patent_url_re.match(url)
    if match:
        patent_id = match.group(1)
        if patent_id not in downloaded_patent_ids:
            unique_urls.add(url)

# Print the unique URLs
for url in sorted(unique_urls):
    print(url)
