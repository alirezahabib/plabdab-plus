import csv
import os
import re

csv_path = 'data.csv'
url_prefix = 'https://www.ncbi.nlm.nih.gov/protein/'
url_prefix_google_patents = 'https://patents.google.com'
unique_urls = set()
missing_patent_files_count = 0
patent_id_pattern = re.compile(r'https://patents\.google\.com/patent/([^/]+)/?')

with open(csv_path, newline='', encoding='utf-8') as csvfile:
    reader = csv.DictReader(csvfile)
    for row in reader:
        url = row.get('url', '')
        if url.startswith(url_prefix):
            unique_urls.add(url)
        # For Google Patents URLs, check for missing HTML file
        if url.startswith(url_prefix_google_patents):
            match = patent_id_pattern.match(url)
            if match:
                patent_id = match.group(1)
                html_path = f'./data_google_patents/patents/{patent_id}.html'
                if not os.path.isfile(html_path):
                    missing_patent_files_count += 1

print(f"Number of unique URLs starting with '{url_prefix}': {len(unique_urls)}")
print(f"Number of rows with url starting with '{url_prefix_google_patents}' and missing patent HTML file: {missing_patent_files_count}")
