import csv

csv_path = 'data.csv'
url_prefix = 'https://patents.google.com'
unique_urls = set()

with open(csv_path, newline='', encoding='utf-8') as csvfile:
    reader = csv.DictReader(csvfile)
    for row in reader:
        url = row.get('url', '')
        if url.startswith(url_prefix):
            unique_urls.add(url)

print(f"Number of unique URLs starting with '{url_prefix}': {len(unique_urls)}")
