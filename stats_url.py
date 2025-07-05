import csv

csv_path = 'data.csv'
url_prefix = 'https://www.ncbi.nlm.nih.gov/protein/'
url_prefix_google_patents = 'https://patents.google.com'
unique_urls = set()
count_google_patents_without_crawl = 0

with open(csv_path, newline='', encoding='utf-8') as csvfile:
    reader = csv.DictReader(csvfile)
    for row in reader:
        url = row.get('url', '')
        if url.startswith(url_prefix):
            unique_urls.add(url)
        # Count for Google Patents URLs without 'google_patents' in crawl
        if url.startswith(url_prefix_google_patents):
            crawl_field = row.get('crawl', '')
            if 'google_patents' not in crawl_field.split(';'):
                count_google_patents_without_crawl += 1

print(f"Number of unique URLs starting with '{url_prefix}': {len(unique_urls)}")
print(f"Number of rows with url starting with '{url_prefix_google_patents}' and 'google_patents' not in crawl field: {count_google_patents_without_crawl}")
