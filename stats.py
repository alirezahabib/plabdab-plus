import csv
from collections import Counter

csv_path = 'data.csv'
url_prefix = 'https://patents.google.com'
count_url = 0
count_url_and_target = 0
count_url_target_model = 0
urls = []

with open(csv_path, newline='', encoding='utf-8') as csvfile:
    reader = csv.DictReader(csvfile)
    for row in reader:
        url = row.get('url', '')
        urls.append(url)
        targets_mentioned = row.get('targets_mentioned', '')
        model = row.get('model', '')
        if url.startswith(url_prefix):
            count_url += 1
            if targets_mentioned != 'Unidentified':
                count_url_and_target += 1
                if model != 'FAILED':
                    count_url_target_model += 1

print(f"Entries with url starting with '{url_prefix}': {count_url}")
print(f"Entries with url starting with '{url_prefix}' and targets_mentioned not 'Unidentified': {count_url_and_target}")
print(f"Entries with url starting with '{url_prefix}', targets_mentioned not 'Unidentified', and model not 'FAILED': {count_url_target_model}")

url_counts = Counter(urls)
duplicates = {url: count for url, count in url_counts.items() if count > 1}
num_unique_duplicates = len(duplicates)

print(f"\nFound {num_unique_duplicates} unique URLs that have duplicates.")

if num_unique_duplicates > 0:
    total_duplicate_entries = sum(count - 1 for count in duplicates.values())
    print(f"This corresponds to {total_duplicate_entries} duplicate rows in total.")

    print("\nFirst 5 duplicate URLs:")
    i = 0
    for url, count in duplicates.items():
        if i >= 5:
            break
        print(f"  {url}: {count}")
        i += 1
