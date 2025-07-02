import csv

csv_path = 'data.csv'
url_prefix = 'https://patents.google.com'
count_url = 0
count_url_and_target = 0
count_url_target_model = 0

with open(csv_path, newline='', encoding='utf-8') as csvfile:
    reader = csv.DictReader(csvfile)
    for row in reader:
        url = row.get('url', '')
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
