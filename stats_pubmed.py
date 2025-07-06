import os
import json

def list_json_files(directory):
    for root, _, files in os.walk(directory):
        for file in files:
            if file.endswith('.json'):
                yield os.path.join(root, file)

def count_missing_pubmed(directory):
    total = 0
    missing_pubmed = 0
    for filepath in list_json_files(directory):
        total += 1
        try:
            with open(filepath, 'r', encoding='utf-8') as f:
                data = json.load(f)
            if 'pubmed' not in data:
                missing_pubmed += 1
        except Exception as e:
            print(f"Error reading {filepath}: {e}")
    return total, missing_pubmed

if __name__ == "__main__":
    directory = './data_ncbi'
    total, missing = count_missing_pubmed(directory)
    print(f"Total JSON files: {total}")
    print(f"Files missing 'pubmed' key: {missing}")
