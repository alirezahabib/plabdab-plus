import pandas as pd
import os
import requests
from urllib.parse import urlparse, parse_qs
from tqdm import tqdm

def main():
    # Create the output directory if it doesn't exist
    output_dir = './data_sabdab'
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    # Load existing PDB IDs from the output directory
    existing_ids = set()
    for filename in os.listdir(output_dir):
        if filename.endswith('.tsv'):
            existing_ids.add(filename.split('.')[0])
    
    print(f'Found {len(existing_ids)} existing files.')

    # Load the data from the CSV file
    try:
        df = pd.read_csv('./data.csv')
    except FileNotFoundError:
        print("Error: './data.csv' not found. Please make sure the file exists.")
        return

    # Iterate over each row in the dataframe
    for index, row in tqdm(df.iterrows(), total=df.shape[0], desc="Processing rows"):
        url = row['url']
        if isinstance(url, str) and url.startswith('https://opig.stats.ox.ac.uk/webapps/newsabdab/sabdab/structureviewer/'):
            parsed_url = urlparse(url)
            query_params = parse_qs(parsed_url.query)
            pdb_id = query_params.get('pdb', [None])[0]

            if pdb_id and pdb_id not in existing_ids:
                download_url = f'https://opig.stats.ox.ac.uk/webapps/sabdab-sabpred/sabdab/summary/{pdb_id}/'
                try:
                    response = requests.get(download_url)
                    response.raise_for_status()  # Raise an exception for bad status codes
                    
                    # The content is HTML, but we will save it as TSV as requested.
                    # Usually summary pages are HTML. If the downloaded content is not what you expect
                    # you might need to find a direct download link for the TSV file.
                    output_path = os.path.join(output_dir, f'{pdb_id}.tsv')
                    with open(output_path, 'wb') as f:
                        f.write(response.content)
                    
                    existing_ids.add(pdb_id)

                except requests.exceptions.RequestException as e:
                    tqdm.write(f'Failed to download {download_url}: {e}')

if __name__ == '__main__':
    main()
