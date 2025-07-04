import pandas as pd
import requests
import os
import time
import logging
import re

# Setup basic logging
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')

def crawl_google_patents():
    """
    Crawls Google Patents URLs from a CSV file, downloads the patent pages,
    and updates the CSV.
    """
    csv_file_path = 'data.csv'
    output_dir = './data_google_patents/patents/'

    os.makedirs(output_dir, exist_ok=True)

    # Collect already downloaded patent IDs based on existing HTML files in the output directory
    try:
        downloaded_patent_ids = {os.path.splitext(fname)[0] for fname in os.listdir(output_dir) if fname.endswith('.html')}
    except FileNotFoundError:
        downloaded_patent_ids = set()

    try:
        df = pd.read_csv(csv_file_path)
    except FileNotFoundError:
        logging.error(f"Error: {csv_file_path} not found.")
        return
    except Exception as e:
        logging.error(f"Error loading {csv_file_path}: {e}")
        return

    # Rename the first column to 'num' if it's not already named 'num'
    if df.columns[0] != 'num':
        logging.info(f"Renaming first column '{df.columns[0]}' to 'num'.")
        df.rename(columns={str(df.columns[0]): 'num'}, inplace=True)

    # Ensure 'crawl' column exists
    if 'crawl' not in df.columns:
        df['crawl'] = ''

    # Fill NaN values in 'crawl' column with empty strings for easier processing
    df['crawl'] = df['crawl'].fillna('')

    crawled_count = 0

    for index, row in df.iterrows():
        url = str(row.get('url', ''))

        # Extract the patent ID from URLs such as https://patents.google.com/patent/US5830675/en
        match = re.match(r'https?://patents\.google\.com/patent/([^/]+)/?', url)
        if not match:
            continue  # Skip rows that are not Google Patent URLs in the expected format

        patent_id = match.group(1)

        # Skip if this patent has already been downloaded in a previous run or earlier in this run
        if patent_id in downloaded_patent_ids:
            continue

        try:
            logging.info(f"Crawling patent_id: {patent_id}, url: {url}")

            headers = {
                'User-Agent': 'Mozilla/5.0 (Windows NT 10.0; Win64; x64) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/91.0.4472.124 Safari/537.36'
            }
            response = requests.get(url, headers=headers, timeout=30)
            response.raise_for_status()

            file_path = os.path.join(output_dir, f"{patent_id}.html")
            with open(file_path, 'w', encoding='utf-8') as f:
                f.write(response.text)

            # Update crawl column
            current_crawl = str(df.loc[index, 'crawl'])
            if 'google_patents' not in current_crawl.split(';'):
                df.loc[index, 'crawl'] = (current_crawl + ';' if current_crawl else '') + 'google_patents'

            crawled_count += 1

            # Add to set so subsequent rows in this run skip downloading the same patent
            downloaded_patent_ids.add(patent_id)

            if crawled_count > 0 and crawled_count % 500 == 0:
                logging.info(f"Crawled {crawled_count} rows, saving CSV...")
                df.to_csv(csv_file_path, index=False)
                logging.info("CSV saved.")

            time.sleep(0.1)

        except requests.exceptions.RequestException as e:
            logging.error(f"Error downloading {url} for num: {row.get('num', 'N/A')}, ID: {row.get('ID', 'N/A')}. Error: {e}")
        except IOError as e:
            logging.error(f"Error saving file for patent_id: {patent_id}, num: {row.get('num', 'N/A')}, ID: {row.get('ID', 'N/A')}. Error: {e}")
        except Exception as e:
            logging.error(f"An error occurred processing patent_id: {patent_id}, num: {row.get('num', 'N/A')}, ID: {row.get('ID', 'N/A')}. Error: {e}")

    if crawled_count > 0:
        logging.info("Crawling finished. Saving final CSV.")
        df.to_csv(csv_file_path, index=False)
        logging.info("Final CSV saved.")

if __name__ == "__main__":
    crawl_google_patents()
