import pandas as pd
import requests
import os
import time
import logging

# Setup basic logging
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')

def crawl_google_patents():
    """
    Crawls Google Patents URLs from a CSV file, downloads the patent pages,
    and updates the CSV.
    """
    csv_file_path = 'data.csv'
    output_dir = './data_google_patents/'

    os.makedirs(output_dir, exist_ok=True)

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
    rows_to_process = []

    for index, row in df.iterrows():
        is_google_patent = str(row.get('url', '')).startswith('https://patents.google.com')
        is_target_identified = row.get('targets_mentioned') != 'Unidentified'
        is_not_crawled = 'google_patents' not in str(row.get('crawl', ''))

        if is_google_patent and is_target_identified and is_not_crawled:
            rows_to_process.append(index)

    if not rows_to_process:
        logging.info("No new patents to crawl.")
        return

    logging.info(f"Found {len(rows_to_process)} patents to crawl.")

    for index in rows_to_process:
        row = df.loc[index]
        try:
            num = row['num']
            url = row['url']
            
            # logging.info(f"Crawling row num: {num}, url: {url}")

            headers = {
                'User-Agent': 'Mozilla/5.0 (Windows NT 10.0; Win64; x64) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/91.0.4472.124 Safari/537.36'
            }
            response = requests.get(url, headers=headers, timeout=30)
            response.raise_for_status()

            file_path = os.path.join(output_dir, f"{num}.html")
            with open(file_path, 'w', encoding='utf-8') as f:
                f.write(response.text)

            current_crawl = df.loc[index, 'crawl']
            if not current_crawl:
                df.loc[index, 'crawl'] = 'google_patents'
            else:
                df.loc[index, 'crawl'] = f"{current_crawl};google_patents"
            
            crawled_count += 1
            
            # logging.info(f"Successfully crawled and saved {file_path}")

            if crawled_count > 0 and crawled_count % 500 == 0:
                logging.info(f"Crawled {crawled_count} rows, saving CSV...")
                df.to_csv(csv_file_path, index=False)
                logging.info("CSV saved.")

            break
            time.sleep(0.2)

        except requests.exceptions.RequestException as e:
            logging.error(f"Error downloading {row.get('url', 'N/A')} for num: {row.get('num', 'N/A')}, ID: {row.get('ID', 'N/A')}. Error: {e}")
        except IOError as e:
            logging.error(f"Error saving file for num: {row.get('num', 'N/A')}, ID: {row.get('ID', 'N/A')}. Error: {e}")
        except Exception as e:
            logging.error(f"An error occurred processing row num: {row.get('num', 'N/A')}, ID: {row.get('ID', 'N/A')}. Error: {e}")

    if crawled_count > 0:
        logging.info("Crawling finished. Saving final CSV.")
        df.to_csv(csv_file_path, index=False)
        logging.info("Final CSV saved.")

if __name__ == "__main__":
    crawl_google_patents()
