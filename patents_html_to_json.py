import os
import json
from bs4 import BeautifulSoup
import concurrent.futures
import time

def parse_patent_html(html_content):
    """
    Parses the HTML content of a Google Patents page to extract key information.

    Args:
        html_content (str): The HTML content of the patent page.

    Returns:
        dict: A dictionary containing the extracted patent information.
    """
    soup = BeautifulSoup(html_content, 'lxml')
    patent_data = {}

    # --- Extract and Clean Title ---
    title_tag = soup.find('h1', itemprop='pageTitle')
    if not title_tag:
        title_tag = soup.find('span', itemprop='title')

    if title_tag:
        full_title = title_tag.text.strip()
        patent_data['title'] = full_title.split('\n')[0].strip()
    
    # --- Extract Abstract ---
    abstract_section = soup.find('section', itemprop='abstract')
    if abstract_section:
        abstract_content = abstract_section.find('div', itemprop='content')
        if abstract_content:
            patent_data['abstract'] = abstract_content.get_text(separator='\n', strip=True)

    # --- Extract Publication Number ---
    pub_num_dd = soup.find('dd', itemprop='publicationNumber')
    if pub_num_dd:
        patent_data['publication_number'] = pub_num_dd.text.strip()

    # --- Extract Inventors ---
    inventors = [inventor.text.strip() for inventor in soup.find_all('dd', itemprop='inventor')]
    if inventors:
        patent_data['inventors'] = inventors

    # --- Extract Assignee ---
    assignee_tag = soup.find('dd', itemprop='assigneeCurrent')
    if assignee_tag:
        patent_data['assignee'] = assignee_tag.text.strip()

    # --- Extract Dates ---
    dates = {}
    priority_date_tag = soup.find('time', itemprop='priorityDate')
    if priority_date_tag and priority_date_tag.has_attr('datetime'):
        dates['priority_date'] = priority_date_tag['datetime']

    filing_date_tag = soup.find('time', itemprop='filingDate')
    if filing_date_tag and filing_date_tag.has_attr('datetime'):
        dates['filing_date'] = filing_date_tag['datetime']

    publication_date_tag = soup.find('time', itemprop='publicationDate')
    if publication_date_tag and publication_date_tag.has_attr('datetime'):
        dates['publication_date'] = publication_date_tag['datetime']
    
    if dates:
        patent_data['dates'] = dates

    # --- Extract Claims ---
    claims_section = soup.find('section', itemprop='claims')
    if claims_section:
        claims = []
        for claim_div in claims_section.find_all('div', class_='claim'):
            claim_text = claim_div.get_text(separator=' ', strip=True)
            claims.append(claim_text)
        if claims:
            patent_data['claims'] = claims

    # --- Extract Description ---
    description_section = soup.find('section', itemprop='description')
    if description_section:
        description_content = description_section.find('div', itemprop='content')
        if description_content:
            patent_data['description'] = description_content.get_text(separator='\n', strip=True)

    return patent_data

def process_file(html_filepath, json_output_dir):
    """
    Reads a single HTML file, parses it, and writes the JSON output.
    Includes error handling for individual file processing.

    Args:
        html_filepath (str): The full path to the input HTML file.
        json_output_dir (str): The directory to save the output JSON file.
    
    Returns:
        str: The path of the file processed, or None if an error occurred.
    """
    try:
        print(f"Processing {html_filepath}")
        json_filename = os.path.splitext(os.path.basename(html_filepath))[0] + ".json"
        json_filepath = os.path.join(json_output_dir, json_filename)

        with open(html_filepath, 'r', encoding='utf-8') as f:
            html_content = f.read()

        patent_data = parse_patent_html(html_content)

        # If parsing results in no data, it might be a corrupted file.
        if not patent_data:
            print(f"Warning: No data extracted from {html_filepath}. It may be empty or malformed. Skipping.")
            return None

        with open(json_filepath, 'w', encoding='utf-8') as f:
            json.dump(patent_data, f, indent=4, ensure_ascii=False)
        
        return html_filepath
    except Exception as e:
        # This will catch errors during file reading or parsing (e.g., from BeautifulSoup)
        print(f"Error processing {html_filepath} due to potential corruption: {e}. Skipping.")
        return None

def convert_html_to_json_parallel(patents_dir, patents_json_dir):
    """
    Converts all HTML patent files in a directory to JSON format using multithreading.

    Args:
        patents_dir (str): The directory containing the HTML patent files.
        patents_json_dir (str): The directory where JSON files will be saved.
    """
    if not os.path.exists(patents_json_dir):
        os.makedirs(patents_json_dir)
        print(f"Created directory: {patents_json_dir}")

    html_files = [os.path.join(patents_dir, f) for f in os.listdir(patents_dir) if f.endswith(".html")]
    
    if not html_files:
        print(f"No HTML files found in {patents_dir}.")
        return

    start_time = time.time()
    
    with concurrent.futures.ThreadPoolExecutor() as executor:
        futures = [executor.submit(process_file, filepath, patents_json_dir) for filepath in html_files]
        
        processed_count = 0
        error_count = 0
        for future in concurrent.futures.as_completed(futures):
            result = future.result()
            if result:
                processed_count += 1
            else:
                error_count += 1

    end_time = time.time()
    print(f"\n--- Conversion Complete ---")
    print(f"Successfully processed: {processed_count}/{len(html_files)} files")
    print(f"Skipped due to errors: {error_count}")
    print(f"Total time: {end_time - start_time:.2f} seconds.")


if __name__ == '__main__':
    patents_directory = './data_google_patents/patents'
    json_output_directory = './data_google_patents/patents_json'

    if not os.path.exists(patents_directory):
        os.makedirs(patents_directory)
        dummy_html_path = os.path.join(patents_directory, 'US10597465.html')
        if not os.path.exists(dummy_html_path):
            with open(dummy_html_path, 'w', encoding='utf-8') as f:
                f.write('<html><head><title>US10597465B2 - Methods for the generation of multispecific and multivalent antibodies \n - Google Patents</title></head><body><section itemprop="abstract"><div itemprop="content">The invention provides novel bispecific monoclonal antibodies...</div></section></body></html>')
            print(f"Created a dummy file for testing: {dummy_html_path}")

    convert_html_to_json_parallel(patents_directory, json_output_directory)
