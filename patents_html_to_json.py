import os
import json
from bs4 import BeautifulSoup

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

    # --- Extract Title ---
    title_tag = soup.find('h1', itemprop='pageTitle')
    if title_tag:
        patent_data['title'] = title_tag.text.strip()
    else:
         title_tag = soup.find('span', itemprop='title')
         if title_tag:
              patent_data['title'] = title_tag.text.strip()


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
    patent_data['inventors'] = inventors

    # --- Extract Assignee ---
    assignee_tag = soup.find('dd', itemprop='assigneeCurrent')
    if assignee_tag:
        patent_data['assignee'] = assignee_tag.text.strip()

    # --- Extract Dates ---
    dates = {}
    priority_date_tag = soup.find('time', itemprop='priorityDate')
    if priority_date_tag:
        dates['priority_date'] = priority_date_tag['datetime']

    filing_date_tag = soup.find('time', itemprop='filingDate')
    if filing_date_tag:
        dates['filing_date'] = filing_date_tag['datetime']

    publication_date_tag = soup.find('time', itemprop='publicationDate')
    if publication_date_tag:
        dates['publication_date'] = publication_date_tag['datetime']
    patent_data['dates'] = dates


    # --- Extract Claims ---
    claims_section = soup.find('section', itemprop='claims')
    if claims_section:
        claims = []
        for claim_div in claims_section.find_all('div', class_='claim'):
            claim_text = claim_div.get_text(separator=' ', strip=True)
            claims.append(claim_text)
        patent_data['claims'] = claims

    # --- Extract Description ---
    description_section = soup.find('section', itemprop='description')
    if description_section:
        description_content = description_section.find('div', itemprop='content')
        if description_content:
            patent_data['description'] = description_content.get_text(separator='\n', strip=True)


    return patent_data

def convert_html_to_json(patents_dir, patents_json_dir):
    """
    Converts all HTML patent files in a directory to JSON format.

    Args:
        patents_dir (str): The directory containing the HTML patent files.
        patents_json_dir (str): The directory where JSON files will be saved.
    """
    if not os.path.exists(patents_json_dir):
        os.makedirs(patents_json_dir)
        print(f"Created directory: {patents_json_dir}")

    for filename in os.listdir(patents_dir):
        if filename.endswith(".html"):
            html_filepath = os.path.join(patents_dir, filename)
            json_filename = os.path.splitext(filename)[0] + ".json"
            json_filepath = os.path.join(patents_json_dir, json_filename)

            print(f"Processing {html_filepath}...")

            with open(html_filepath, 'r', encoding='utf-8') as f:
                html_content = f.read()

            patent_data = parse_patent_html(html_content)

            with open(json_filepath, 'w', encoding='utf-8') as f:
                json.dump(patent_data, f, indent=4, ensure_ascii=False)

            # print(f"Successfully converted to {json_filepath}")

if __name__ == '__main__':
    # Define the input and output directories
    patents_directory = './data_google_patents/patents'
    json_output_directory = './data_google_patents/patents_json'

    # Create the input directory and a dummy file for testing if they don't exist
    if not os.path.exists(patents_directory):
        os.makedirs(patents_directory)
        dummy_html_path = os.path.join(patents_directory, 'US10597465.html')
        if not os.path.exists(dummy_html_path):
             # NOTE: You should replace this with the actual content of your HTML file.
             # This is a placeholder to make the script runnable.
            with open(dummy_html_path, 'w', encoding='utf-8') as f:
                f.write('<html><head><title>US10597465B2 - Methods for the generation of multispecific and multivalent antibodies - Google Patents</title></head><body><section itemprop="abstract"><div itemprop="content">The invention provides novel bispecific monoclonal antibodies...</div></section></body></html>')
            print(f"Created a dummy file for testing: {dummy_html_path}")


    convert_html_to_json(patents_directory, json_output_directory)
