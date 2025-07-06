import os
import json
import time
import google.generativeai as genai

# --- Configuration ---
# IMPORTANT: Set your Gemini API key as an environment variable
# export GEMINI_API_KEY="YOUR_API_KEY"
try:
    API_KEY = os.environ["GEMINI_API_KEY"]
except KeyError:
    raise ValueError("GEMINI_API_KEY environment variable not set. Please export your key.")

genai.configure(api_key=API_KEY)

MODEL_NAME = 'gemini-2.5-flash'
INPUT_DIR = './data_ncbi'
OUTPUT_DIR = './data_ncbi_gemini'
MAX_FILES_TO_PROCESS = 10

# --- Define the required JSON schema for validation ---
REQUIRED_SCHEMA = {
  "primary_target": (str, type(None)),
  "disease_context": (str, type(None)),
  "development_stage": str,
  "efficacy_indicators": list,
  "author_conclusion_strength": str
}
ALLOWED_DEV_STAGES = {'In-Vitro', 'Cell-Based-Assay', 'Animal-Model', 'Human-Trial', 'Approved-Therapeutic', 'Unspecified'}
ALLOWED_CONCLUSION_STRENGTHS = {'Exploratory', 'Promising_Pre-clinical', 'Clinical_Efficacy_Demonstrated', 'Unspecified'}


def get_files_to_process(input_dir, output_dir):
    """
    Lists UIDs that are in input_dir but not in output_dir.
    """
    if not os.path.isdir(input_dir):
        print(f"Error: Input directory not found: {input_dir}")
        return []

    os.makedirs(output_dir, exist_ok=True)

    source_files = {os.path.splitext(f)[0] for f in os.listdir(input_dir) if f.endswith(".json")}
    processed_files = {os.path.splitext(f)[0] for f in os.listdir(output_dir) if f.endswith(".json")}

    files_to_process = sorted(list(source_files - processed_files))
    print(f"Found {len(source_files)} total source files.")
    print(f"Found {len(processed_files)} already processed files.")
    print(f"{len(files_to_process)} new files to process.")

    return files_to_process


def extract_pubmed_data(filepath, uid):
    """
    Extracts title and abstract from the specific PubMed JSON structure.
    Returns (title, abstract) or (None, None) if not found.
    """
    try:
        with open(filepath, 'r', encoding='utf-8') as f:
            data = json.load(f)

        if "pubmed" not in data or not isinstance(data["pubmed"], dict):
            print(f"Skipping {uid}: 'pubmed' key missing or not a dictionary.")
            return None, None
            
        # The UID is a key within the 'pubmed' object
        if uid not in data["pubmed"]:
            # Sometimes the key is slightly different. Let's try to find it.
            if len(data["pubmed"].keys()) == 1:
                uid = list(data["pubmed"].keys())[0]
            else:
                 print(f"Skipping {uid}: UID key not found directly in 'pubmed' object.")
                 return None, None


        article = data.get("pubmed", {}).get(uid, {}).get("PubmedArticleSet", {}).get("PubmedArticle", {}).get("MedlineCitation", {}).get("Article", {})
        if not article:
            print(f"Skipping {uid}: Could not navigate to Article object.")
            return None, None

        title = article.get("ArticleTitle")
        abstract = article.get("Abstract", {}).get("AbstractText")

        if not title or not abstract:
            print(f"Skipping {uid}: Missing Title or Abstract.")
            return None, None

        return str(title), str(abstract)

    except json.JSONDecodeError:
        print(f"Skipping {uid}: Could not decode JSON from {filepath}")
        return None, None
    except Exception as e:
        print(f"An unexpected error occurred extracting data for {uid}: {e}")
        return None, None


def build_prompt(title, abstract):
    """Constructs the prompt for the Gemini API call."""
    return f"""
Your task is to act as a data extraction engine. Analyze the following PubMed abstract and return a single, valid JSON object with the specified schema. Do not include any text, explanations, or markdown formatting before or after the JSON object.

Title:
{title}

Abstract Text:
{abstract}

Required JSON Schema:
{{
  "primary_target": "The main antigen or molecule the antibody binds to. Null if not specified.",
  "disease_context": "The primary disease or condition being studied. Null if not specified.",
  "development_stage": "Classify the research stage. Choose one: 'In-Vitro', 'Cell-Based-Assay', 'Animal-Model', 'Human-Trial', 'Approved-Therapeutic', 'Unspecified'.",
  "efficacy_indicators": [
    "Create a list of specific keywords or short phrases from the text that demonstrate effectiveness (e.g., 'IC50 of 5nM', 'tumor regression', 'potent neutralization')."
  ],
  "author_conclusion_strength": "Classify the strength of the conclusion. Choose one: 'Exploratory' (describes a basic finding), 'Promising_Pre-clinical' (shows strong results in cells/animals and suggests therapeutic potential), 'Clinical_Efficacy_Demonstrated' (reports positive human trial results), or 'Unspecified'."
}}
"""


def validate_json_response(data):
    """Validates if the received data matches the required schema."""
    if not isinstance(data, dict):
        print("Validation Error: Response is not a dictionary.")
        return False

    if set(data.keys()) != set(REQUIRED_SCHEMA.keys()):
        print(f"Validation Error: JSON keys mismatch. Got {set(data.keys())}, expected {set(REQUIRED_SCHEMA.keys())}")
        return False

    for key, expected_type in REQUIRED_SCHEMA.items():
        if not isinstance(data.get(key), expected_type):
            print(f"Validation Error: Key '{key}' has wrong type. Got {type(data.get(key))}, expected {expected_type}")
            return False
        if key == 'efficacy_indicators' and not all(isinstance(item, str) for item in data[key]):
             print(f"Validation Error: Not all items in 'efficacy_indicators' are strings.")
             return False


    if data.get('development_stage') not in ALLOWED_DEV_STAGES:
        print(f"Validation Error: 'development_stage' has invalid value '{data.get('development_stage')}'")
        return False

    if data.get('author_conclusion_strength') not in ALLOWED_CONCLUSION_STRENGTHS:
        print(f"Validation Error: 'author_conclusion_strength' has invalid value '{data.get('author_conclusion_strength')}'")
        return False

    return True


def process_file(uid):
    """Processes a single file: reads, calls Gemini, validates, and saves."""
    input_filepath = os.path.join(INPUT_DIR, f"{uid}.json")
    output_filepath = os.path.join(OUTPUT_DIR, f"{uid}.json")

    print(f"--- Processing UID: {uid} ---")
    title, abstract = extract_pubmed_data(input_filepath, uid)

    if not title or not abstract:
        return False

    prompt = build_prompt(title, abstract)

    try:
        model = genai.GenerativeModel(
            MODEL_NAME,
            generation_config={"response_mime_type": "application/json"}
        )
        response = model.generate_content(prompt)

        try:
            gemini_json = json.loads(response.text)
        except json.JSONDecodeError:
            print(f"Error: Failed to decode JSON from Gemini for UID {uid}.")
            print(f"Received text: {response.text[:500]}...")
            return False

        if not validate_json_response(gemini_json):
            print(f"Error: Gemini response for UID {uid} failed schema validation.")
            print(f"Received JSON: {gemini_json}")
            return False

        with open(output_filepath, 'w', encoding='utf-8') as f:
            json.dump(gemini_json, f, indent=4)

        print(f"Successfully processed and saved UID {uid}")
        return True

    except Exception as e:
        print(f"An unexpected error occurred while processing UID {uid}: {e}")
        return False


def main():
    """Main function to orchestrate the processing workflow."""
    files_to_process = get_files_to_process(INPUT_DIR, OUTPUT_DIR)

    if not files_to_process:
        print("No new files to process. Exiting.")
        return

    files_for_this_run = files_to_process[:MAX_FILES_TO_PROCESS]
    print(f"\nWill process the first {len(files_for_this_run)} files...")

    success_count = 0
    failure_count = 0

    for uid in files_for_this_run:
        if process_file(uid):
            success_count += 1
        else:
            failure_count += 1
        time.sleep(1) # Simple rate limiting

    print("\n--- Processing Complete ---")
    print(f"Successfully processed: {success_count}")
    print(f"Failed: {failure_count}")


if __name__ == "__main__":
    main()
