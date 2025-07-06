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
PATENTS_JSON_DIR = './data_google_patents/patents_json'
GEMINI_OUTPUT_DIR = './data_google_patents/patents_gemini'
MAX_PATENTS_TO_PROCESS = 10

# --- Define the required JSON schema for validation ---
REQUIRED_SCHEMA = {
  "claimed_target": (str, type(None)),
  "claimed_application": (str, type(None)),
  "novelty_type": str,
  "evidence_summary": (str, type(None)),
  "is_therapeutic_claim": bool
}
ALLOWED_NOVELTY_TYPES = {'New-Sequence', 'New-Use-for-Existing-Antibody', 'New-Method', 'Unspecified'}


def get_patents_to_process(patents_json_dir, gemini_output_dir):
    """
    Lists patents that are in patents_json_dir but not in gemini_output_dir.
    """
    if not os.path.isdir(patents_json_dir):
        print(f"Error: Input directory not found: {patents_json_dir}")
        return []

    os.makedirs(gemini_output_dir, exist_ok=True)

    source_patents = {os.path.splitext(f)[0] for f in os.listdir(patents_json_dir) if f.endswith(".json")}
    processed_patents = {os.path.splitext(f)[0] for f in os.listdir(gemini_output_dir) if f.endswith(".json")}

    patents_to_process = sorted(list(source_patents - processed_patents)) # Sort for consistent order
    print(f"Found {len(source_patents)} total patents in source.")
    print(f"Found {len(processed_patents)} already processed patents.")
    print(f"{len(patents_to_process)} new patents to process.")

    return patents_to_process


def build_prompt(patent_data):
    """
    Constructs the prompt for the Gemini API call.
    Removes the 'description' field from the patent data as requested.
    """
    if 'description' in patent_data:
        del patent_data['description']

    patent_text = json.dumps(patent_data, indent=2)

    return f"""
Your task is to act as a data extraction engine. Analyze the following patent text. In your analysis, give the highest priority to information found under the headings 'Examples' and 'Summary of the Invention'.

Return a single, valid JSON object with the specified schema. Do not include any text, explanations, or markdown formatting before or after the JSON object.

Patent Text:

{patent_text}

Required JSON Schema:
{{
  "claimed_target": "The main antigen or molecule the invention is for. Null if not specified.",
  "claimed_application": "The specific disease or condition the invention is claimed to treat or diagnose. Null if not specified.",
  "novelty_type": "Classify the core invention. Choose one: 'New-Sequence', 'New-Use-for-Existing-Antibody', 'New-Method', 'Unspecified'.",
  "evidence_summary": "A brief, one-sentence summary of the experimental results or data presented in the 'Examples' section to support the claims.",
  "is_therapeutic_claim": "Return true if the text makes a direct claim about treating a disease, otherwise return false."
}}
"""


def validate_json_response(data):
    """
    Validates if the received data matches the required schema.
    """
    if not isinstance(data, dict):
        print(f"Validation Error: Response is not a dictionary.")
        return False

    if set(data.keys()) != set(REQUIRED_SCHEMA.keys()):
        print(f"Validation Error: JSON keys do not match schema. Got {set(data.keys())}, expected {set(REQUIRED_SCHEMA.keys())}")
        return False

    for key, expected_type in REQUIRED_SCHEMA.items():
        if not isinstance(data.get(key), expected_type):
            print(f"Validation Error: Key '{key}' has wrong type. Got {type(data.get(key))}, expected {expected_type}")
            return False

    if data.get('novelty_type') not in ALLOWED_NOVELTY_TYPES:
        print(f"Validation Error: 'novelty_type' has invalid value '{data.get('novelty_type')}'")
        return False

    return True


def process_patent(patent_id):
    """
    Processes a single patent: reads, calls Gemini, validates, and saves.
    """
    input_filepath = os.path.join(PATENTS_JSON_DIR, f"{patent_id}.json")
    output_filepath = os.path.join(GEMINI_OUTPUT_DIR, f"{patent_id}.json")

    try:
        print(f"Processing patent: {patent_id}")
        with open(input_filepath, 'r', encoding='utf-8') as f:
            patent_data = json.load(f)

        prompt = build_prompt(patent_data)

        model = genai.GenerativeModel(
            MODEL_NAME,
            generation_config={"response_mime_type": "application/json"}
        )
        response = model.generate_content(prompt)

        try:
            gemini_json = json.loads(response.text)
        except json.JSONDecodeError:
            print(f"Error: Failed to decode JSON from Gemini for patent {patent_id}.")
            print(f"Received text: {response.text[:500]}...")
            return False

        if not validate_json_response(gemini_json):
            print(f"Error: Gemini response for patent {patent_id} failed schema validation.")
            print(f"Received JSON: {gemini_json}")
            return False

        with open(output_filepath, 'w', encoding='utf-8') as f:
            json.dump(gemini_json, f, indent=4)

        print(f"Successfully processed and saved patent {patent_id}")
        return True

    except FileNotFoundError:
        print(f"Error: Patent file not found: {input_filepath}")
        return False
    except Exception as e:
        print(f"An unexpected error occurred while processing patent {patent_id}: {e}")
        return False


def main():
    """
    Main function to orchestrate the patent processing workflow.
    """
    patents_to_process = get_patents_to_process(PATENTS_JSON_DIR, GEMINI_OUTPUT_DIR)

    if not patents_to_process:
        print("No new patents to process. Exiting.")
        return

    patents_for_this_run = patents_to_process[:MAX_PATENTS_TO_PROCESS]
    print(f"\nWill process the first {len(patents_for_this_run)} patents...")

    success_count = 0
    failure_count = 0

    for patent_id in patents_for_this_run:
        if process_patent(patent_id):
            success_count += 1
        else:
            failure_count += 1
        time.sleep(1)  # Simple rate limiting to avoid overwhelming the API

    print("\n--- Processing Complete ---")
    print(f"Successfully processed: {success_count}")
    print(f"Failed to process: {failure_count}")


if __name__ == "__main__":
    main()
