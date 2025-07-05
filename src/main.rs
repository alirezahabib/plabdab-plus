use glob::glob;
use scraper::{Html, Selector};
use serde::Serialize;
use std::fs;
use std::path::Path;
use tokio::task;

/// Represents the key dates associated with a patent.
#[derive(Serialize, Debug)]
struct PatentDates {
    filing_date: String,
    publication_date: String,
    priority_date: String,
}

/// Represents the structure of a patent document.
/// This struct will be serialized into JSON format.
#[derive(Serialize, Debug)]
struct Patent {
    title: String,
    publication_number: String,
    abstract_text: String,
    inventors: Vec<String>,
    assignees: Vec<String>,
    dates: PatentDates,
    description: String,
    claims: Vec<String>,
}

/// A helper function to extract the text content from the first element matching a selector.
///
/// # Arguments
/// * `document` - A reference to the parsed HTML document.
/// * `selector_str` - The CSS selector to match.
///
/// # Returns
/// * `String` - The trimmed text content, or an empty string if not found.
fn extract_text(document: &Html, selector_str: &str) -> String {
    let selector = Selector::parse(selector_str).expect("Failed to parse selector");
    document
        .select(&selector)
        .next()
        .map_or(String::new(), |e| e.text().collect::<String>().trim().to_string())
}

/// A helper function to extract an attribute from the first element matching a selector.
///
/// # Arguments
/// * `document` - A reference to the parsed HTML document.
/// * `selector_str` - The CSS selector to match.
/// * `attr` - The name of the attribute to extract.
///
/// # Returns
/// * `String` - The attribute value, or an empty string if not found.
fn extract_attribute(document: &Html, selector_str: &str, attr: &str) -> String {
    let selector = Selector::parse(selector_str).expect("Failed to parse selector");
    document
        .select(&selector)
        .next()
        .and_then(|e| e.value().attr(attr))
        .map_or(String::new(), |s| s.to_string())
}

/// A helper function to extract all text content from elements matching a selector.
///
/// # Arguments
/// * `document` - A reference to the parsed HTML document.
/// * `selector_str` - The CSS selector to match.
///
/// # Returns
/// * `Vec<String>` - A vector of trimmed text content for each matched element.
fn extract_all_text(document: &Html, selector_str: &str) -> Vec<String> {
    let selector = Selector::parse(selector_str).expect("Failed to parse selector");
    document
        .select(&selector)
        .map(|e| e.text().collect::<String>().trim().to_string())
        .collect()
}

/// Parses an HTML document and extracts detailed patent information.
///
/// # Arguments
///
/// * `html_content` - A string slice that holds the HTML content of the patent page.
///
/// # Returns
///
/// * `Result<Patent, Box<dyn std::error::Error>>` - A `Result` containing the parsed `Patent` struct or an error.
fn parse_patent_html(html_content: &str) -> Result<Patent, Box<dyn std::error::Error>> {
    let document = Html::parse_document(html_content);

    // --- Basic Information ---
    let title = extract_text(&document, "span[itemprop='title']");
    let abstract_text = extract_text(&document, "div.abstract");
    let publication_number = extract_text(&document, "dd[itemprop='publicationNumber']");

    // --- People & Companies ---
    let inventors = extract_all_text(&document, "dd[itemprop='inventor']");
    let assignees = extract_all_text(&document, "dd[itemprop='assigneeCurrent']");

    // --- Dates ---
    let filing_date = extract_attribute(&document, "time[itemprop='filingDate']", "datetime");
    let publication_date = extract_attribute(&document, "time[itemprop='publicationDate']", "datetime");
    let priority_date = extract_attribute(&document, "time[itemprop='priorityDate']", "datetime");
    let dates = PatentDates {
        filing_date,
        publication_date,
        priority_date,
    };

    // --- Full Description ---
    // The description consists of multiple headings and paragraphs. We select all of them
    // and join them together with newlines to preserve some structure.
    let description_selector = Selector::parse("section[itemprop='description'] .description > *").unwrap();
    let description_parts: Vec<String> = document
        .select(&description_selector)
        .map(|element| {
            let text = element.text().collect::<String>().trim().to_string();
            // Add extra newlines after headings for better readability in the JSON output.
            if element.value().name() == "heading" {
                format!("\n\n{}\n", text)
            } else {
                text
            }
        })
        .collect();
    let description = description_parts.join("\n\n");

    // --- Claims ---
    let claim_selector = Selector::parse("div.claim-text").unwrap();
    let claims: Vec<String> = document
        .select(&claim_selector)
        .map(|e| e.text().collect::<String>().trim().to_string())
        .collect();

    Ok(Patent {
        title,
        publication_number,
        abstract_text,
        inventors,
        assignees,
        dates,
        description,
        claims,
    })
}

/// The main entry point for the asynchronous program.
#[tokio::main]
async fn main() -> Result<(), Box<dyn std::error::Error>> {
    let input_dir = "./data_google_patents/patents";
    let output_dir = "./data_google_patents/patents_json";

    // Create the output directory if it doesn't exist
    fs::create_dir_all(output_dir)?;

    let mut tasks = vec![];

    // Iterate over all HTML files in the input directory
    for entry in glob(&format!("{}/*.html", input_dir))? {
        let path = entry?;
        let output_path = Path::new(output_dir)
            .join(path.file_stem().unwrap())
            .with_extension("json");

        // Spawn a new asynchronous task for each file
        let task = task::spawn(async move {
            println!("Processing: {}", path.display());
            let html_content = fs::read_to_string(&path).unwrap();
            match parse_patent_html(&html_content) {
                Ok(patent) => {
                    let json_content = serde_json::to_string_pretty(&patent).unwrap();
                    fs::write(&output_path, json_content).unwrap();
                    println!("Successfully converted {} to {}", path.display(), output_path.display());
                }
                Err(e) => {
                    eprintln!("Failed to parse {}: {}", path.display(), e);
                }
            }
        });
        tasks.push(task);
    }

    // Wait for all tasks to complete
    for task in tasks {
        task.await?;
    }

    Ok(())
}
