use glob::glob;
use scraper::{Html, Selector};
use serde::Serialize;
use std::fs;
use std::path::Path;
use tokio::task;

/// Represents the structure of a patent document.
/// This struct will be serialized into JSON format.
#[derive(Serialize, Debug)]
struct Patent {
    title: String,
    abstract_text: String,
    publication_number: String,
    inventors: Vec<String>,
    assignees: Vec<String>,
    claims: Vec<String>,
}

/// Parses an HTML document and extracts patent information.
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

    // Selector for the patent title
    let title_selector = Selector::parse("span[itemprop='title']").unwrap();
    let title = document
        .select(&title_selector)
        .next()
        .map_or(String::new(), |e| e.text().collect::<String>().trim().to_string());

    // Selector for the abstract
    let abstract_selector = Selector::parse("div.abstract").unwrap();
    let abstract_text = document
        .select(&abstract_selector)
        .next()
        .map_or(String::new(), |e| e.text().collect::<String>().trim().to_string());

    // Selector for the publication number
    let publication_number_selector = Selector::parse("dd[itemprop='publicationNumber']").unwrap();
    let publication_number = document
        .select(&publication_number_selector)
        .next()
        .map_or(String::new(), |e| e.text().collect::<String>().trim().to_string());

    // Selector for inventors
    let inventor_selector = Selector::parse("dd[itemprop='inventor']").unwrap();
    let inventors: Vec<String> = document
        .select(&inventor_selector)
        .map(|e| e.text().collect::<String>().trim().to_string())
        .collect();

    // Selector for assignees
    let assignee_selector = Selector::parse("dd[itemprop='assigneeCurrent']").unwrap();
    let assignees: Vec<String> = document
        .select(&assignee_selector)
        .map(|e| e.text().collect::<String>().trim().to_string())
        .collect();

    // Selector for claims
    let claim_selector = Selector::parse("div.claim-text").unwrap();
    let claims: Vec<String> = document
        .select(&claim_selector)
        .map(|e| e.text().collect::<String>().trim().to_string())
        .collect();

    Ok(Patent {
        title,
        abstract_text,
        publication_number,
        inventors,
        assignees,
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
