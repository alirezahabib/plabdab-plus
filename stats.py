import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import numpy as np

def save_plot_to_pdf(fig, pdf, title):
    """Saves the given figure to the PDF."""
    fig.suptitle(title, fontsize=16)
    pdf.savefig(fig, bbox_inches='tight')

def main():
    """
    Main function to read data, calculate statistics, and generate plots.
    """
    try:
        df = pd.read_csv('./data.csv')
    except FileNotFoundError:
        print("Error: data.csv not found. Please make sure the file is in the root directory.")
        return

    # Ensure the 'url' column exists
    if 'url' not in df.columns:
        print("Error: 'url' column not found in data.csv.")
        return

    # Define URL prefixes
    url_prefixes = {
        'Google Patents': 'https://patents.google.com',
        'NCBI Protein': 'https://www.ncbi.nlm.nih.gov/protein',
        'SAbDab Viewer': 'https://opig.stats.ox.ac.uk/webapps/newsabdab/sabdab/structureviewer/?pdb=',
        'TheraSAbDab': 'https://opig.stats.ox.ac.uk/webapps/newsabdab/therasabdab/search/?therapeutic='
    }

    stats = {}

    # Calculate stats for each prefix
    for name, prefix in url_prefixes.items():
        # Count rows with the prefix
        total_count = df['url'].str.startswith(prefix, na=False).sum()
        # Filter rows with the prefix
        prefixed_urls = df.loc[df['url'].str.startswith(prefix, na=False), 'url']
        # Count unique URLs
        unique_count = prefixed_urls.nunique()
        stats[name] = {'total': total_count, 'unique': unique_count}

    # Calculate "other" URLs
    all_prefixed_mask = pd.Series([False] * len(df))
    for prefix in url_prefixes.values():
        all_prefixed_mask |= df['url'].str.startswith(prefix, na=False)

    other_urls = df.loc[~all_prefixed_mask, 'url']
    stats['Other'] = {'total': len(other_urls), 'unique': other_urls.nunique()}

    # Prepare data for plotting
    labels = list(stats.keys())
    total_counts = [stats[label]['total'] for label in labels]
    unique_counts = [stats[label]['unique'] for label in labels]

    # Print statistics to console
    print("URL Statistics:")
    for label in labels:
        print(f"  {label}:")
        print(f"    Total rows: {stats[label]['total']}")
        print(f"    Unique URLs: {stats[label]['unique']}")

    # Create plots
    with PdfPages('url_statistics.pdf') as pdf:
        # Plot for total counts
        fig1, ax1 = plt.subplots(figsize=(10, 7))
        bars1 = ax1.bar(labels, total_counts, color='skyblue')
        ax1.set_ylabel('Total Count')
        ax1.set_title('Total URL Counts by Category')
        ax1.set_xticklabels(labels, rotation=45, ha='right')
        for bar in bars1:
            yval = bar.get_height()
            ax1.text(bar.get_x() + bar.get_width()/2.0, yval, int(yval), va='bottom') # va: vertical alignment
        plt.tight_layout(rect=[0, 0.03, 1, 0.95])
        save_plot_to_pdf(fig1, pdf, 'Total URL Counts')

        # Plot for unique counts
        fig2, ax2 = plt.subplots(figsize=(10, 7))
        bars2 = ax2.bar(labels, unique_counts, color='lightgreen')
        ax2.set_ylabel('Unique Count')
        ax2.set_title('Unique URL Counts by Category')
        ax2.set_xticklabels(labels, rotation=45, ha='right')
        for bar in bars2:
            yval = bar.get_height()
            ax2.text(bar.get_x() + bar.get_width()/2.0, yval, int(yval), va='bottom')
        plt.tight_layout(rect=[0, 0.03, 1, 0.95])
        save_plot_to_pdf(fig2, pdf, 'Unique URL Counts')

    print("\nGraphs saved to url_statistics.pdf")

if __name__ == '__main__':
    main()
