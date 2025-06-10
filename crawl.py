#!/usr/bin/env python3
import os
import json
import time
import argparse
import pandas as pd
import xmltodict
from Bio import Entrez

# ——— Defaults ———
DATA_CSV    = "data.csv"
OUTDIR      = "data"
BASE_LIMIT  = 0.34    # ~3 req/sec without an API key
API_LIMIT   = 0.1     # ~10 req/sec with an API key

def fetch_json(func, **kwargs):
    handle = func(**kwargs)
    obj = json.load(handle)
    handle.close()
    time.sleep(fetch_json.delay)
    return obj

def fetch_text(func, **kwargs):
    handle = func(**kwargs)
    txt = handle.read()
    handle.close()
    time.sleep(fetch_text.delay)
    return txt

def main():
    p = argparse.ArgumentParser(
        description="Download raw NCBI API data for proteins + linked PubMed.")
    p.add_argument("--email", required=True,
                   help="Your email (required by NCBI Entrez)")
    p.add_argument("--api_key", default=None,
                   help="Your NCBI API key (optional, raises rate limit)")
    p.add_argument("--outdir", default=OUTDIR,
                   help="Directory to save JSON files")
    p.add_argument("--csv", default=DATA_CSV,
                   help="Input CSV filename")
    args = p.parse_args()

    # — Entrez setup —
    Entrez.email = args.email
    Entrez.tool  = "protein_batch_crawler"
    if args.api_key:
        Entrez.api_key = args.api_key
        delay = API_LIMIT
        print(f"Using API key, rate limit set to ~10 req/sec.")
    else:
        delay = BASE_LIMIT
        print(f"No API key, rate limit set to ~3 req/sec.")

    fetch_json.delay = delay
    fetch_text.delay = delay

    # 1) Load CSV (including that unnamed first column)
    df = pd.read_csv(args.csv, dtype=str)
    firstcol = df.columns[0]

    # 2) Ensure 'crawl' column exists
    if "crawl" not in df.columns:
        df["crawl"] = ""
    else:
        df["crawl"] = df["crawl"].fillna("")

    os.makedirs(args.outdir, exist_ok=True)
    success_count = 0

    for idx, row in df.iterrows():
        url     = str(row.get("url",""))
        crawled = str(row["crawl"])
        if url.startswith("https://www.ncbi.nlm.nih.gov/protein/") and "ncbi" not in crawled:
            key     = row[firstcol]
            prot_id = row["ID"]

            allraw = {}
            got_pubmed = False
            try:
                # ESummary → JSON
                allraw["esummary"] = fetch_json(
                    Entrez.esummary,
                    db="protein",
                    id=prot_id,
                    retmode="json"
                )

                # ELink → JSON (to find PubMed IDs)
                elink_json = fetch_json(
                    Entrez.elink,
                    dbfrom="protein",
                    db="pubmed",
                    id=prot_id,
                    retmode="json"
                )
                allraw["elink"] = elink_json

                # Extract PMIDs
                pmids = []
                for linkset in elink_json.get("linksets", []):
                    for ldb in linkset.get("linksetdbs", []):
                        if ldb.get("dbto") == "pubmed":
                            pmids.extend(ldb.get("links", []))

                # Fetch & convert PubMed to JSON if present
                if pmids:
                    got_pubmed = True
                    allraw["pubmed"] = {}
                    for pmid in pmids:
                        xml_text = fetch_text(
                            Entrez.efetch,
                            db="pubmed",
                            id=pmid,
                            rettype="xml",
                            retmode="xml"
                        )
                        parsed = xmltodict.parse(xml_text)
                        allraw["pubmed"][pmid] = parsed

                # Write out data/{key}.json
                outpath = os.path.join(args.outdir, f"{key}.json")
                with open(outpath, "w") as f:
                    json.dump(allraw, f, indent=2)

                # Update crawl flags
                flags = ["ncbi"]
                if got_pubmed:
                    flags.append("pubmed")
                df.at[idx, "crawl"] = ",".join(flags)

                success_count += 1
                if success_count % 1000 == 0:
                    df.to_csv(args.csv, index=False)
                    print(f"[{success_count}] crawled — checkpoint saved.")

            except Exception as e:
                print(f"ERROR fetching {prot_id} ({key}): {e}")
            break

    # Final save
    df.to_csv(args.csv, index=False)
    print(f"Done. Total new crawls: {success_count}")

if __name__ == "__main__":
    main()
