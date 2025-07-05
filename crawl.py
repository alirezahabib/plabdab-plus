#!/usr/bin/env python3
import os
import re
import json
import time
import argparse
import pandas as pd
import xmltodict
from Bio import Entrez

# ——— Defaults ———
DATA_CSV    = "data.csv"
OUTDIR      = "data_ncbi"
BASE_LIMIT  = 0.34    # ~3 req/sec without an API key
API_LIMIT   = 0.1     # ~10 req/sec with an API key

# pre-compile a regex that matches any control char except newline, carriage return, tab
_CTRLS = re.compile(r'[\x00-\x08\x0B\x0C\x0E-\x1F\x7F]')

def _sanitize(raw_text: str) -> str:
    """Strip out unprintable control chars that break json.loads."""
    return _CTRLS.sub('', raw_text)

def fetch_json(func, **kwargs):
    """
    Call an Entrez JSON endpoint, sanitize out bad chars, parse, close, rate-limit.
    Falls back to a second sanitization pass if needed.
    """
    handle = func(**kwargs)
    raw = handle.read()
    handle.close()
    time.sleep(fetch_json.delay)

    try:
        # first try strict=False
        return json.loads(raw, strict=False)
    except json.JSONDecodeError:
        # second, strip out any remaining control chars
        cleaned = _sanitize(raw)
        return json.loads(cleaned, strict=False)

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

    # Build a lookup of already-crawled accession.version → pubmed_present
    already_crawled = {}
    for fname in os.listdir(args.outdir):
        if fname.endswith(".json"):
            acc_ver = fname[:-5]
            pubmed_present = False
            try:
                with open(os.path.join(args.outdir, fname)) as fp:
                    j = json.load(fp)
                    pubmed_present = bool(j.get("pubmed"))
            except Exception:
                # corrupt json or other issues – treat as not having pubmed
                pass
            already_crawled[acc_ver] = pubmed_present

    success_count = 0

    for idx, row in df.iterrows():
        url     = str(row.get("url", ""))
        # Skip rows without an NCBI protein link
        if not url.startswith("https://www.ncbi.nlm.nih.gov/protein/"):
            continue

        # Extract accession.version string from URL (part after last '/')
        acc_ver = url.rstrip("/").split("/")[-1].split("?")[0]

        # If we've already processed this accession.version, just update crawl col
        if acc_ver in already_crawled:
            flags = ["ncbi"]
            if already_crawled[acc_ver]:
                flags.append("pubmed")
            df.at[idx, "crawl"] = ";".join(flags)  # type: ignore[assignment]
            continue

        allraw = {"accession_version": acc_ver}
        got_pubmed = False

        try:
            # -------------------- UID lookup --------------------
            search_term = f"{acc_ver.replace('_', '.')}[Accession]"
            handle = Entrez.esearch(db="protein", term=search_term, retmode="xml")
            rec = Entrez.read(handle)
            handle.close()
            if isinstance(rec, dict):
                uid_list = rec.get("IdList", [])  # type: ignore[attr-defined]
            else:
                uid_list = []

            if not uid_list:
                # No valid uids – don't tag with ncbi
                df.at[idx, "crawl"] = ""  # type: ignore[assignment]
                continue

            allraw["uid_list"] = uid_list

            # -------------------- ESummary --------------------
            allraw["esummary"] = fetch_json(
                Entrez.esummary,
                db="protein",
                id=",".join(uid_list),
                retmode="json"
            )

            # -------------------- ELink (PubMed) --------------------
            elink_json = fetch_json(
                Entrez.elink,
                dbfrom="protein",
                db="pubmed",
                id=",".join(uid_list),
                retmode="json"
            )
            allraw["elink"] = elink_json

            pmids = []
            for linkset in elink_json.get("linksets", []):  # type: ignore[attr-defined]
                for ldb in linkset.get("linksetdbs", []):  # type: ignore[attr-defined]
                    if ldb.get("dbto") == "pubmed":
                        pmids.extend(ldb.get("links", []))

            if pmids:
                got_pubmed = True
                allraw["pubmed"] = {}
                for pmid in set(pmids):
                    try:
                        xml_text = fetch_text(
                            Entrez.efetch,
                            db="pubmed",
                            id=pmid,
                            rettype="xml",
                            retmode="xml"
                        )
                        parsed = xmltodict.parse(xml_text)
                        allraw["pubmed"][pmid] = parsed
                    except Exception as e:
                        print(f"ERROR fetching PubMed {pmid} for {acc_ver}: {e}")

            # -------------------- Write JSON --------------------
            outpath = os.path.join(args.outdir, f"{acc_ver}.json")
            with open(outpath, "w") as f:
                json.dump(allraw, f, indent=2)

            # Cache so we don't duplicate work in this run
            already_crawled[acc_ver] = got_pubmed

            # Update crawl column
            flags = ["ncbi"]
            if got_pubmed:
                flags.append("pubmed")
            df.at[idx, "crawl"] = ";".join(flags)  # type: ignore[assignment]

            success_count += 1
            if success_count % 1000 == 0:
                df.to_csv(args.csv, index=False)
                print(f"[{success_count}] crawled — checkpoint saved.")

        except Exception as e:
            print(f"ERROR processing {acc_ver}: {e}")

    # Final save
    df.to_csv(args.csv, index=False)
    print(f"Done. Total new crawls: {success_count}")

if __name__ == "__main__":
    main()
