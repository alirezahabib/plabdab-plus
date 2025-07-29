#!/usr/bin/env python3
"""Fetch full-text articles for PubMed records referenced in data_ncbi JSON dumps.

This script walks over every .json file previously created by crawl.py under
./data_ncbi.  Each json contains a top-level key ``pubmed`` whose value is a
mapping of PMID → parsed XML (produced by ``xmltodict``).

For every PMID we attempt the following, in order:

1. Use Entrez.elink to see if the article is in PubMed Central (PMC).  If a
   matching PMC id is found, we fetch the NXML full-text with
   ``Entrez.efetch(db="pmc", rettype="full", retmode="text")``.
2. If no PMC version is available, fall back to ``Entrez.efetch`` from the
   PubMed database (rettype="medline", retmode="text") which at least returns
   the abstract.

The fetched text is written to ``<outdir>/<identifier>/fulltext.txt`` where
``identifier`` is, in order of preference: DOI (sanitised so it is a safe
path component) or the PMID/PMCID.

A small rate-limit is applied so we play nicely with NCBI servers.  Supply an
NCBI API key to increase the throughput.
"""
from __future__ import annotations

import argparse
import json
import re
import sys
import time
from pathlib import Path
from typing import Optional

import xmltodict
from Bio import Entrez

# -------- Rate limits (matching crawl.py) --------
BASE_LIMIT = 0.34   # ~3 req/sec without key
API_LIMIT = 0.06    # ~20 req/sec with key

# -------- Helpers --------
_ctrls = re.compile(r"[\x00-\x08\x0B\x0C\x0E-\x1F\x7F]")

def _sanitize_control_chars(raw: str) -> str:
    """Remove control characters that break parsers."""
    return _ctrls.sub("", raw)

def _sanitize_filename(s: str) -> str:
    """Make a string safe to use as a single path component."""
    return re.sub(r"[^A-Za-z0-9._-]", "_", s)

def _extract_doi(pubmed_xml: dict) -> Optional[str]:
    """Navigate the xmltodict structure to pull out the DOI if present."""
    try:
        article_ids = (
            pubmed_xml["PubmedArticleSet"]["PubmedArticle"]["MedlineCitation"]["Article"][
                "ArticleIdList"
            ]["ArticleId"]
        )
    except Exception:
        return None

    # ArticleId can be dict or list[dict]
    if isinstance(article_ids, dict):
        article_ids = [article_ids]

    for aid in article_ids:
        if aid.get("@IdType") == "doi":
            return str(aid.get("#text", "")).strip()
    return None

def _efetch_text(**kwargs) -> str:
    """Wrapper around Entrez.efetch that returns decoded string and respects delay."""
    handle = Entrez.efetch(**kwargs)
    raw: str = handle.read()
    handle.close()
    time.sleep(_efetch_text.delay)  # type: ignore[attr-defined]
    return raw

# Will be patched at runtime
_efetch_text.delay = BASE_LIMIT  # type: ignore[attr-defined]

# -------- Full-text retrieval --------

def fetch_full_text(pmid: str) -> tuple[Optional[str], str]:
    """Attempt to retrieve full-text for *pmid*.

    Returns (identifier, text). *identifier* is DOI, PMC id, or PMID.
    If no text is obtained, returns (None, "").
    """
    # 1) Try to find PMC id linked to this PMID
    try:
        handle = Entrez.elink(dbfrom="pubmed", db="pmc", id=pmid, retmode="json")
        raw = handle.read()
        handle.close()
        time.sleep(_efetch_text.delay)
        elink_json = json.loads(raw)
    except Exception as exc:
        print(f"  - Failed elink PMC lookup for PMID {pmid}: {exc}", file=sys.stderr)
        elink_json = None

    pmcid: Optional[str] = None
    if elink_json:
        try:
            linksets = elink_json.get("linksets", [])
            if linksets:
                ldbs = linksets[0].get("linksetdbs", [])
                for ldb in ldbs:
                    if ldb.get("dbto") == "pmc":
                        links = ldb.get("links", [])
                        if links:
                            pmcid = links[0]
                            break
        except Exception:
            pass

    # Fetch from PMC if possible
    if pmcid:
        try:
            text = _efetch_text(db="pmc", id=pmcid, rettype="full", retmode="text")
            if text and len(text.strip()) > 500:  # heuristic for "real" full-text
                return pmcid, text
        except Exception as exc:
            print(f"  - PMC fetch failed for {pmcid}: {exc}", file=sys.stderr)

    # Fall back to PubMed MEDLINE (abstract at minimum)
    try:
        text = _efetch_text(db="pubmed", id=pmid, rettype="medline", retmode="text")
        if text:
            return pmid, text
    except Exception as exc:
        print(f"  - PubMed fetch failed for PMID {pmid}: {exc}", file=sys.stderr)

    return None, ""

# -------- Main --------

def main() -> None:
    p = argparse.ArgumentParser(description="Fetch full-text for PubMed records present in data_ncbi json files.")
    p.add_argument("--email", required=True, help="NCBI requires an email address for Entrez usage.")
    p.add_argument("--api_key", help="NCBI API key (optional). Increases rate limit.")
    p.add_argument("--indir", default="data_ncbi", help="Directory containing json files from crawl.py")
    p.add_argument("--outdir", default="full_texts", help="Directory to place fetched full-texts")

    args = p.parse_args()

    # Entrez configuration
    Entrez.email = args.email
    Entrez.tool = "protein_pubmed_fulltext_crawler"
    if args.api_key:
        Entrez.api_key = args.api_key
        delay = API_LIMIT
        print("Using API key. Rate limit ~10 req/sec.")
    else:
        delay = BASE_LIMIT
        print("No API key. Rate limit ~3 req/sec.")

    _efetch_text.delay = delay  # type: ignore[attr-defined]

    indir = Path(args.indir)
    outdir = Path(args.outdir)
    if not indir.is_dir():
        print(f"Input directory '{indir}' does not exist.", file=sys.stderr)
        sys.exit(1)
    outdir.mkdir(parents=True, exist_ok=True)

    json_files = sorted(indir.glob("*.json"))
    print(f"Found {len(json_files)} JSON files to process.")

    processed = 0
    saved = 0

    for jpath in json_files:
        processed += 1
        try:
            with open(jpath) as fp:
                j = json.load(fp)
        except Exception as exc:
            print(f"[WARN] Skipping {jpath.name}: cannot parse JSON ({exc})")
            continue

        pubmed_map = j.get("pubmed")
        if not pubmed_map:
            continue

        for pmid, xml_dict in pubmed_map.items():
            doi = _extract_doi(xml_dict)
            identifier_source, text = fetch_full_text(pmid)
            if not text:
                continue

            identifier = doi or identifier_source or f"PMID_{pmid}"
            safe_id = _sanitize_filename(identifier)
            art_dir = outdir / safe_id
            art_dir.mkdir(parents=True, exist_ok=True)
            outfile = art_dir / "fulltext.txt"
            if outfile.exists():
                # Already have it; skip
                continue
            try:
                with open(outfile, "w", encoding="utf-8") as out_fp:
                    out_fp.write(text)
                saved += 1
            except Exception as exc:
                print(f"[ERROR] Failed to write full-text for {safe_id}: {exc}")

        if processed % 500 == 0:
            print(f"Processed {processed}/{len(json_files)} files, saved {saved} full-texts.")

    print(f"Done. Total files processed: {processed}. Full-texts saved: {saved}.")

if __name__ == "__main__":
    main() 