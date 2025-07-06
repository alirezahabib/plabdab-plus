from Bio import Entrez
from time import sleep

Entrez.email = "your_email@example.com"

def accession_to_uid(acc_ver):
    """Convert an accession.version string (e.g. 'AF416909_1') into the NCBI UID."""
    # NCBI expects accession.version with a dot, not underscore
    acc_ver = acc_ver.replace("_", ".")
    handle = Entrez.esearch(
        db="protein",
        term=f"{acc_ver}[Accession]",
        retmode="xml"
    )
    rec = Entrez.read(handle)
    handle.close()
    idlist = rec.get("IdList", [])
    if not idlist:
        raise ValueError(f"No UID found for accession {acc_ver}"))
    return idlist[0]

def fetch_summaries(accession_ids):
    """
    Given a list of accession.version strings (underscore or dot),
    look up their UIDs and return the esummary record list.
    """
    # Step 1: map accessions → UIDs
    uids = []
    for acc in accession_ids:
        try:
            uid = accession_to_uid(acc)
            uids.append(uid)
        except Exception as e:
            print(f"Warning: {e}")
        # be kind to NCBI
        sleep(0.4)

    if not uids:
        return []

    # Step 2: batch‐fetch summaries
    id_str = ",".join(uids)
    with Entrez.esummary(db="protein", id=id_str, retmode="xml") as handle:
        summaries = Entrez.read(handle)
    return summaries

if __name__ == "__main__":
    accessions = ["AF416909_1", "AAB18787"]
    summaries = fetch_summaries(accessions)
    for info in summaries:
        print(f"> {info['Title']}  (Length: {info['Length']} aa, UID: {info['Id']})")

