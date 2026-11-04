import csv
import os
import time
import requests
from Bio.Seq import Seq

BASE_URL = "https://rest.ensembl.org"
QC_FILE = "02_metadata/ATP7B_sequence_QC.tsv"
PROTEIN_DIR = "01_raw_sequences/protein"
CDS_DIR = "01_raw_sequences/cds"

os.makedirs(PROTEIN_DIR, exist_ok=True)
os.makedirs(CDS_DIR, exist_ok=True)

session = requests.Session()
session.headers.update({
    "Content-Type": "application/json",
    "User-Agent": "ATP7B-MSc-project/1.0"
})


def request_json(url, attempts=6):
    for attempt in range(1, attempts + 1):
        try:
            r = session.get(url, timeout=90)

            if r.status_code == 200:
                return r.json()

            if r.status_code in (429, 500, 502, 503, 504):
                wait = min(60, 5 * (2 ** (attempt - 1)))
                print(f"    HTTP {r.status_code}; retry {attempt}/{attempts}; waiting {wait}s")
                time.sleep(wait)
                continue

            r.raise_for_status()

        except requests.RequestException as e:
            if attempt == attempts:
                raise

            wait = min(60, 5 * (2 ** (attempt - 1)))
            print(f"    Request error; retry {attempt}/{attempts}; waiting {wait}s")
            time.sleep(wait)

    raise RuntimeError("Request failed after all retries")


def request_text(url, attempts=6):
    for attempt in range(1, attempts + 1):
        try:
            r = session.get(url, timeout=90)

            if r.status_code == 200:
                return r.text.strip()

            if r.status_code in (429, 500, 502, 503, 504):
                wait = min(60, 5 * (2 ** (attempt - 1)))
                print(f"    HTTP {r.status_code}; retry {attempt}/{attempts}; waiting {wait}s")
                time.sleep(wait)
                continue

            r.raise_for_status()

        except requests.RequestException:
            if attempt == attempts:
                raise

            wait = min(60, 5 * (2 ** (attempt - 1)))
            print(f"    Request error; retry {attempt}/{attempts}; waiting {wait}s")
            time.sleep(wait)

    raise RuntimeError("Request failed after all retries")


def clean_sequence(text):
    lines = text.splitlines()
    seq = "".join(
        line.strip()
        for line in lines
        if line.strip() and not line.startswith(">")
    )
    return seq.upper()


def write_fasta(path, header, sequence):
    with open(path, "w") as f:
        f.write(f">{header}\n")
        for i in range(0, len(sequence), 80):
            f.write(sequence[i:i+80] + "\n")


# Read the original QC table.
with open(QC_FILE, newline="") as f:
    rows = list(csv.DictReader(f, delimiter="\t"))

failed = [
    row for row in rows
    if row["qc_status"] == "FAIL"
    and "Translated CDS does not match protein" not in row["error"]
]

print("=" * 70)
print("ATP7B SEQUENCE RECOVERY")
print("=" * 70)
print(f"Transient/API failures to retry: {len(failed)}")
print()

recovered = 0
still_failed = 0

for i, row in enumerate(failed, 1):

    species = row["species"]
    gene_id = row["gene_id"]
    expected_protein = row["expected_protein_id"]

    print(f"[{i}/{len(failed)}] {species}")
    print(f"    Gene: {gene_id}")
    print(f"    Expected protein: {expected_protein}")

    try:
        # Get expanded gene annotation.
        lookup_url = f"{BASE_URL}/lookup/id/{gene_id}?expand=1"
        data = request_json(lookup_url)

        transcript_id = None

        for transcript in data.get("Transcript", []):
            translation = transcript.get("Translation")

            if translation and translation.get("id") == expected_protein:
                transcript_id = transcript.get("id")
                break

        if transcript_id is None:
            raise RuntimeError(
                f"Could not find transcript corresponding to {expected_protein}"
            )

        print(f"    Transcript: {transcript_id}")

        # Retrieve protein.
        protein_url = (
            f"{BASE_URL}/sequence/id/{transcript_id}?type=protein"
        )
        protein = clean_sequence(request_text(protein_url))

        # Retrieve CDS.
        cds_url = (
            f"{BASE_URL}/sequence/id/{transcript_id}?type=cds"
        )
        cds = clean_sequence(request_text(cds_url))

        if not protein:
            raise RuntimeError("Empty protein sequence")

        if not cds:
            raise RuntimeError("Empty CDS sequence")

        if len(cds) % 3 != 0:
            raise RuntimeError(
                f"CDS length {len(cds)} is not divisible by 3"
            )

        translated = str(Seq(cds).translate())

        # Remove terminal stop only for comparison.
        translated_no_stop = translated.rstrip("*")

        if translated_no_stop != protein.rstrip("*"):
            raise RuntimeError("Translated CDS does not match protein")

        # Check internal stop codons.
        if "*" in translated_no_stop:
            raise RuntimeError("Internal stop codon detected")

        protein_file = os.path.join(
            PROTEIN_DIR, f"{species}_ATP7B_protein.fasta"
        )
        cds_file = os.path.join(
            CDS_DIR, f"{species}_ATP7B_CDS.fasta"
        )

        write_fasta(
            protein_file,
            f"{species}|{gene_id}|{expected_protein}",
            protein
        )

        write_fasta(
            cds_file,
            f"{species}|{gene_id}|{transcript_id}",
            cds
        )

        recovered += 1
        print(
            f"    RECOVERED: protein={len(protein)} aa, "
            f"CDS={len(cds)} nt"
        )

    except Exception as e:
        still_failed += 1
        print(f"    FAILED: {e}")

    # Be polite to Ensembl.
    time.sleep(2)

print()
print("=" * 70)
print("RECOVERY SUMMARY")
print("=" * 70)
print(f"Recovered:       {recovered}")
print(f"Still failed:    {still_failed}")
print("=" * 70)
