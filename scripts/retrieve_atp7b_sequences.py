import csv
import os
import time
import requests
from Bio.Seq import Seq

BASE_URL = "https://rest.ensembl.org"

METADATA_FILE = "02_metadata/ATP7B_one2one_metadata.tsv"
PROTEIN_DIR = "01_raw_sequences/protein"
CDS_DIR = "01_raw_sequences/cds"
QC_FILE = "02_metadata/ATP7B_sequence_QC.tsv"

os.makedirs(PROTEIN_DIR, exist_ok=True)
os.makedirs(CDS_DIR, exist_ok=True)

session = requests.Session()

JSON_HEADERS = {
    "Content-Type": "application/json",
    "Accept": "application/json"
}

FASTA_HEADERS = {
    "Content-Type": "text/plain",
    "Accept": "text/plain"
}


def get(url, headers, attempts=3):
    for attempt in range(1, attempts + 1):
        try:
            response = session.get(
                url,
                headers=headers,
                timeout=60
            )

            if response.status_code == 429:
                wait = int(response.headers.get("Retry-After", "5"))
                print(f"    Rate limited. Waiting {wait} seconds...")
                time.sleep(wait)
                continue

            response.raise_for_status()
            return response

        except requests.RequestException as e:
            if attempt == attempts:
                raise

            wait = attempt * 3
            print(f"    Request failed: {e}")
            print(f"    Retrying in {wait} seconds...")
            time.sleep(wait)

    raise RuntimeError("Request failed after retries")


def get_gene(gene_id):
    url = f"{BASE_URL}/lookup/id/{gene_id}?expand=1"
    return get(url, JSON_HEADERS).json()


def find_matching_transcript(gene_data, protein_id):

    for transcript in gene_data.get("Transcript", []):

        translation = transcript.get("Translation")

        if translation:
            translation_id = translation.get("id")

            if translation_id == protein_id:
                return transcript

    return None


def get_sequence(transcript_id, sequence_type):

    url = (
        f"{BASE_URL}/sequence/id/"
        f"{transcript_id}?type={sequence_type}"
    )

    return get(url, FASTA_HEADERS).text.strip()


def clean_sequence(fasta_text):

    return "".join(
        line.strip()
        for line in fasta_text.splitlines()
        if not line.startswith(">")
    ).upper()


with open(METADATA_FILE, newline="") as f:
    records = list(
        csv.DictReader(f, delimiter="\t")
    )


qc_results = []

print("=" * 70)
print(f"ATP7B SEQUENCE RETRIEVAL")
print(f"Orthologs to process: {len(records)}")
print("=" * 70)


for i, record in enumerate(records, start=1):

    species = record["species"]
    gene_id = record["ensembl_gene_id"]
    expected_protein_id = record["ensembl_protein_id"]

    print(f"\n[{i}/{len(records)}] {species}")
    print(f"    Gene:    {gene_id}")
    print(f"    Protein: {expected_protein_id}")

    result = {
        "species": species,
        "taxon_id": record["taxon_id"],
        "gene_id": gene_id,
        "expected_protein_id": expected_protein_id,
        "transcript_id": "",
        "protein_length": "",
        "cds_length": "",
        "cds_divisible_by_3": "",
        "translation_matches_protein": "",
        "internal_stop": "",
        "qc_status": "FAIL",
        "error": ""
    }

    try:

        # --------------------------------------------------
        # 1. Retrieve gene record
        # --------------------------------------------------

        gene_data = get_gene(gene_id)

        transcript = find_matching_transcript(
            gene_data,
            expected_protein_id
        )

        if transcript is None:
            raise ValueError(
                "Matching transcript for expected protein not found"
            )

        transcript_id = transcript["id"]
        result["transcript_id"] = transcript_id

        print(f"    Transcript: {transcript_id}")

        # --------------------------------------------------
        # 2. Retrieve protein + CDS
        # --------------------------------------------------

        protein_fasta = get_sequence(
            transcript_id,
            "protein"
        )

        cds_fasta = get_sequence(
            transcript_id,
            "cds"
        )

        protein = clean_sequence(protein_fasta)
        cds = clean_sequence(cds_fasta)

        result["protein_length"] = len(protein)
        result["cds_length"] = len(cds)

        print(f"    Protein length: {len(protein)} aa")
        print(f"    CDS length:     {len(cds)} nt")

        # --------------------------------------------------
        # 3. CDS length check
        # --------------------------------------------------

        divisible = len(cds) % 3 == 0
        result["cds_divisible_by_3"] = divisible

        if not divisible:
            raise ValueError(
                "CDS length is not divisible by 3"
            )

        # --------------------------------------------------
        # 4. Translate CDS
        # --------------------------------------------------

        translated = str(
            Seq(cds).translate(to_stop=False)
        )

        translated_no_stop = translated.rstrip("*")

        matches = translated_no_stop == protein

        result["translation_matches_protein"] = matches

        if not matches:
            raise ValueError(
                "Translated CDS does not match protein"
            )

        # --------------------------------------------------
        # 5. Internal stop check
        # --------------------------------------------------

        internal_stop = "*" in translated[:-1]

        result["internal_stop"] = internal_stop

        if internal_stop:
            raise ValueError(
                "Internal stop codon detected"
            )

        # --------------------------------------------------
        # 6. Save protein
        # --------------------------------------------------

        safe_species = species.replace(" ", "_")

        protein_file = (
            f"{PROTEIN_DIR}/"
            f"{safe_species}_ATP7B.fasta"
        )

        cds_file = (
            f"{CDS_DIR}/"
            f"{safe_species}_ATP7B_CDS.fasta"
        )

        with open(protein_file, "w") as f:
            f.write(
                f">{species}|{expected_protein_id}\n"
            )
            f.write(protein + "\n")

        with open(cds_file, "w") as f:
            f.write(
                f">{species}|{transcript_id}\n"
            )
            f.write(cds + "\n")

        result["qc_status"] = "PASS"

        print("    QC: PASS")

    except Exception as e:

        result["error"] = str(e)

        print(f"    QC: FAIL — {e}")

    qc_results.append(result)

    # Respect Ensembl REST service
    time.sleep(0.3)


# ------------------------------------------------------
# Write QC table
# ------------------------------------------------------

fields = [
    "species",
    "taxon_id",
    "gene_id",
    "expected_protein_id",
    "transcript_id",
    "protein_length",
    "cds_length",
    "cds_divisible_by_3",
    "translation_matches_protein",
    "internal_stop",
    "qc_status",
    "error"
]

with open(QC_FILE, "w", newline="") as f:

    writer = csv.DictWriter(
        f,
        fieldnames=fields,
        delimiter="\t"
    )

    writer.writeheader()
    writer.writerows(qc_results)


passed = sum(
    r["qc_status"] == "PASS"
    for r in qc_results
)

failed = len(qc_results) - passed

print()
print("=" * 70)
print("FINAL QC SUMMARY")
print("=" * 70)
print(f"Total orthologs: {len(qc_results)}")
print(f"Passed QC:       {passed}")
print(f"Failed QC:       {failed}")
print(f"QC table:        {QC_FILE}")
print("=" * 70)
