import json
import csv

input_file = "02_metadata/ATP7B_orthologues_full.json"
output_file = "02_metadata/ATP7B_one2one_metadata.tsv"

with open(input_file, "r") as f:
    data = json.load(f)

homologies = data["data"][0]["homologies"]

one2one = [
    h for h in homologies
    if h["type"] == "ortholog_one2one"
]

fields = [
    "species",
    "taxon_id",
    "ensembl_gene_id",
    "ensembl_protein_id",
    "orthology_type",
    "percent_identity",
    "percent_positive"
]

with open(output_file, "w", newline="") as f:
    writer = csv.DictWriter(
        f,
        fieldnames=fields,
        delimiter="\t"
    )

    writer.writeheader()

    for h in one2one:
        target = h["target"]

        writer.writerow({
            "species": target["species"],
            "taxon_id": target["taxon_id"],
            "ensembl_gene_id": target["id"],
            "ensembl_protein_id": target["protein_id"],
            "orthology_type": h["type"],
            "percent_identity": target["perc_id"],
            "percent_positive": target["perc_pos"]
        })

print(f"Total one-to-one orthologs: {len(one2one)}")
print(f"Metadata written to: {output_file}")
