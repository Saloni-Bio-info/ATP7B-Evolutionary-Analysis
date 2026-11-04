from Bio import AlignIO
import pandas as pd

alignment = AlignIO.read(
    "04_protein_alignment/ATP7B_72_filtered_aligned.fasta",
    "fasta"
)

# Identify human reference sequence
human_record = None

for record in alignment:
    if "P35670" in record.id:
        human_record = record
        break

if human_record is None:
    raise ValueError("Human ATP7B reference P35670 not found!")

human_sequence = str(human_record.seq)

# Load conservation results
df = pd.read_csv(
    "05_conservation/ATP7B_position_conservation.tsv",
    sep="\t"
)

human_position = 0
mapping = []

for i, residue in enumerate(human_sequence):

    if residue not in "-.":
        human_position += 1
        mapping.append((i + 1, human_position, residue))
    else:
        mapping.append((i + 1, None, "-"))

mapping_df = pd.DataFrame(
    mapping,
    columns=[
        "Alignment_Position",
        "Human_ATP7B_Position",
        "Human_Residue"
    ]
)

# Merge with conservation data
merged = df.merge(
    mapping_df,
    on="Alignment_Position"
)

merged.to_csv(
    "05_conservation/ATP7B_human_mapped_conservation.tsv",
    sep="\t",
    index=False
)

# Select highly conserved positions
highly_conserved = merged[
    (merged["Occupancy_Percent"] >= 70) &
    (merged["Residue_Conservation_Percent"] >= 90) &
    (merged["Human_Residue"] != "-")
]

highly_conserved.to_csv(
    "05_conservation/ATP7B_highly_conserved_human_positions.tsv",
    sep="\t",
    index=False
)

print("Human reference:", human_record.id)
print("Human protein length:", human_position)
print("Mapped alignment positions:", len(mapping_df))
print("Highly conserved positions mapped to human:", len(highly_conserved))

print("\nFirst 20 mapped highly conserved positions:")
print(
    highly_conserved[
        [
            "Alignment_Position",
            "Human_ATP7B_Position",
            "Human_Residue",
            "Occupancy_Percent",
            "Residue_Conservation_Percent"
        ]
    ].head(20).to_string(index=False)
)
