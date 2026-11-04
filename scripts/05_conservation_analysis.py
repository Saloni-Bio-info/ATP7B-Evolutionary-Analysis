from Bio import AlignIO
from collections import Counter
import pandas as pd

# Load alignment
alignment = AlignIO.read(
    "04_protein_alignment/ATP7B_72_filtered_aligned.fasta",
    "fasta"
)

total_sequences = len(alignment)
alignment_length = alignment.get_alignment_length()

results = []

for i in range(alignment_length):

    column = alignment[:, i]

    residues = [aa for aa in column if aa not in "-."]
    counts = Counter(residues)

    non_gap = len(residues)
    gaps = total_sequences - non_gap

    if non_gap > 0:
        most_common_residue, most_common_count = counts.most_common(1)[0]
        residue_conservation = (most_common_count / non_gap) * 100
    else:
        most_common_residue = "-"
        most_common_count = 0
        residue_conservation = 0

    occupancy = (non_gap / total_sequences) * 100

    results.append({
        "Alignment_Position": i + 1,
        "Most_Common_Residue": most_common_residue,
        "Most_Common_Count": most_common_count,
        "Non_Gap_Sequences": non_gap,
        "Gap_Count": gaps,
        "Occupancy_Percent": round(occupancy, 2),
        "Residue_Conservation_Percent": round(residue_conservation, 2)
    })

df = pd.DataFrame(results)

# Save all positions
df.to_csv(
    "05_conservation/ATP7B_position_conservation.tsv",
    sep="\t",
    index=False
)

# Consider positions with at least 70% occupancy
reliable = df[df["Occupancy_Percent"] >= 70].copy()

# Summary statistics
fully_conserved = len(
    reliable[
        (reliable["Residue_Conservation_Percent"] == 100) &
        (reliable["Occupancy_Percent"] == 100)
    ]
)

highly_conserved = len(
    reliable[reliable["Residue_Conservation_Percent"] >= 90]
)

variable = len(
    reliable[reliable["Residue_Conservation_Percent"] < 90]
)

print("Total sequences:", total_sequences)
print("Alignment length:", alignment_length)
print("Positions with >=70% occupancy:", len(reliable))
print("Positions with <70% occupancy:", alignment_length - len(reliable))
print("Fully conserved positions across all sequences:", fully_conserved)
print("Highly conserved positions (>=90% identity):", highly_conserved)
print("Variable positions (<90% identity):", variable)

print("\nTop 20 conserved positions (>=70% occupancy):")

print(
    reliable.sort_values(
        ["Residue_Conservation_Percent", "Occupancy_Percent"],
        ascending=False
    ).head(20).to_string(index=False)
)
