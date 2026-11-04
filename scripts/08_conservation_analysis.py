from Bio import AlignIO
from collections import Counter
import pandas as pd
import matplotlib.pyplot as plt

alignment_file = "04_protein_alignment/ATP7B_72_filtered_aligned.fasta"

aln = AlignIO.read(alignment_file, "fasta")

human = next(
    r for r in aln
    if "P35670" in r.description
)

human_seq = str(human.seq)

results = []
human_position = 0

for i in range(aln.get_alignment_length()):

    residues = [str(r.seq)[i] for r in aln]
    nongap = [aa for aa in residues if aa != "-"]

    if human_seq[i] != "-":
        human_position += 1
        human_residue = human_seq[i]
    else:
        human_residue = "-"

    if not nongap:
        continue

    counts = Counter(nongap)
    most_common_residue, count = counts.most_common(1)[0]

    conservation = count / len(nongap) * 100
    gap_percentage = (len(residues) - len(nongap)) / len(residues) * 100

    results.append({
        "Alignment_Position": i + 1,
        "Human_Position": human_position if human_residue != "-" else None,
        "Human_Amino_Acid": human_residue,
        "Most_Common_Amino_Acid": most_common_residue,
        "Conservation_Percent": round(conservation, 2),
        "Sequences_With_Residue": len(nongap),
        "Gap_Percent": round(gap_percentage, 2)
    })

df = pd.DataFrame(results)

df.to_csv(
    "09_conservation/ATP7B_conservation.tsv",
    sep="\t",
    index=False
)

# Highly conserved positions with at least 70% sequence coverage
high = df[
    (df["Conservation_Percent"] >= 90) &
    (df["Sequences_With_Residue"] >= 0.70 * len(aln))
]

high.to_csv(
    "09_conservation/ATP7B_highly_conserved_positions.tsv",
    sep="\t",
    index=False
)

plt.figure(figsize=(14, 5))
plt.plot(
    df["Alignment_Position"],
    df["Conservation_Percent"]
)

plt.axhline(90, linestyle="--", label="90% conservation")
plt.xlabel("Alignment position")
plt.ylabel("Conservation (%)")
plt.title("ATP7B Protein Conservation Across 72 Species")
plt.ylim(0, 105)
plt.legend()
plt.tight_layout()

plt.savefig(
    "12_figures/ATP7B_conservation.png",
    dpi=300
)

print("Sequences:", len(aln))
print("Alignment positions analysed:", len(df))
print("Highly conserved positions:", len(high))
print("Mean conservation:", round(df["Conservation_Percent"].mean(), 2), "%")
print("Mean gap percentage:", round(df["Gap_Percent"].mean(), 2), "%")
print("\nTop 10 highly conserved positions:")
print(high.head(10).to_string(index=False))
print("\nConservation analysis completed.")
