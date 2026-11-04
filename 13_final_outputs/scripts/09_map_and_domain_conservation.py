import pandas as pd
from pathlib import Path

source = "09_conservation/ATP7B_conservation.tsv"

df = pd.read_csv(source, sep="\t")

# Map alignment positions to human ATP7B coordinates
df = df.rename(columns={
    "Human_Position": "Human_ATP7B_Position",
    "Human_Amino_Acid": "Human_Residue",
    "Conservation_Percent": "Residue_Conservation_Percent"
})

# Retain positions represented in human ATP7B
mapped = df[df["Human_ATP7B_Position"].notna()].copy()

mapped.to_csv(
    "09_conservation/ATP7B_human_mapped_conservation.tsv",
    sep="\t", index=False
)

# Identify highly conserved human positions
high = mapped[
    (mapped["Residue_Conservation_Percent"] >= 90) &
    (mapped["Sequences_With_Residue"] >= 0.70 * 72)
].copy()

high.to_csv(
    "09_conservation/ATP7B_highly_conserved_human_positions.tsv",
    sep="\t", index=False
)

# HMA domain coordinates from the existing project annotation
domains = {
    "HMA1": (58, 124),
    "HMA2": (143, 209),
    "HMA3": (257, 323),
    "HMA4": (359, 425),
    "HMA5": (488, 554),
    "HMA6": (564, 630)
}

results = []

for name, (start, end) in domains.items():

    region = mapped[
        mapped["Human_ATP7B_Position"].between(start, end)
    ]

    reliable = region[
        region["Sequences_With_Residue"] >= 0.70 * 72
    ]

    results.append({
        "Domain": name,
        "Human_Start": start,
        "Human_End": end,
        "Positions_Analysed": len(reliable),
        "Mean_Conservation_Percent": round(
            reliable["Residue_Conservation_Percent"].mean(), 2
        ),
        "Highly_Conserved_Positions_GE90": sum(
            reliable["Residue_Conservation_Percent"] >= 90
        ),
        "Fully_Conserved_Positions": sum(
            (reliable["Residue_Conservation_Percent"] == 100) &
            (reliable["Sequences_With_Residue"] == 72)
        )
    })

result = pd.DataFrame(results)

result.to_csv(
    "09_conservation/ATP7B_HMA_domain_conservation.tsv",
    sep="\t", index=False
)

print("\nHUMAN ATP7B MAPPING")
print("Mapped human positions:", len(mapped))
print("Highly conserved human positions:", len(high))

print("\nHMA DOMAIN CONSERVATION")
print(result.to_string(index=False))

print("\nAll output files saved successfully.")
