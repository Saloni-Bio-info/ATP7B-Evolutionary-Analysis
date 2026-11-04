import pandas as pd

df = pd.read_csv(
    "05_conservation/ATP7B_human_mapped_conservation.tsv",
    sep="\t"
)

# UniProt P35670 domain annotations
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

    region = df[
        (df["Human_ATP7B_Position"] >= start) &
        (df["Human_ATP7B_Position"] <= end)
    ].copy()

    # Apply the same occupancy threshold
    reliable = region[region["Occupancy_Percent"] >= 70]

    results.append({
        "Domain": name,
        "Human_Start": start,
        "Human_End": end,
        "Total_Positions": len(region),
        "Positions_At_Least_70_Occupancy": len(reliable),
        "Mean_Conservation_Percent": round(
            reliable["Residue_Conservation_Percent"].mean(), 2
        ),
        "Highly_Conserved_Positions_GE90": sum(
            reliable["Residue_Conservation_Percent"] >= 90
        ),
        "Fully_Conserved_Positions": sum(
            (reliable["Residue_Conservation_Percent"] == 100) &
            (reliable["Occupancy_Percent"] == 100)
        )
    })

result_df = pd.DataFrame(results)

result_df.to_csv(
    "05_conservation/ATP7B_HMA_domain_conservation.tsv",
    sep="\t",
    index=False
)

print("\nHMA DOMAIN CONSERVATION")
print(result_df.to_string(index=False))
