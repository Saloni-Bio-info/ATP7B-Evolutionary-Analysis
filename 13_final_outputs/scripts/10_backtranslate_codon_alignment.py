from Bio import AlignIO, SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from pathlib import Path

protein_file = "04_protein_alignment/ATP7B_72_filtered_aligned.fasta"
cds_dir = Path("01_raw_sequences/cds")
output = "05_codon_alignment/ATP7B_71_codon_aligned.fasta"

protein_alignment = AlignIO.read(protein_file, "fasta")

codon_records = []
errors = []

for record in protein_alignment:

    species = record.id.split("|")[0]

    # Human reference excluded: no human CDS currently available
    if species == "sp":
        continue

    cds_file = cds_dir / f"{species}_ATP7B_CDS.fasta"

    if not cds_file.exists():
        errors.append(f"{species}: CDS file missing")
        continue

    protein = str(record.seq).upper()
    cds_record = SeqIO.read(cds_file, "fasta")
    cds = str(cds_record.seq).upper()

    if len(cds) % 3 != 0:
        errors.append(f"{species}: CDS length not divisible by 3")
        continue

    translated = str(Seq(cds).translate()).upper()

    # Remove terminal stop codon if present
    if translated.endswith("*"):
        translated = translated[:-1]
        cds = cds[:-3]

    if translated != protein.replace("-", ""):
        errors.append(f"{species}: protein-CDS translation mismatch")
        continue

    codons = []
    position = 0

    for aa in protein:
        if aa == "-":
            codons.append("---")
        else:
            codon = cds[position:position+3]
            codons.append(codon)
            position += 3

    if position != len(cds):
        errors.append(f"{species}: unused CDS nucleotides")
        continue

    codon_records.append(
        SeqRecord(
            Seq("".join(codons)),
            id=record.id,
            description=""
        )
    )

if errors:
    print("ERRORS:")
    for error in errors:
        print(error)
    raise SystemExit("Codon alignment was not saved.")

SeqIO.write(codon_records, output, "fasta")

print("Codon alignment completed.")
print("Sequences:", len(codon_records))
print("Codon alignment length:", len(codon_records[0].seq) // 3)
print("Nucleotide alignment length:", len(codon_records[0].seq))
print("Output:", output)
