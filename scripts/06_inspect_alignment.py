from Bio import AlignIO

alignment = AlignIO.read(
    "data/processed/coi_aligned.fasta",
    "fasta"
)

print("Species | Leading gaps | Trailing gaps | Internal gaps")

for record in alignment:
    sequence = str(record.seq)

    first_base = next(
        (i for i, base in enumerate(sequence) if base != "-"),
        len(sequence)
    )

    last_base = next(
        (i for i, base in enumerate(reversed(sequence)) if base != "-"),
        len(sequence)
    )

    leading = sequence[:first_base].count("-")
    trailing = sequence[len(sequence) - last_base:].count("-") if last_base else 0
    internal = sequence.count("-") - leading - trailing

    print(
        record.id,
        "|", leading,
        "|", trailing,
        "|", internal
    )
