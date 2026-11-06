from Bio import AlignIO

alignment = AlignIO.read(
    "data/processed/coi_aligned.fasta",
    "fasta"
)

sequences = [str(record.seq).upper() for record in alignment]

shared_positions = []

for position, column in enumerate(zip(*sequences), start=1):
    if all(base in "ATGC" for base in column):
        shared_positions.append(position)

print("Total alignment positions:", alignment.get_alignment_length())
print("Number of species:", len(alignment))
print("Shared nucleotide positions:", len(shared_positions))

if shared_positions:
    print("First shared position:", shared_positions[0])
    print("Last shared position:", shared_positions[-1])

    blocks = []
    start = previous = shared_positions[0]

    for position in shared_positions[1:]:
        if position != previous + 1:
            blocks.append((start, previous))
            start = position
        previous = position

    blocks.append((start, previous))

    print("\nContinuous shared regions:")

    for start, end in blocks:
        print(f"{start}-{end} ({end - start + 1} bp)")
else:
    print("No shared nucleotide positions found.")
