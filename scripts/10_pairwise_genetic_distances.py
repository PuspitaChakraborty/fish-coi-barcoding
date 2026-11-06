from Bio import AlignIO
from itertools import combinations
import csv
import math
from pathlib import Path

INPUT = Path("data/processed/coi_aligned_trimmed.fasta")
OUTPUT = Path("results/pairwise_genetic_distances.csv")

alignment = AlignIO.read(INPUT, "fasta")

if len(alignment) != 15:
    raise ValueError(f"Expected 15 sequences, found {len(alignment)}")

lengths = {len(record.seq) for record in alignment}

if lengths != {558}:
    raise ValueError(f"Expected all sequences to be 558 bp, found {lengths}")

sequences = {}

for record in alignment:
    species = record.id
    sequence = str(record.seq).upper()

    if set(sequence) - set("ATGC"):
        raise ValueError(f"Invalid nucleotide found in {species}")

    sequences[species] = sequence

results = []

for (species1, seq1), (species2, seq2) in combinations(sequences.items(), 2):

    differences = 0
    transitions = 0
    transversions = 0

    for base1, base2 in zip(seq1, seq2):

        if base1 == base2:
            continue

        differences += 1

        if (base1, base2) in [
            ("A", "G"), ("G", "A"),
            ("C", "T"), ("T", "C")
        ]:
            transitions += 1
        else:
            transversions += 1

    p_distance = differences / 558
    P = transitions / 558
    Q = transversions / 558

    term1 = 1 - 2 * P - Q
    term2 = 1 - 2 * Q

    if term1 > 0 and term2 > 0:
        k2p_distance = (
            -0.5 * math.log(term1)
            -0.25 * math.log(term2)
        )
    else:
        k2p_distance = ""

    results.append({
        "species_1": species1,
        "species_2": species2,
        "compared_sites": 558,
        "differences": differences,
        "transitions": transitions,
        "transversions": transversions,
        "p_distance": round(p_distance, 6),
        "K2P_distance": (
            round(k2p_distance, 6)
            if k2p_distance != "" else ""
        )
    })

OUTPUT.parent.mkdir(parents=True, exist_ok=True)

with open(OUTPUT, "w", newline="") as file:
    writer = csv.DictWriter(file, fieldnames=results[0].keys())
    writer.writeheader()
    writer.writerows(results)

if len(results) != 105:
    raise ValueError(f"Expected 105 comparisons, found {len(results)}")

print("Pairwise genetic distance analysis completed.")
print(f"Sequences analysed: {len(sequences)}")
print(f"Pairwise comparisons: {len(results)}")
print(f"Compared sites per pair: 558")
print(f"Results saved: {OUTPUT}")
