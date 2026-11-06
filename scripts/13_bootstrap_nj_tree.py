import random
import math
from pathlib import Path
from collections import Counter

from Bio import AlignIO, Phylo
from Bio.Phylo.TreeConstruction import (
    DistanceMatrix,
    DistanceTreeConstructor
)

INPUT = Path("data/processed/coi_aligned_trimmed.fasta")
OUTPUT = Path("results/nj_bootstrap_tree.nwk")

REPLICATES = 1000
random.seed(42)

alignment = AlignIO.read(INPUT, "fasta")
names = [record.id for record in alignment]
sequences = [str(record.seq).upper() for record in alignment]

if len(names) != 15:
    raise ValueError("Expected 15 species.")

length = len(sequences[0])

if length != 558 or any(len(s) != length for s in sequences):
    raise ValueError("Unexpected alignment length.")

if any(set(s) - set("ATGC") for s in sequences):
    raise ValueError("Invalid nucleotide found.")

constructor = DistanceTreeConstructor()


def calculate_tree(sampled_sequences):
    lower = []

    for i in range(len(names)):
        row = []

        for j in range(i + 1):
            if i == j:
                row.append(0.0)
                continue

            transitions = 0
            transversions = 0

            for a, b in zip(sampled_sequences[i], sampled_sequences[j]):
                if a == b:
                    continue

                if (a, b) in [
                    ("A", "G"), ("G", "A"),
                    ("C", "T"), ("T", "C")
                ]:
                    transitions += 1
                else:
                    transversions += 1

            P = transitions / length
            Q = transversions / length

            term1 = 1 - 2 * P - Q
            term2 = 1 - 2 * Q

            if term1 <= 0 or term2 <= 0:
                raise ValueError("Undefined K2P distance.")

            distance = (
                -0.5 * math.log(term1)
                -0.25 * math.log(term2)
            )

            row.append(distance)

        lower.append(row)

    dm = DistanceMatrix(names, lower)
    return constructor.nj(dm)


def get_splits(tree):
    all_names = set(names)
    splits = set()

    for clade in tree.get_nonterminals():
        descendants = {
            tip.name for tip in clade.get_terminals()
        }

        other_side = all_names - descendants

        if len(descendants) > 1 and len(other_side) > 1:
            side1 = tuple(sorted(descendants))
            side2 = tuple(sorted(other_side))
            splits.add(min(side1, side2))

    return splits


# Original tree
original_tree = calculate_tree(sequences)

# Bootstrap resampling
support = Counter()

for replicate in range(REPLICATES):

    columns = [
        random.randrange(length)
        for _ in range(length)
    ]

    sampled_sequences = [
        "".join(seq[i] for i in columns)
        for seq in sequences
    ]

    bootstrap_tree = calculate_tree(sampled_sequences)
    support.update(get_splits(bootstrap_tree))

    if (replicate + 1) % 100 == 0:
        print(f"Completed {replicate + 1}/{REPLICATES}")

# Assign support values
assigned = 0

for clade in original_tree.get_nonterminals():
    descendants = {
        tip.name for tip in clade.get_terminals()
    }

    other_side = set(names) - descendants

    if len(descendants) > 1 and len(other_side) > 1:

        side1 = tuple(sorted(descendants))
        side2 = tuple(sorted(other_side))
        split = min(side1, side2)

        clade.confidence = round(
            support[split] / REPLICATES * 100, 1
        )

        # Ensure bootstrap values are written as labels in Newick
        clade.name = None
        assigned += 1

if assigned == 0:
    raise ValueError("No bootstrap support values were assigned.")

OUTPUT.parent.mkdir(parents=True, exist_ok=True)
Phylo.write(original_tree, OUTPUT, "newick")

print("\nBootstrap analysis completed.")
print(f"Replicates: {REPLICATES}")
print(f"Internal branches with support: {assigned}")
print(f"Tree saved: {OUTPUT}")
