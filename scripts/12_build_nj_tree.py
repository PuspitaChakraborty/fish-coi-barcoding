import csv
from pathlib import Path
from Bio.Phylo.TreeConstruction import DistanceMatrix, DistanceTreeConstructor
from Bio import Phylo

INPUT = Path("results/k2p_distance_matrix.csv")
OUTPUT = Path("results/nj_tree.nwk")

with open(INPUT, newline="") as file:
    rows = list(csv.reader(file))

species = rows[0][1:]
matrix = []

for row in rows[1:]:
    matrix.append([float(value) for value in row[1:]])

if len(species) != 15 or len(matrix) != 15:
    raise ValueError("Expected a 15-species distance matrix.")

# Validate matrix dimensions and symmetry
for i in range(15):
    if len(matrix[i]) != 15:
        raise ValueError("Invalid matrix dimensions.")
    for j in range(15):
        if abs(matrix[i][j] - matrix[j][i]) > 1e-9:
            raise ValueError("Distance matrix is not symmetrical.")

# Biopython requires only the lower triangle, including diagonal
lower_triangle = [
    matrix[i][:i + 1]
    for i in range(15)
]

distance_matrix = DistanceMatrix(
    names=species,
    matrix=lower_triangle
)

constructor = DistanceTreeConstructor()
tree = constructor.nj(distance_matrix)

# Save tree in Newick format
OUTPUT.parent.mkdir(parents=True, exist_ok=True)
Phylo.write(tree, OUTPUT, "newick")

print("Neighbor-Joining tree constructed successfully.")
print(f"Species included: {len(species)}")
print(f"Tree saved: {OUTPUT}")

print("\nTree structure:")
Phylo.draw_ascii(tree)
