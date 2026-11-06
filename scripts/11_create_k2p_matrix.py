import csv
from pathlib import Path

INPUT = Path("results/pairwise_genetic_distances.csv")
OUTPUT = Path("results/k2p_distance_matrix.csv")

with open(INPUT, newline="") as file:
    rows = list(csv.DictReader(file))

species = sorted({
    row["species_1"] for row in rows
} | {
    row["species_2"] for row in rows
})

matrix = {
    name: {other: 0.0 for other in species}
    for name in species
}

for row in rows:
    sp1 = row["species_1"]
    sp2 = row["species_2"]
    distance = row["K2P_distance"]

    if distance == "":
        raise ValueError(f"Undefined K2P distance: {sp1} vs {sp2}")

    distance = float(distance)

    matrix[sp1][sp2] = distance
    matrix[sp2][sp1] = distance

if len(species) != 15 or len(rows) != 105:
    raise ValueError("Unexpected number of species or comparisons")

OUTPUT.parent.mkdir(parents=True, exist_ok=True)

with open(OUTPUT, "w", newline="") as file:
    writer = csv.writer(file)
    writer.writerow(["species"] + species)

    for name in species:
        writer.writerow(
            [name] + [matrix[name][other] for other in species]
        )

print("K2P distance matrix created successfully.")
print(f"Species: {len(species)}")
print(f"Matrix dimensions: {len(species)} x {len(species)}")
print(f"Output saved: {OUTPUT}")
