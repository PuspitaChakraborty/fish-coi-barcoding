from pathlib import Path
import matplotlib.pyplot as plt
from Bio import Phylo

# Input and output
INPUT = "results/nj_bootstrap_tree.nwk"
OUTPUT = Path("figures/phylogenetic_tree.png")

# Read the bootstrap NJ tree
tree = Phylo.read(INPUT, "newick")

# Format species labels
for clade in tree.find_clades():
    if clade.name:
        clade.name = clade.name.replace("_", " ")

# Create high-resolution figure
fig, ax = plt.subplots(figsize=(13, 8), dpi=300)

# Draw phylogenetic tree with bootstrap support
Phylo.draw(
    tree,
    axes=ax,
    do_show=False,
    show_confidence=True,
    label_func=lambda node: node.name if node.is_terminal() else ""
)

# Remove unnecessary vertical axis information
ax.set_ylabel("")
ax.set_yticks([])
ax.tick_params(axis="y", left=False)

# Improve title and horizontal axis
ax.set_title(
    "COI-Based Neighbor-Joining Phylogenetic Tree of 15 Freshwater Fish Taxa",
    fontsize=13,
    fontweight="bold",
    pad=18
)

ax.set_xlabel(
    "Kimura 2-Parameter (K2P) Genetic Distance",
    fontsize=11
)

# Improve layout
ax.tick_params(axis="x", labelsize=9)
plt.tight_layout()

# Ensure output folder exists
OUTPUT.parent.mkdir(parents=True, exist_ok=True)

# Save high-resolution image
plt.savefig(
    OUTPUT,
    dpi=300,
    bbox_inches="tight",
    facecolor="white"
)

plt.close()

print("Phylogenetic tree rendered successfully.")
print(f"Output: {OUTPUT}")
print("Resolution: 300 DPI")
