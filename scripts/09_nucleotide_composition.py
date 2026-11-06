from Bio import AlignIO
import csv

INPUT_FILE = "data/processed/coi_aligned_trimmed.fasta"
OUTPUT_FILE = "results/nucleotide_composition.csv"

alignment = AlignIO.read(INPUT_FILE, "fasta")

results = []

for record in alignment:
    sequence = str(record.seq).upper()
    length = len(sequence)

    counts = {
        base: sequence.count(base)
        for base in "ATGC"
    }

    gc_content = (
        (counts["G"] + counts["C"]) / length
    ) * 100

    at_content = (
        (counts["A"] + counts["T"]) / length
    ) * 100

    results.append({
        "species": record.id,
        "sequence_length": length,
        "A_count": counts["A"],
        "T_count": counts["T"],
        "G_count": counts["G"],
        "C_count": counts["C"],
        "A_percent": round(counts["A"] / length * 100, 2),
        "T_percent": round(counts["T"] / length * 100, 2),
        "G_percent": round(counts["G"] / length * 100, 2),
        "C_percent": round(counts["C"] / length * 100, 2),
        "GC_percent": round(gc_content, 2),
        "AT_percent": round(at_content, 2)
    })

with open(OUTPUT_FILE, "w", newline="", encoding="utf-8") as f:
    writer = csv.DictWriter(f, fieldnames=results[0].keys())
    writer.writeheader()
    writer.writerows(results)

print("Nucleotide composition analysis completed.")
print("Species analysed:", len(results))
print("Alignment length:", alignment.get_alignment_length())
print("Results saved:", OUTPUT_FILE)
