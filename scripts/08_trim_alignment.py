from Bio import AlignIO
from Bio.Align import MultipleSeqAlignment

INPUT_FILE = "data/processed/coi_aligned.fasta"
OUTPUT_FILE = "data/processed/coi_aligned_trimmed.fasta"

START = 82
END = 639

alignment = AlignIO.read(INPUT_FILE, "fasta")

trimmed_records = []

for record in alignment:
    trimmed_sequence = record.seq[START - 1:END]

    if len(trimmed_sequence) != 558:
        raise ValueError(f"Unexpected length: {record.id}")

    if "-" in trimmed_sequence:
        raise ValueError(f"Gap found in shared region: {record.id}")

    new_record = record[:]
    new_record.seq = trimmed_sequence
    trimmed_records.append(new_record)

trimmed_alignment = MultipleSeqAlignment(trimmed_records)

AlignIO.write(trimmed_alignment, OUTPUT_FILE, "fasta")

print("Trimming completed successfully.")
print("Number of sequences:", len(trimmed_alignment))
print("Trimmed alignment length:", trimmed_alignment.get_alignment_length())
print("Output saved:", OUTPUT_FILE)
