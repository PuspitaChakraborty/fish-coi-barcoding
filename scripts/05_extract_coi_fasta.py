from pathlib import Path
import csv

from Bio import SeqIO
from Bio.SeqRecord import SeqRecord


PROJECT = Path(__file__).resolve().parent.parent
RAW_DIR = PROJECT / "data" / "raw"
ACCESSIONS_FILE = PROJECT / "data" / "accessions.csv"
QC_FILE = PROJECT / "results" / "sequence_quality_control.csv"

OUTPUT_DIR = PROJECT / "data" / "processed"
FASTA_OUTPUT = OUTPUT_DIR / "coi_sequences.fasta"
SUMMARY_OUTPUT = PROJECT / "results" / "coi_extraction_summary.csv"

VALID_BASES = set("ACGT")


def find_coi_cds(record):
    for feature in record.features:
        if feature.type != "CDS":
            continue

        gene = " ".join(feature.qualifiers.get("gene", [])).upper()
        product = " ".join(feature.qualifiers.get("product", [])).lower()

        if (
            gene in {"COI", "CO1", "COX1"}
            or "cytochrome oxidase subunit 1" in product
            or "cytochrome c oxidase subunit i" in product
        ):
            return feature

    return None


def main():
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

    with QC_FILE.open(newline="", encoding="utf-8") as f:
        qc_status = {
            row["accession"]: row["status"]
            for row in csv.DictReader(f)
        }

    with ACCESSIONS_FILE.open(newline="", encoding="utf-8") as f:
        accessions = list(csv.DictReader(f))

    fasta_records = []
    summaries = []

    for item in accessions:
        accession = item["accession"]
        species = item["target_species"]

        if qc_status.get(accession) not in {"PASS", "PASS_WITH_NOTE"}:
            raise ValueError(f"{accession}: QC status is not acceptable")

        gb_file = RAW_DIR / f"{accession}.gb"

        record = SeqIO.read(gb_file, "genbank")

        if record.id != accession:
            raise ValueError(
                f"Accession mismatch: expected {accession}, found {record.id}"
            )

        feature = find_coi_cds(record)

        if feature is None:
            raise ValueError(f"{accession}: COI CDS not found")

        cds = feature.extract(record.seq).upper()

        codon_start = int(
            feature.qualifiers.get("codon_start", ["1"])[0]
        )

        coding_sequence = cds[codon_start - 1:]

        trailing_bases_removed = len(coding_sequence) % 3

        if trailing_bases_removed:
            coding_sequence = coding_sequence[:-trailing_bases_removed]

        if not coding_sequence:
            raise ValueError(f"{accession}: Empty extracted sequence")

        if set(str(coding_sequence)) - VALID_BASES:
            raise ValueError(f"{accession}: Ambiguous nucleotides detected")

        fasta_id = f"{species.replace(' ', '_')}_{accession}"

        fasta_records.append(
            SeqRecord(
                coding_sequence,
                id=fasta_id,
                description=f"target_species={species} accession={accession}",
            )
        )

        summaries.append({
            "target_species": species,
            "accession": accession,
            "original_cds_length": len(cds),
            "codon_start": codon_start,
            "leading_bases_removed": codon_start - 1,
            "trailing_bases_removed": trailing_bases_removed,
            "final_sequence_length": len(coding_sequence),
        })

        print(
            f"EXTRACTED {accession} | {species} | "
            f"{len(coding_sequence)} bp"
        )

    SeqIO.write(fasta_records, FASTA_OUTPUT, "fasta")

    with SUMMARY_OUTPUT.open("w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=summaries[0].keys())
        writer.writeheader()
        writer.writerows(summaries)

    print("\nExtraction completed.")
    print(f"Sequences extracted: {len(fasta_records)}")
    print(f"FASTA saved: {FASTA_OUTPUT}")
    print(f"Summary saved: {SUMMARY_OUTPUT}")


if __name__ == "__main__":
    main()
