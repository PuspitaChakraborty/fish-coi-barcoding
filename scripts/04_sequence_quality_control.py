from pathlib import Path
import csv

from Bio import SeqIO

PROJECT = Path(__file__).resolve().parent.parent
RAW_DIR = PROJECT / "data" / "raw"
OUTPUT = PROJECT / "results" / "sequence_quality_control.csv"

VALID_BASES = set("ACGT")


def get_cds(record):
    """Find the annotated COI/COX1 CDS feature."""
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


def analyze_record(path):
    record = SeqIO.read(path, "genbank")
    feature = get_cds(record)

    if feature is None:
        return {
            "accession": record.id,
            "status": "REVIEW",
            "message": "COI/COX1 CDS feature not found",
        }

    sequence = feature.extract(record.seq).upper()

    codon_start = int(feature.qualifiers.get("codon_start", ["1"])[0])
    transl_table = int(feature.qualifiers.get("transl_table", ["1"])[0])

    coding_sequence = sequence[codon_start - 1:]

    usable_length = len(coding_sequence) - (len(coding_sequence) % 3)
    translatable_sequence = coding_sequence[:usable_length]

    ambiguous_count = sum(
        base not in VALID_BASES for base in str(sequence)
    )

    generated_translation = str(
        translatable_sequence.translate(table=transl_table)
    )

    annotated_translation = "".join(
        feature.qualifiers.get("translation", [])
    ).replace(" ", "").replace("\n", "")

    generated_normalized = generated_translation.rstrip("*")
    annotated_normalized = annotated_translation.rstrip("*")

    internal_stops = generated_normalized.count("*")
    translation_match = generated_normalized == annotated_normalized

    warnings = []

    if ambiguous_count:
        warnings.append("ambiguous_nucleotides")

    if internal_stops:
        warnings.append("internal_stop_codons")

    if not translation_match:
        warnings.append("translation_mismatch")

    if len(coding_sequence) % 3:
        warnings.append("incomplete_terminal_codon")

    critical_warnings = [
        warning for warning in warnings
        if warning != "incomplete_terminal_codon"
    ]

    if critical_warnings:
        status = "REVIEW"
    elif "incomplete_terminal_codon" in warnings:
        status = "PASS_WITH_NOTE"
    else:
        status = "PASS"

    return {
        "accession": record.id,
        "organism": record.annotations.get("organism", ""),
        "sequence_length_bp": len(sequence),
        "codon_start": codon_start,
        "translation_table": transl_table,
        "ambiguous_nucleotides": ambiguous_count,
        "translatable_length_bp": usable_length,
        "generated_protein_length": len(generated_translation),
        "annotated_protein_length": len(annotated_translation),
        "internal_stop_codons": internal_stops,
        "translation_match": translation_match,
        "status": status,
        "warnings": "; ".join(warnings),
    }


def main():
    OUTPUT.parent.mkdir(parents=True, exist_ok=True)
    results = []

    for path in sorted(RAW_DIR.glob("*.gb")):
        try:
            result = analyze_record(path)
        except Exception as error:
            result = {
                "accession": path.stem,
                "status": "ERROR",
                "message": str(error),
            }

        results.append(result)

        print(
            f"{result['status']:14} "
            f"{result['accession']} "
            f"{result.get('organism', '')} "
            f"{result.get('warnings', result.get('message', ''))}"
        )

    if not results:
        raise SystemExit("ERROR: No GenBank files found.")

    fieldnames = list(dict.fromkeys(
        key for result in results for key in result
    ))

    with OUTPUT.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(results)

    passed = sum(r["status"] == "PASS" for r in results)
    passed_with_note = sum(
        r["status"] == "PASS_WITH_NOTE" for r in results
    )
    review = sum(r["status"] == "REVIEW" for r in results)
    errors = sum(r["status"] == "ERROR" for r in results)

    print("\nSequence QC summary")
    print(f"Total records: {len(results)}")
    print(f"PASS: {passed}")
    print(f"PASS_WITH_NOTE: {passed_with_note}")
    print(f"REVIEW: {review}")
    print(f"ERROR: {errors}")
    print(f"Report saved: {OUTPUT}")


if __name__ == "__main__":
    main()
