#!/usr/bin/env python3

import csv
import re
from pathlib import Path
from Bio import SeqIO

ROOT = Path(__file__).resolve().parent.parent
CSV_FILE = ROOT / "data" / "accessions.csv"
RAW_DIR = ROOT / "data" / "raw"
REPORT = ROOT / "results" / "genbank_validation.csv"

def clean(value):
    return re.sub(r"\s+", " ", str(value or "")).strip()

def validate(row):
    accession = row["accession"]
    expected_organism = clean(row["genbank_organism"])
    expected_length = int(row["expected_coi_length_bp"])
    file_path = RAW_DIR / f"{accession}.gb"

    checks = {
        "accession_match": False,
        "organism": False,
        "mitochondrial": False,
        "COX1_gene": False,
        "protein_product": False,
        "length": False,
    }

    observed_organism = ""
    observed_length = ""
    message = ""

    try:
        if not file_path.exists():
            raise FileNotFoundError("GenBank file missing")

        record = SeqIO.read(file_path, "genbank")

        checks["accession_match"] = record.id == accession
        observed_length = len(record.seq)
        checks["length"] = observed_length == expected_length

        source_features = [
            feature for feature in record.features
            if feature.type == "source"
        ]

        for feature in source_features:
            qualifiers = feature.qualifiers
            observed_organism = clean(
                qualifiers.get("organism", [""])[0]
            )

            if observed_organism == expected_organism:
                checks["organism"] = True

            organelle = clean(
                qualifiers.get("organelle", [""])[0]
            ).lower()

            if "mitochond" in organelle:
                checks["mitochondrial"] = True

        gene_features = [
            feature for feature in record.features
            if feature.type in ("gene", "CDS")
        ]

        for feature in gene_features:
            gene = clean(
                feature.qualifiers.get("gene", [""])[0]
            ).upper()

            if gene in ("COX1", "COI"):
                checks["COX1_gene"] = True

            product = clean(
                feature.qualifiers.get("product", [""])[0]
            ).lower()

            if (
                "cytochrome oxidase subunit 1" in product
                or "cytochrome c oxidase subunit i" in product
            ):
                checks["protein_product"] = True

    except Exception as error:
        message = str(error)

    passed = all(checks.values())
    status = "PASS" if passed else "FAIL"

    failed_checks = [
        name for name, result in checks.items() if not result
    ]

    return {
        "accession": accession,
        "expected_organism": expected_organism,
        "observed_organism": observed_organism,
        "expected_length": expected_length,
        "observed_length": observed_length,
        **checks,
        "status": status,
        "failed_checks": "; ".join(failed_checks),
        "message": message,
    }


def main():
    with open(CSV_FILE, newline="") as handle:
        rows = list(csv.DictReader(handle))

    results = [validate(row) for row in rows]

    REPORT.parent.mkdir(parents=True, exist_ok=True)

    with open(REPORT, "w", newline="") as handle:
        writer = csv.DictWriter(
            handle, fieldnames=results[0].keys()
        )
        writer.writeheader()
        writer.writerows(results)

    for result in results:
        print(
            f"{result['status']:4} "
            f"{result['accession']} "
            f"{result['observed_organism']} "
            f"({result['observed_length']} bp)"
        )

        if result["failed_checks"]:
            print("     Failed:", result["failed_checks"])

        if result["message"]:
            print("     Error:", result["message"])

    passed = sum(r["status"] == "PASS" for r in results)

    print(f"\nTotal records: {len(results)}")
    print(f"Passed: {passed}")
    print(f"Failed: {len(results) - passed}")
    print(f"Report saved: {REPORT}")


if __name__ == "__main__":
    main()
