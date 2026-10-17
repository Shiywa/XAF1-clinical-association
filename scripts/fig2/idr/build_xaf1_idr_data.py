#!/usr/bin/env python3
"""Build an alignment-aware, maximal mammalian XAF1 disorder dataset.

The input is the current NCBI Datasets mammalian ortholog set for human XAF1
(GeneID 54739).  One representative protein is selected per non-human species
against canonical human XAF1 (NP_059993.2/Q6GPH4), then IUPred2 long-mode
scores are projected through a multiple sequence alignment onto human residue
coordinates.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import importlib.util
import json
import os
from pathlib import Path
import re
import subprocess
import sys
from zipfile import ZipFile

import numpy as np


HUMAN_ACCESSION = "NP_059993.2"
HUMAN_TAXID = "9606"
HUMAN_IDR_START = 189
HUMAN_IDR_END = 257
STANDARD_AA = set("ACDEFGHIKLMNPQRSTVWY")


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def parse_fasta_text(text: str) -> dict[str, dict[str, str]]:
    records: dict[str, dict[str, str]] = {}
    header: str | None = None
    sequence: list[str] = []
    for raw in text.splitlines():
        if raw.startswith(">"):
            if header is not None:
                accession = header.split()[0]
                records[accession] = {"header": header, "sequence": "".join(sequence)}
            header, sequence = raw[1:].strip(), []
        else:
            sequence.append(raw.strip().upper())
    if header is not None:
        accession = header.split()[0]
        records[accession] = {"header": header, "sequence": "".join(sequence)}
    return records


def parse_fasta_file(path: Path) -> dict[str, str]:
    parsed = parse_fasta_text(path.read_text())
    return {key: value["sequence"] for key, value in parsed.items()}


def wrap_fasta(sequence: str, width: int = 70) -> str:
    return "\n".join(sequence[i : i + width] for i in range(0, len(sequence), width))


def write_fasta(records: list[tuple[str, str]], path: Path) -> None:
    with path.open("w") as handle:
        for name, sequence in records:
            handle.write(f">{name}\n{wrap_fasta(sequence)}\n")


def load_ncbi_packages(package_dir: Path):
    proteins: dict[str, dict[str, str]] = {}
    genes: dict[str, dict] = {}
    products: dict[str, dict] = {}
    package_hashes: dict[str, str] = {}

    for package in sorted(package_dir.glob("batch_*.zip")):
        package_hashes[package.name] = sha256(package)
        with ZipFile(package) as archive:
            assert archive.testzip() is None
            base = "ncbi_dataset/data/"
            for line in archive.read(base + "data_report.jsonl").decode().splitlines():
                if line.strip():
                    row = json.loads(line)
                    genes[str(row["geneId"])] = row
            for line in archive.read(base + "product_report.jsonl").decode().splitlines():
                if line.strip():
                    row = json.loads(line)
                    products[str(row["geneId"])] = row
            proteins.update(parse_fasta_text(archive.read(base + "protein.faa").decode()))
    return proteins, genes, products, package_hashes


def product_metadata(products: dict[str, dict]) -> dict[str, dict]:
    metadata: dict[str, dict] = {}
    for gene_id, gene in products.items():
        for transcript in gene.get("transcripts", []):
            protein = transcript.get("protein")
            if not protein:
                continue
            accession = protein.get("accessionVersion")
            if not accession:
                continue
            row = {
                "accession": accession,
                "gene_id": gene_id,
                "tax_id": str(gene["taxId"]),
                "taxname": gene["taxname"],
                "common_name": gene.get("commonName", ""),
                "symbol": gene.get("symbol", ""),
                "protein_name": protein.get("name", ""),
                "protein_length_reported": protein.get("length", ""),
                "isoform_name": protein.get("isoformName", ""),
                "transcript_accession": transcript.get("accessionVersion", ""),
                "transcript_name": transcript.get("name", ""),
                "transcript_type": transcript.get("type", ""),
                "ensembl_protein": protein.get("ensemblProtein", ""),
                "ensembl_transcript": transcript.get("ensemblTranscript", ""),
            }
            # The same protein can appear in more than one assembly/transcript entry.
            # Prefer a curated transcript record if one is available.
            old = metadata.get(accession)
            if old is None or (
                old["transcript_type"] != "PROTEIN_CODING"
                and row["transcript_type"] == "PROTEIN_CODING"
            ):
                metadata[accession] = row
    return metadata


def run_mafft(mafft: Path, input_fasta: Path, output_fasta: Path, high_accuracy: bool) -> list[str]:
    command = [str(mafft)]
    if high_accuracy:
        command += ["--localpair", "--maxiterate", "1000"]
    else:
        command += ["--auto"]
    command += ["--thread", "-1", "--quiet", str(input_fasta)]
    environment = os.environ.copy()
    # The conda package wrapper retains its build-time prefix. Point it to the
    # binaries installed beside the selected executable.
    conda_libexec = mafft.parent.parent / "libexec" / "mafft"
    if conda_libexec.exists():
        environment["MAFFT_BINARIES"] = str(conda_libexec)
    else:
        environment.pop("MAFFT_BINARIES", None)
    with output_fasta.open("w") as output:
        subprocess.run(command, check=True, stdout=output, env=environment)
    return command


def alignment_metrics(human_aligned: str, candidate_aligned: str) -> dict[str, float]:
    human_position = 0
    overlap = 0
    matches = 0
    idr_overlap = 0
    for human_aa, candidate_aa in zip(human_aligned, candidate_aligned):
        if human_aa != "-":
            human_position += 1
            if candidate_aa != "-":
                overlap += 1
                matches += human_aa == candidate_aa
                if HUMAN_IDR_START <= human_position <= HUMAN_IDR_END:
                    idr_overlap += 1
    human_length = sum(aa != "-" for aa in human_aligned)
    idr_length = HUMAN_IDR_END - HUMAN_IDR_START + 1
    return {
        "human_coverage": overlap / human_length,
        "pairwise_identity": matches / overlap if overlap else 0.0,
        "human_idr_coverage": idr_overlap / idr_length,
        "overlap_residues": overlap,
        "idr_overlap_residues": idr_overlap,
    }


def load_iupred(iupred_source: Path):
    module_path = iupred_source / "iupred" / "iupred2a.py"
    spec = importlib.util.spec_from_file_location("iupred2a_local", module_path)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"Cannot load IUPred2 implementation from {module_path}")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module.iupred, module_path


def project_scores(human_aligned: str, candidate_aligned: str, scores: list[float]):
    human_position = 0
    candidate_position = 0
    projected: dict[int, float] = {}
    for human_aa, candidate_aa in zip(human_aligned, candidate_aligned):
        if human_aa != "-":
            human_position += 1
        if candidate_aa != "-":
            candidate_position += 1
        if human_aa != "-" and candidate_aa != "-":
            projected[human_position] = scores[candidate_position - 1]
    return projected


def write_tsv(path: Path, rows: list[dict], fieldnames: list[str]) -> None:
    with path.open("w", newline="", encoding="utf-8-sig") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--package-dir", type=Path, required=True)
    parser.add_argument("--gene-report", type=Path, required=True)
    parser.add_argument("--iupred-source", type=Path, required=True)
    parser.add_argument("--mafft", type=Path, required=True)
    parser.add_argument("--outdir", type=Path, required=True)
    parser.add_argument("--min-human-coverage", type=float, default=0.70)
    parser.add_argument("--min-idr-coverage", type=float, default=0.80)
    parser.add_argument("--min-identity", type=float, default=0.20)
    args = parser.parse_args()

    args.outdir.mkdir(parents=True, exist_ok=True)
    proteins, genes, products, package_hashes = load_ncbi_packages(args.package_dir)
    metadata = product_metadata(products)
    human_sequence = proteins[HUMAN_ACCESSION]["sequence"]
    if len(human_sequence) != 301:
        raise RuntimeError(f"Unexpected canonical human XAF1 length: {len(human_sequence)}")

    candidates: list[dict] = []
    for accession, row in metadata.items():
        if accession not in proteins or row["tax_id"] == HUMAN_TAXID:
            continue
        sequence = proteins[accession]["sequence"]
        nonstandard = sum(aa not in STANDARD_AA for aa in sequence)
        candidates.append(
            {
                **row,
                "sequence": sequence,
                "sequence_length": len(sequence),
                "nonstandard_residues": nonstandard,
                "nonstandard_fraction": nonstandard / len(sequence),
                "curated_refseq": accession.startswith("NP_"),
            }
        )

    initial_fasta = args.outdir / "all_candidate_proteins.fasta"
    write_fasta(
        [(HUMAN_ACCESSION, human_sequence)]
        + [(row["accession"], row["sequence"]) for row in candidates],
        initial_fasta,
    )
    initial_alignment = args.outdir / "all_candidate_proteins_mafft.fasta"
    initial_mafft_command = run_mafft(args.mafft, initial_fasta, initial_alignment, False)
    aligned = parse_fasta_file(initial_alignment)
    human_aligned = aligned[HUMAN_ACCESSION]

    by_taxon: dict[str, list[dict]] = {}
    for row in candidates:
        metrics = alignment_metrics(human_aligned, aligned[row["accession"]])
        row.update(metrics)
        row["selection_score"] = (
            0.45 * row["human_coverage"]
            + 0.30 * row["pairwise_identity"]
            + 0.25 * row["human_idr_coverage"]
            + (0.015 if row["curated_refseq"] else 0.0)
            - 0.02 * row["nonstandard_fraction"]
        )
        by_taxon.setdefault(row["tax_id"], []).append(row)

    candidate_rows: list[dict] = []
    representatives: list[dict] = []
    for tax_id, rows in sorted(by_taxon.items(), key=lambda item: item[1][0]["taxname"]):
        ranked = sorted(
            rows,
            key=lambda row: (
                row["selection_score"],
                row["curated_refseq"],
                row["human_idr_coverage"],
                row["human_coverage"],
                -abs(row["sequence_length"] - len(human_sequence)),
                row["accession"],
            ),
            reverse=True,
        )
        for rank, row in enumerate(ranked, 1):
            row["rank_within_species"] = rank
            candidate_rows.append(row)
        best = dict(ranked[0])
        reasons: list[str] = []
        if best["human_coverage"] < args.min_human_coverage:
            reasons.append("human_coverage")
        if best["human_idr_coverage"] < args.min_idr_coverage:
            reasons.append("human_IDR_coverage")
        if best["pairwise_identity"] < args.min_identity:
            reasons.append("pairwise_identity")
        if best["nonstandard_fraction"] > 0.02:
            reasons.append("nonstandard_residues")
        best["included"] = not reasons
        best["exclusion_reason"] = ";".join(reasons)
        representatives.append(best)

    included = [row for row in representatives if row["included"]]
    excluded = [row for row in representatives if not row["included"]]
    if len(included) < 150:
        raise RuntimeError(f"Only {len(included)} species passed QC; inspect thresholds and inputs")

    metric_fields = [
        "rank_within_species",
        "included",
        "exclusion_reason",
        "tax_id",
        "taxname",
        "common_name",
        "gene_id",
        "accession",
        "transcript_accession",
        "transcript_type",
        "isoform_name",
        "sequence_length",
        "curated_refseq",
        "nonstandard_residues",
        "human_coverage",
        "human_idr_coverage",
        "pairwise_identity",
        "overlap_residues",
        "idr_overlap_residues",
        "selection_score",
    ]
    write_tsv(args.outdir / "all_candidate_metrics.tsv", candidate_rows, metric_fields)
    write_tsv(args.outdir / "selected_species.tsv", included, metric_fields)
    write_tsv(args.outdir / "excluded_species.tsv", excluded, metric_fields)

    selected_fasta = args.outdir / "selected_representatives.fasta"
    write_fasta(
        [(HUMAN_ACCESSION, human_sequence)]
        + [(row["accession"], row["sequence"]) for row in included],
        selected_fasta,
    )
    final_alignment = args.outdir / "selected_representatives_mafft_linsi.fasta"
    final_mafft_command = run_mafft(args.mafft, selected_fasta, final_alignment, True)
    final_aligned = parse_fasta_file(final_alignment)
    final_human_aligned = final_aligned[HUMAN_ACCESSION]

    iupred, iupred_module_path = load_iupred(args.iupred_source)
    score_rows: list[dict] = []
    projected_by_position: dict[int, list[float]] = {
        position: [] for position in range(1, len(human_sequence) + 1)
    }
    for row in included:
        scores, _ = iupred(row["sequence"], mode="long")
        projected = project_scores(final_human_aligned, final_aligned[row["accession"]], scores)
        for position, score in projected.items():
            projected_by_position[position].append(score)
            score_rows.append(
                {
                    "tax_id": row["tax_id"],
                    "taxname": row["taxname"],
                    "accession": row["accession"],
                    "human_position": position,
                    "iupred2_long_score": f"{score:.8f}",
                }
            )
    human_scores, _ = iupred(human_sequence, mode="long")

    summary_rows: list[dict] = []
    for position in range(1, len(human_sequence) + 1):
        values = np.asarray(projected_by_position[position], dtype=float)
        summary_rows.append(
            {
                "human_position": position,
                "human_residue": human_sequence[position - 1],
                "human_iupred2_long": f"{human_scores[position - 1]:.8f}",
                "species_covered": len(values),
                "species_total": len(included),
                "coverage_fraction": f"{len(values) / len(included):.8f}",
                "cross_species_mean": f"{np.mean(values):.8f}",
                "cross_species_median": f"{np.median(values):.8f}",
                "cross_species_q10": f"{np.quantile(values, 0.10):.8f}",
                "cross_species_q90": f"{np.quantile(values, 0.90):.8f}",
                "fraction_score_gt_0.5": f"{np.mean(values > 0.5):.8f}",
            }
        )
    write_tsv(
        args.outdir / "position_summary_human_coordinates.tsv",
        summary_rows,
        list(summary_rows[0]),
    )
    write_tsv(
        args.outdir / "projected_species_scores.tsv",
        score_rows,
        ["tax_id", "taxname", "accession", "human_position", "iupred2_long_score"],
    )

    gene_report_rows = [json.loads(line) for line in args.gene_report.read_text().splitlines() if line]
    manifest = {
        "analysis_date": "2026-10-17",
        "input": {
            "source": "NCBI Datasets mammalian ortholog set for human XAF1 GeneID 54739",
            "ncbi_datasets_version": "18.39.0",
            "gene_records": len(gene_report_rows),
            "unique_taxa": len({str(row["tax_id"]) for row in gene_report_rows}),
            "nonhuman_candidate_taxa": len(by_taxon),
            "gene_report_sha256": sha256(args.gene_report),
            "package_sha256": package_hashes,
        },
        "representative_selection": {
            "human_reference": HUMAN_ACCESSION,
            "human_reference_length": len(human_sequence),
            "human_idr_coordinates": [HUMAN_IDR_START, HUMAN_IDR_END],
            "protein_candidates": len(candidates),
            "nonhuman_species_included": len(included),
            "nonhuman_species_excluded": len(excluded),
            "thresholds": {
                "minimum_human_coverage": args.min_human_coverage,
                "minimum_human_IDR_coverage": args.min_idr_coverage,
                "minimum_pairwise_identity": args.min_identity,
                "maximum_nonstandard_fraction": 0.02,
            },
            "selection_score": "0.45*human_coverage + 0.30*identity + 0.25*IDR_coverage + 0.015 if NP_ - 0.02*nonstandard_fraction",
        },
        "alignment": {
            "software": "MAFFT v7.526",
            "initial_command": initial_mafft_command,
            "final_command": final_mafft_command,
            "coordinates": "IUPred scores projected through final alignment onto human NP_059993.2 residue positions",
        },
        "disorder_prediction": {
            "algorithm": "IUPred2 long mode",
            "python_package": "saezlab/iupred 0.1.2",
            "source_url": "https://github.com/saezlab/iupred",
            "source_archive_sha256": "acd6fd223b368e178bbbef2c0d231174cebf0f630317133d8c7160cd811ff3f8",
            "module_path": str(iupred_module_path),
            "module_sha256": sha256(iupred_module_path),
        },
        "summary": {
            "line": "median across non-human mammals",
            "ribbon": "10th-90th percentile across non-human mammals",
            "human_in_cross_species_summary": False,
        },
    }
    (args.outdir / "analysis_manifest.json").write_text(
        json.dumps(manifest, ensure_ascii=False, indent=2)
    )
    print(json.dumps(manifest, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
