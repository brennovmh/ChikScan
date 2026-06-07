#!/usr/bin/env python3

import argparse
import csv
import gzip
import sys


def open_text(path):
    if str(path).endswith(".gz"):
        return gzip.open(path, "rt")
    return open(path)


def parse_fasta(path):
    records = []
    current = None
    chunks = []
    with open_text(path) as handle:
        for line in handle:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if current is not None:
                    records.append((current, "".join(chunks).upper()))
                current = line[1:]
                chunks = []
            else:
                chunks.append(line)
    if current is not None:
        records.append((current, "".join(chunks).upper()))
    if not records:
        raise ValueError(f"{path}: no FASTA records found")
    return records


def parse_fastq_sequences(paths, max_reads):
    sequences = []
    for path in paths:
        with open_text(path) as handle:
            for line_number, line in enumerate(handle):
                if line_number % 4 == 1:
                    sequences.append(line.strip().upper())
                    if max_reads and len(sequences) >= max_reads:
                        return sequences
    if not sequences:
        raise ValueError("No FASTQ reads found")
    return sequences


def kmers(sequence, kmer_size):
    for index in range(0, len(sequence) - kmer_size + 1):
        kmer = sequence[index:index + kmer_size]
        if set(kmer) <= {"A", "C", "G", "T"}:
            yield kmer


def reference_kmers(sequence, kmer_size):
    return set(kmers(sequence, kmer_size))


def score_reference(reads, ref_kmers, kmer_size):
    total = 0
    matches = 0
    matched_reads = 0
    for read in reads:
        read_total = 0
        read_matches = 0
        for kmer in kmers(read, kmer_size):
            read_total += 1
            if kmer in ref_kmers:
                read_matches += 1
        total += read_total
        matches += read_matches
        if read_matches:
            matched_reads += 1
    identity = matches / total if total else 0.0
    matched_read_fraction = matched_reads / len(reads) if reads else 0.0
    return matches, total, matched_reads, identity, matched_read_fraction


def classify_confidence(best, runner_up, min_identity, min_matched_read_fraction, min_score_margin):
    if best["score"] == 0:
        return "no_match"
    if best["kmer_identity_raw"] < min_identity:
        return "low"
    if best["matched_read_fraction_raw"] < min_matched_read_fraction:
        return "low"
    if runner_up is not None:
        total = best["total_kmers_raw"]
        margin = (best["score"] - runner_up["score"]) / total if total else 0.0
        if margin < min_score_margin:
            return "low"
    return "high"


def write_fasta(header, sequence, output):
    with open(output, "w") as handle:
        handle.write(f">{header}\n")
        for index in range(0, len(sequence), 80):
            handle.write(sequence[index:index + 80] + "\n")


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--sample-id", required=True)
    parser.add_argument("--reads", nargs="+", required=True)
    parser.add_argument("--references", required=True)
    parser.add_argument("--selected-fasta", required=True)
    parser.add_argument("--output", required=True)
    parser.add_argument("--kmer-size", type=int, default=31)
    parser.add_argument("--max-reads", type=int, default=10000)
    parser.add_argument("--min-kmer-identity", type=float, default=0.05)
    parser.add_argument("--min-matched-read-fraction", type=float, default=0.10)
    parser.add_argument("--min-score-margin", type=float, default=0.01)
    args = parser.parse_args()

    references = parse_fasta(args.references)
    reads = parse_fastq_sequences(args.reads, args.max_reads)
    rows = []

    for reference_index, (header, sequence) in enumerate(references):
        ref_id = header.split()[0].split("|")[0]
        matches, total, matched_reads, identity, matched_read_fraction = score_reference(
            reads,
            reference_kmers(sequence, args.kmer_size),
            args.kmer_size,
        )
        row = {
            "sample_id": args.sample_id,
            "reference_id": ref_id,
            "score": matches,
            "score_margin": "0.000000",
            "total_kmers": total,
            "matched_reads": matched_reads,
            "read_count": len(reads),
            "kmer_identity": f"{identity:.6f}",
            "matched_read_fraction": f"{matched_read_fraction:.6f}",
            "confidence": "not_selected",
            "selected": "false",
            "selection_status": "not_selected",
            "selection_note": "",
            "reference_index": reference_index,
            "header": header,
            "sequence": sequence,
            "kmer_identity_raw": identity,
            "matched_read_fraction_raw": matched_read_fraction,
            "total_kmers_raw": total,
        }
        rows.append(row)

    rows.sort(
        key=lambda row: (
            row["score"],
            row["matched_reads"],
            row["kmer_identity_raw"],
            -row["reference_index"],
        ),
        reverse=True,
    )
    best = rows[0]
    runner_up = rows[1] if len(rows) > 1 else None
    confidence = classify_confidence(
        best,
        runner_up,
        args.min_kmer_identity,
        args.min_matched_read_fraction,
        args.min_score_margin,
    )
    best["confidence"] = confidence
    best["selected"] = "true"
    best["selection_status"] = "selected"
    if runner_up is not None:
        margin = (best["score"] - runner_up["score"]) / best["total_kmers_raw"] if best["total_kmers_raw"] else 0.0
        best["score_margin"] = f"{margin:.6f}"
    if confidence == "low":
        best["selection_note"] = "Selected reference did not meet one or more confidence thresholds"
    elif confidence == "no_match":
        best["selection_note"] = "No exact read k-mer matches met the reference panel; selected fallback reference"

    for row in rows:
        if row is not best and best["total_kmers_raw"]:
            margin = (best["score"] - row["score"]) / best["total_kmers_raw"]
            row["score_margin"] = f"{margin:.6f}"

    write_fasta(best["header"], best["sequence"], args.selected_fasta)
    with open(args.output, "w", newline="") as handle:
        fieldnames = [
            "sample_id",
            "reference_id",
            "score",
            "score_margin",
            "total_kmers",
            "matched_reads",
            "read_count",
            "kmer_identity",
            "matched_read_fraction",
            "confidence",
            "selected",
            "selection_status",
            "selection_note",
        ]
        writer = csv.DictWriter(
            handle,
            fieldnames=fieldnames,
        )
        writer.writeheader()
        for row in rows:
            writer.writerow({key: row[key] for key in fieldnames})


if __name__ == "__main__":
    try:
        main()
    except Exception as error:
        print(f"ERROR: {error}", file=sys.stderr)
        sys.exit(1)
