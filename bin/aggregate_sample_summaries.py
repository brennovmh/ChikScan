#!/usr/bin/env python3

import argparse
import csv
import sys


def read_rows(paths):
    fieldnames = None
    rows = []

    for path in sorted(paths):
        with open(path, newline="") as handle:
            reader = csv.DictReader(handle)
            if reader.fieldnames is None:
                raise ValueError(f"{path}: missing CSV header")
            if fieldnames is None:
                fieldnames = reader.fieldnames
            elif reader.fieldnames != fieldnames:
                raise ValueError(f"{path}: CSV header does not match previous summaries")
            rows.extend(reader)

    if fieldnames is None:
        raise ValueError("No summary CSV files were provided")

    return fieldnames, rows


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--summary", nargs="+", required=True)
    parser.add_argument("--output", required=True)
    args = parser.parse_args()

    fieldnames, rows = read_rows(args.summary)
    with open(args.output, "w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)


if __name__ == "__main__":
    try:
        main()
    except Exception as error:
        print(f"ERROR: {error}", file=sys.stderr)
        sys.exit(1)
