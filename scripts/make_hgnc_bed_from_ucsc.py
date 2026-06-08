#!/usr/bin/env python3
"""Build a BED-like HGNC annotation table from the UCSC HGNC BigBed track.

Inputs:
  1) UCSC HGNC BigBed converted to BED using bigBedToBed
  2) Official HGNC complete set TSV

Output columns:
  chrom, start, end, hgnc_symbol, hgnc_id, hgnc_name,
  locus_group, locus_type, status
"""

import argparse
import csv
import gzip
import re
import sys
from typing import Dict, Iterable, TextIO

HGNC_RE = re.compile(r"HGNC:\d+")


def open_text(path: str, mode: str = "rt") -> TextIO:
    if path.endswith(".gz"):
        return gzip.open(path, mode)
    return open(path, mode)


def load_hgnc_complete(path: str) -> Dict[str, Dict[str, str]]:
    records: Dict[str, Dict[str, str]] = {}

    with open_text(path, "rt") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        required = ["hgnc_id", "symbol", "name", "locus_group", "locus_type", "status"]
        missing = [x for x in required if x not in (reader.fieldnames or [])]
        if missing:
            raise ValueError(
                "HGNC complete set is missing expected columns: {}".format(
                    ", ".join(missing)
                )
            )

        for row in reader:
            hgnc_id = (row.get("hgnc_id") or "").strip()
            if not hgnc_id:
                continue
            records[hgnc_id] = {
                "symbol": (row.get("symbol") or "").strip(),
                "name": (row.get("name") or "").strip(),
                "locus_group": (row.get("locus_group") or "").strip(),
                "locus_type": (row.get("locus_type") or "").strip(),
                "status": (row.get("status") or "").strip(),
            }

    return records


def find_hgnc_id(fields: Iterable[str]) -> str:
    for value in fields:
        match = HGNC_RE.search(value)
        if match:
            return match.group(0)
    return ""


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Convert UCSC HGNC BigBed-derived BED into a BED-like table with official HGNC metadata."
    )
    parser.add_argument("--ucsc-hgnc-bed", required=True)
    parser.add_argument("--hgnc-complete-set", required=True)
    parser.add_argument("--out", required=True)
    args = parser.parse_args()

    hgnc_meta = load_hgnc_complete(args.hgnc_complete_set)

    n_total = 0
    n_with_id = 0
    n_joined = 0

    with open_text(args.ucsc_hgnc_bed, "rt") as in_handle, open_text(args.out, "wt") as out:
        for line in in_handle:
            if not line.strip() or line.startswith("#") or line.startswith("track"):
                continue

            fields = line.rstrip("\n").split("\t")
            if len(fields) < 3:
                continue

            n_total += 1
            chrom, start, end = fields[0], fields[1], fields[2]
            hgnc_id = find_hgnc_id(fields[3:])

            if not hgnc_id:
                continue

            n_with_id += 1
            meta = hgnc_meta.get(hgnc_id, {})
            if meta:
                n_joined += 1

            symbol = meta.get("symbol", "") or hgnc_id
            name = meta.get("name", "")
            locus_group = meta.get("locus_group", "")
            locus_type = meta.get("locus_type", "")
            status = meta.get("status", "")

            out.write(
                "\t".join(
                    [
                        chrom,
                        start,
                        end,
                        symbol,
                        hgnc_id,
                        name,
                        locus_group,
                        locus_type,
                        status,
                    ]
                )
                + "\n"
            )

    sys.stderr.write(
        "HGNC BED conversion complete: raw_records={}; with_hgnc_id={}; joined_to_complete_set={}\n".format(
            n_total, n_with_id, n_joined
        )
    )

    if n_total == 0:
        raise RuntimeError("No rows were read from the UCSC HGNC BED file.")
    if n_with_id == 0:
        raise RuntimeError(
            "No HGNC IDs were found in the UCSC HGNC BED file. "
            "Check whether bigBedToBed produced the expected track fields."
        )
    if n_joined == 0:
        raise RuntimeError(
            "No UCSC HGNC intervals joined to the official HGNC complete set. "
            "Check the HGNC ID format and hgnc_complete_set.txt columns."
        )


if __name__ == "__main__":
    main()
