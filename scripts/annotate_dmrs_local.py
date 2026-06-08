#!/usr/bin/env python3
"""Annotate comb-p DMR regions using local refGene, CpG island, and HGNC BED files.

The output keeps all original comb-p region columns and appends local annotation
columns. The first three columns must be chrom/start/end.
"""

import argparse
import gzip
from collections import defaultdict
from typing import Dict, Iterable, List, TextIO, Tuple


def open_text(path: str, mode: str = "rt") -> TextIO:
    if path.endswith(".gz"):
        return gzip.open(path, mode)
    return open(path, mode)


def split_values(values: Iterable[str]) -> str:
    clean = sorted({v for v in values if v not in ["", ".", "NA", None]})
    return ";".join(clean) if clean else "NA"


def read_bed(path: str, kind: str) -> Dict[str, List[dict]]:
    records: Dict[str, List[dict]] = defaultdict(list)

    with open_text(path, "rt") as handle:
        for line in handle:
            if not line.strip() or line.startswith("#") or line.startswith("track"):
                continue

            fields = line.rstrip("\n").split("\t")
            if len(fields) < 3:
                continue

            try:
                start = int(fields[1])
                end = int(fields[2])
            except ValueError:
                continue

            records[fields[0]].append(
                {
                    "chrom": fields[0],
                    "start": start,
                    "end": end,
                    "fields": fields,
                    "kind": kind,
                }
            )

    for chrom in records:
        records[chrom].sort(key=lambda x: (x["start"], x["end"]))

    return records


def read_regions(path: str) -> Tuple[List[str], List[dict]]:
    header = None
    regions = []

    with open_text(path, "rt") as handle:
        for line in handle:
            if not line.strip():
                continue

            raw = line.rstrip("\n")
            if raw.startswith("#"):
                candidate = raw.lstrip("#").split("\t")
                if len(candidate) >= 3:
                    header = candidate
                continue

            fields = raw.split("\t")
            if len(fields) < 3:
                continue

            try:
                start = int(fields[1])
                end = int(fields[2])
            except ValueError:
                # Skip accidental non-comment header lines.
                continue

            regions.append(
                {
                    "chrom": fields[0],
                    "start": start,
                    "end": end,
                    "fields": fields,
                }
            )

    if not regions:
        raise RuntimeError("No DMR regions were read from {}".format(path))

    if header is None:
        max_cols = max(len(r["fields"]) for r in regions)
        header = ["chrom", "start", "end"] + [
            "combp_col{}".format(i) for i in range(4, max_cols + 1)
        ]

    return header, regions


def interval_distance(region: dict, ann: dict) -> int:
    if ann["end"] > region["start"] and ann["start"] < region["end"]:
        return 0
    if ann["end"] <= region["start"]:
        return region["start"] - ann["end"]
    return ann["start"] - region["end"]


def find_overlaps(region: dict, annotation_by_chrom: Dict[str, List[dict]]) -> List[dict]:
    hits = []
    for ann in annotation_by_chrom.get(region["chrom"], []):
        if ann["start"] >= region["end"]:
            break
        if ann["end"] > region["start"] and ann["start"] < region["end"]:
            hits.append(ann)
    return hits


def find_nearest(region: dict, annotation_by_chrom: Dict[str, List[dict]]) -> Tuple[List[dict], str]:
    anns = annotation_by_chrom.get(region["chrom"], [])
    if not anns:
        return [], "NA"

    best_dist = None
    best = []

    for ann in anns:
        dist = interval_distance(region, ann)
        if best_dist is None or dist < best_dist:
            best_dist = dist
            best = [ann]
        elif dist == best_dist:
            best.append(ann)

    return best, str(best_dist)


def field(rec: dict, idx: int) -> str:
    return rec["fields"][idx] if len(rec["fields"]) > idx else ""


def main() -> None:
    parser = argparse.ArgumentParser(description="Annotate comb-p DMR regions using local annotation files.")
    parser.add_argument("--regions", required=True)
    parser.add_argument("--hgnc-bed", required=True)
    parser.add_argument("--refgene-bed", required=True)
    parser.add_argument("--cpg-bed", required=True)
    parser.add_argument("--manifest", required=True)
    parser.add_argument("--out", required=True)
    args = parser.parse_args()

    region_header, regions = read_regions(args.regions)
    hgnc = read_bed(args.hgnc_bed, "hgnc")
    refgene = read_bed(args.refgene_bed, "refgene")
    cpg = read_bed(args.cpg_bed, "cpg")

    extra_header = [
        "hgnc_overlap_symbols",
        "hgnc_overlap_ids",
        "hgnc_overlap_names",
        "hgnc_overlap_locus_groups",
        "hgnc_overlap_locus_types",
        "hgnc_nearest_symbols",
        "hgnc_nearest_ids",
        "hgnc_nearest_distance_bp",
        "refgene_overlap_symbols",
        "refgene_nearest_symbols",
        "refgene_nearest_distance_bp",
        "cpg_island_overlap_names",
        "annotation_manifest",
    ]

    with open_text(args.out, "wt") as out:
        # Prefix with # so tools that expect BED-like comment headers can still skip it.
        out.write("#" + "\t".join(region_header + extra_header) + "\n")

        for region in regions:
            hgnc_hits = find_overlaps(region, hgnc)
            refgene_hits = find_overlaps(region, refgene)
            cpg_hits = find_overlaps(region, cpg)
            hgnc_near, hgnc_near_dist = find_nearest(region, hgnc)
            refgene_near, refgene_near_dist = find_nearest(region, refgene)

            extra = [
                split_values(field(x, 3) for x in hgnc_hits),
                split_values(field(x, 4) for x in hgnc_hits),
                split_values(field(x, 5) for x in hgnc_hits),
                split_values(field(x, 6) for x in hgnc_hits),
                split_values(field(x, 7) for x in hgnc_hits),
                split_values(field(x, 3) for x in hgnc_near),
                split_values(field(x, 4) for x in hgnc_near),
                hgnc_near_dist,
                split_values(field(x, 3) for x in refgene_hits),
                split_values(field(x, 3) for x in refgene_near),
                refgene_near_dist,
                split_values(field(x, 3) for x in cpg_hits),
                args.manifest,
            ]

            out.write("\t".join(region["fields"] + extra) + "\n")


if __name__ == "__main__":
    main()
