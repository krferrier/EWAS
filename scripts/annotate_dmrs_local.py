#!/usr/bin/env python3
"""Annotate comb-p DMR regions using local refGene, CpG island, and HGNC BED files.

The output keeps all original comb-p region columns and appends local annotation
columns. The first three columns must be chrom/start/end.
"""

import argparse
import gzip
import re
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


def signed_interval_distance(region: dict, ann: dict) -> int:
    """Return signed bp distance from region to annotation.

    0 means overlap. Negative means the annotation is upstream/left of the
    region. Positive means the annotation is downstream/right of the region.
    This matches the sign convention seen in comb-p/cruzdb annotation output.
    """
    if ann["end"] > region["start"] and ann["start"] < region["end"]:
        return 0
    if ann["end"] <= region["start"]:
        return ann["end"] - region["start"]
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

    best_abs = None
    best = []
    best_signed = []

    for ann in anns:
        dist = signed_interval_distance(region, ann)
        abs_dist = abs(dist)
        if best_abs is None or abs_dist < best_abs:
            best_abs = abs_dist
            best = [ann]
            best_signed = [dist]
        elif abs_dist == best_abs:
            best.append(ann)
            best_signed.append(dist)

    # Usually this is a single value. If there are ties, preserve all signed distances.
    dist_string = ";".join(str(x) for x in sorted(set(best_signed)))
    return best, dist_string


def cpg_feature_from_distance(distance_string: str, shore_dist: int = 3000) -> str:
    """Mimic cruzdb.models.cpgIslandExt.distance().

    The cruzdb implementation returns:
      * (0, "island") for overlapping/touching CpG islands
      * (distance, "shore") when distance <= shore_dist
      * (distance, "") otherwise

    This function accepts signed distances from this script and classifies using
    the absolute distance. For tied nearest islands, distances may be semicolon
    separated.
    """
    if distance_string in ["", "NA", None]:
        return ""

    features = []
    for token in str(distance_string).split(";"):
        token = token.strip()
        if not token:
            continue
        try:
            dist = abs(int(token))
        except ValueError:
            continue

        if dist == 0:
            features.append("island")
        elif dist <= shore_dist:
            features.append("shore")
        else:
            features.append("")

    if not features:
        return ""

    # Match comb-p/cruzdb style reasonably closely while avoiding redundant
    # repeated labels in tie cases. Preserve empty if all nearest islands are
    # outside the shore window.
    non_empty = sorted({x for x in features if x})
    return ";".join(non_empty) if non_empty else ""


def field(rec: dict, idx: int) -> str:
    return rec["fields"][idx] if len(rec["fields"]) > idx else ""


def validate_cpg_bed(cpg_by_chrom: Dict[str, List[dict]], path: str) -> None:
    """Fail fast if the cached cpgIslandExt BED was built with whitespace splitting.

    UCSC cpgIslandExt names look like "CpG: 47". A previous version of the
    Snakemake rule used awk's default whitespace separator and converted these
    names to just "CpG:". That creates a syntactically valid BED file but
    incorrect annotations, so detect it explicitly.
    """
    names = []
    for records in cpg_by_chrom.values():
        for rec in records:
            names.append(field(rec, 3))
            if len(names) >= 200:
                break
        if len(names) >= 200:
            break

    if not names:
        raise RuntimeError(
            "No CpG island records were read from {}".format(path)
        )

    good = sum(1 for name in names if re.match(r"^CpG:\s+\d+$", name or ""))
    truncated = sum(1 for name in names if (name or "").strip() == "CpG:")

    if good == 0 or truncated > good:
        example = ", ".join(repr(x) for x in names[:5])
        raise RuntimeError(
            "The cached cpgIslandExt BED appears to have truncated CpG island names. "
            "Expected values like 'CpG: 47', but saw examples: {}. "
            "Rebuild the cache after changing the cpgIslandExt awk command to use "
            "FS=OFS='\t', then delete cpgIslandExt.bed.gz and rerun.".format(example)
        )


def validate_hgnc_bed(hgnc_by_chrom: Dict[str, List[dict]], path: str) -> None:
    """Basic sanity check for the HGNC BED-like table.

    HGNC locus types can legitimately contain spaces, for example
    "RNA, long non-coding". They should still occupy one tab-delimited column.
    """
    checked = 0
    bad = []
    for records in hgnc_by_chrom.values():
        for rec in records:
            checked += 1
            if len(rec["fields"]) < 9:
                bad.append(rec["fields"])
            if checked >= 200:
                break
        if checked >= 200:
            break

    if bad:
        raise RuntimeError(
            "The HGNC BED-like file {} has rows with fewer than 9 tab-delimited columns. "
            "First bad row: {}".format(path, bad[0])
        )


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

    validate_hgnc_bed(hgnc, args.hgnc_bed)
    validate_cpg_bed(cpg, args.cpg_bed)

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
        # comb-p-like CpG island nearest annotation columns
        "cpgIslandExt_name",
        "cpgIslandExt_distance",
        "cpgIslandExt_feature",
        "annotation_manifest",
    ]

    with open_text(args.out, "wt") as out:
        # Prefix with # so tools that expect BED-like comment headers can still skip it.
        out.write("#" + "\t".join(region_header + extra_header) + "\n")

        for region in regions:
            hgnc_hits = find_overlaps(region, hgnc)
            refgene_hits = find_overlaps(region, refgene)
            hgnc_near, hgnc_near_dist = find_nearest(region, hgnc)
            refgene_near, refgene_near_dist = find_nearest(region, refgene)
            cpg_near, cpg_near_dist = find_nearest(region, cpg)

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
                split_values(field(x, 3) for x in cpg_near),
                cpg_near_dist,
                cpg_feature_from_distance(cpg_near_dist),
                args.manifest,
            ]

            out.write("\t".join(region["fields"] + extra) + "\n")


if __name__ == "__main__":
    main()