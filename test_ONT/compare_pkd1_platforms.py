#!/usr/bin/env python3
"""Compare Illumina and ONT evidence for two HG002 PKD1 SNVs (GRCh38).

Uses individual primary alignments; phasing requires one alignment to cover
both sites. Does not infer phase by joining paired-end mates or other reads.
"""

import argparse
import collections
import csv
import os
import sys

try:
    import pysam
except ImportError:
    sys.exit("Install pysam first: python -m pip install pysam")

SITES = ((2125707, "T", "C"), (2135583, "G", "T"))
LABELS = ("REF-REF", "ALT-REF", "REF-ALT", "ALT-ALT")


def observed_base(read, pos0, min_bq):
    if not (read.reference_start <= pos0 < read.reference_end):
        return None, "not_spanned"
    if read.query_sequence is None:
        return None, "no_sequence"
    for qpos, rpos in read.get_aligned_pairs(matches_only=False):
        if rpos != pos0:
            continue
        if qpos is None:
            return None, "deletion_or_skip"
        if read.query_qualities is not None and read.query_qualities[qpos] < min_bq:
            return None, "low_baseq"
        base = read.query_sequence[qpos].upper()
        return (base, "ok") if base in "ACGT" else (None, "ambiguous_base")
    return None, "not_aligned"


def analyze(path, platform, args):
    site_counts = [collections.Counter(), collections.Counter()]
    phase_counts = collections.Counter()
    exclusions = collections.Counter()
    phase_rows = []
    alignment_total = 0
    passed_alignment_filter = 0
    mode = "rc" if path.lower().endswith(".cram") else "rb"
    with pysam.AlignmentFile(path, mode, reference_filename=args.reference) as bam:
        if args.chrom not in bam.references:
            raise ValueError(f"{args.chrom} absent from {path}; check contig names/build")
        for read in bam.fetch(args.chrom, SITES[0][0] - 1, SITES[1][0]):
            alignment_total += 1
            if read.is_unmapped or read.is_secondary or read.is_supplementary:
                exclusions["non_primary"] += 1
                continue
            if read.is_duplicate or read.is_qcfail:
                exclusions["duplicate_or_qcfail"] += 1
                continue
            if read.mapping_quality < args.min_mapq:
                exclusions["low_mapq"] += 1
                continue
            passed_alignment_filter += 1
            calls = [observed_base(read, pos - 1, args.min_baseq)
                     for pos, _, _ in SITES]
            for i, (base, status) in enumerate(calls):
                if status != "ok":
                    site_counts[i]["excluded_" + status] += 1
                elif base == SITES[i][1]:
                    site_counts[i]["REF"] += 1
                elif base == SITES[i][2]:
                    site_counts[i]["ALT"] += 1
                else:
                    site_counts[i]["OTHER"] += 1

            if any(status != "ok" for _, status in calls):
                exclusions["no_two_high_quality_bases"] += 1
                continue
            bases = (calls[0][0], calls[1][0])
            if any(bases[i] not in SITES[i][1:] for i in range(2)):
                exclusions["other_base_on_spanning_read"] += 1
                continue
            label = ("ALT" if bases[0] == SITES[0][2] else "REF") + "-" + (
                "ALT" if bases[1] == SITES[1][2] else "REF")
            phase_counts[label] += 1
            phase_rows.append((platform, read.query_name, read.mapping_quality,
                               *bases, label))
    return dict(platform=platform, path=path, sites=site_counts,
                phases=phase_counts, exclusions=exclusions,
                phase_rows=phase_rows, total=alignment_total,
                filtered=passed_alignment_filter)


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--illumina", required=True, help="Illumina BAM/CRAM with index")
    p.add_argument("--ont", required=True, help="ONT BAM/CRAM with index")
    p.add_argument("--reference", required=True, help="CRAM-compatible GRCh38 FASTA")
    p.add_argument("--chrom", default="chr16")
    p.add_argument("--min-mapq", type=int, default=20)
    p.add_argument("--min-baseq", type=int, default=10)
    p.add_argument("--prefix", default="pkd1_comparison")
    args = p.parse_args()
    if args.min_mapq < 0 or args.min_baseq < 0:
        p.error("quality thresholds must be non-negative")
    results = []
    try:
        for platform, path in (("Illumina", args.illumina), ("ONT", args.ont)):
            results.append(analyze(path, platform, args))
    except (OSError, ValueError, KeyError) as exc:
        sys.exit(f"Error reading indexed alignments: {exc}")

    summary = args.prefix + ".sites.tsv"
    with open(summary, "w", newline="") as f:
        w = csv.writer(f, delimiter="\t")
        w.writerow(("platform", "chrom", "pos", "ref", "alt", "ref_reads",
                    "alt_reads", "other_reads", "usable_depth", "alt_fraction"))
        for result in results:
            for (pos, ref, alt), counts in zip(SITES, result["sites"]):
                depth = counts["REF"] + counts["ALT"] + counts["OTHER"]
                vaf = counts["ALT"] / depth if depth else None
                w.writerow((result["platform"], args.chrom, pos, ref, alt,
                            counts["REF"], counts["ALT"], counts["OTHER"], depth,
                            f"{vaf:.4f}" if vaf is not None else "NA"))

    phase_file = args.prefix + ".phase.tsv"
    with open(phase_file, "w", newline="") as f:
        w = csv.writer(f, delimiter="\t")
        w.writerow(("platform", "ref_ref", "alt_ref", "ref_alt", "alt_alt",
                    "cis_consistent", "trans_consistent", "total_informative"))
        for result in results:
            c = result["phases"]
            w.writerow((result["platform"], *(c[k] for k in LABELS),
                        c["REF-REF"] + c["ALT-ALT"],
                        c["ALT-REF"] + c["REF-ALT"], sum(c.values())))

    reads_file = args.prefix + ".reads.tsv"
    with open(reads_file, "w", newline="") as f:
        w = csv.writer(f, delimiter="\t")
        w.writerow(("platform", "read_name", "mapq", "base_2125707",
                    "base_2135583", "haplotype"))
        for result in results:
            w.writerows(result["phase_rows"])

    print(f"Sites: {args.chrom}:{SITES[0][0]} T>C; {args.chrom}:{SITES[1][0]} G>T")
    print(f"Filters: MAPQ >= {args.min_mapq}, base quality >= {args.min_baseq}; primary nonduplicate reads")
    for result in results:
        c = result["phases"]
        print(f"\n{result['platform']}: alignments in region={result['total']}; alignments after basic filters={result['filtered']}")
        for (pos, _, _), counts in zip(SITES, result["sites"]):
            print(f"  {pos}: REF={counts['REF']}, ALT={counts['ALT']}, OTHER={counts['OTHER']}")
        print("  Spanning reads:", ", ".join(f"{k}={c[k]}" for k in LABELS))
        print("  Exclusions:", dict(result["exclusions"]))
    print(f"\nCreated {summary}, {phase_file}, {reads_file}")
    print("Interpretation cautions: primary alignments are not independent molecules in every library;")
    print("MAPQ=60 is not proof of PKD1 specificity; no phase if no spanning reads.")
    print("Single-read evidence does not test mate-pair, trio or statistical phasing.")


if __name__ == "__main__":
    main()
