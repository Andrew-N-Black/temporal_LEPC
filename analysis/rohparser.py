#!/usr/bin/env python3
"""
rohparser.py -- parse one sample's bcftools-roh "RG" lines (already split
from the joint run's output, one file per sample -- see roh_analyses.sh
Step 6) into ROH length classes and F(ROH): >100kb-1Mb and >1Mb, matching
the classes bcftools roh's own documentation uses.

Vendored and corrected from Andrew-N-Black/LEPC-popgen's rohparser.py.
Fixes relative to that version:

  - Removed a division by a hardcoded cohort-size constant (num_sam) that
    was applied to a SINGLE sample's own ROH region counts/lengths before
    computing F(ROH). An individual's F(ROH) is (that individual's ROH
    segment lengths summed) / (genome length) -- there is no cohort-size
    term in that definition. As written, the division silently deflated
    every reported per-sample F(ROH) value by whatever num_sam happened to
    be hardcoded to (506, left over from a 433-sample cohort's copy of
    this script) -- e.g. a true F(ROH) of 0.02 would have been reported as
    0.02/506 ~= 0.00004.
  - Removed the hardcoded, cluster/project-specific absolute paths
    (path_to_directory, ref_index_file) in favor of CLI arguments. The
    previous version's path_to_directory was hardcoded to a DIFFERENT
    project's directory (GROUSE/nexus) than the one invoking it
    (GROUSE/old_vs_new) and was never actually patched by the caller (only
    ref_index_file was) -- since it's a hardcoded absolute path, `cd`-ing
    into the right directory before invoking this script didn't help.
  - Added --exclude-scaffolds so Z-linked (or other) scaffolds can be
    dropped from both the ROH sum and the genome-length denominator.
    Needed because females are hemizygous for Z, which biases F(ROH)
    upward if Z-linked "ROH" segments and Z's length are left in.

RG line format (bcftools roh, tab-separated):
  RG  Sample  Chromosome  Start  End  Length(bp)  NumMarkers  Quality

Usage:
  rohparser.py <sample>ROH.txt --fai ref.fna.fai \
      [--exclude-scaffolds NW_026294758.1,NW_026294813.1] \
      [--min-quality 30] [-o <sample>ROH.txt_results.txt]
"""
import argparse


def parse_args():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("roh_file", help="Per-sample bcftools-roh RG lines (e.g. <sample>ROH.txt)")
    ap.add_argument("--fai", required=True,
                    help="Reference .fai, for the genome-length denominator")
    ap.add_argument("--exclude-scaffolds", default="",
                    help="Comma-separated scaffold names to drop from both the ROH sum "
                         "and the genome length (e.g. the Z scaffolds)")
    ap.add_argument("--min-quality", type=float, default=30,
                    help="Minimum RG mean fwd-bwd phred quality to keep a region [30]")
    ap.add_argument("-o", "--output", default=None,
                    help="Output path [<roh_file>_results.txt]")
    return ap.parse_args()


def main():
    args = parse_args()
    excluded = {s.strip() for s in args.exclude_scaffolds.split(",") if s.strip()}
    out_path = args.output or f"{args.roh_file}_results.txt"

    print(f"Sample file: {args.roh_file}")

    # Genome length denominator, excluding any --exclude-scaffolds
    len_ref = 0.0
    excluded_len = 0.0
    with open(args.fai) as fh:
        for line in fh:
            field = line.rstrip("\n").split("\t")
            name, length = field[0], float(field[1])
            if name in excluded:
                excluded_len += length
            else:
                len_ref += length
    print(f"Genome length used for F(ROH), excluding {sorted(excluded) or 'nothing'}: "
          f"{len_ref:.0f} bp ({excluded_len:.0f} bp excluded)")

    # Sum ROH regions by length class, dropping excluded scaffolds and
    # regions below the quality threshold
    num_roh_100kb = num_roh_1mb = 0
    len_roh_100kb = len_roh_1mb = 0.0
    n_excluded_regions = 0
    n_low_qual = 0
    with open(args.roh_file) as fh:
        for line in fh:
            if not line.startswith("RG"):
                continue
            field = line.rstrip("\n").split("\t")
            chrom, length, quality = field[2], float(field[5]), float(field[7])
            if chrom in excluded:
                n_excluded_regions += 1
                continue
            if quality < args.min_quality:
                n_low_qual += 1
                continue
            if 100000 < length <= 1000000:
                num_roh_100kb += 1
                len_roh_100kb += length
            elif length > 1000000:
                num_roh_1mb += 1
                len_roh_1mb += length

    num_roh_tot = num_roh_100kb + num_roh_1mb
    len_roh_tot = len_roh_100kb + len_roh_1mb
    f_roh_1mb = len_roh_1mb / len_ref if len_ref else float("nan")
    f_roh_tot = len_roh_tot / len_ref if len_ref else float("nan")

    print(f"Regions on excluded scaffolds skipped: {n_excluded_regions}; "
          f"below quality {args.min_quality}: {n_low_qual}")

    result = (
        f"Number of ROH > 100kb and <= 1mb: {num_roh_100kb}\t"
        f"Length of ROH > 100kb and <= 1mb: {len_roh_100kb:.0f}\n"
        f"Number of ROH > 1mb: {num_roh_1mb}\t"
        f"Length of ROH > 1mb: {len_roh_1mb:.0f}\n"
        f"Number of ROH in total: {num_roh_tot}\t"
        f"Length of ROH in total: {len_roh_tot:.0f}\n"
        f"F(ROH) > 100kb: {f_roh_tot:.6f}\tF(ROH) > 1mb: {f_roh_1mb:.6f}\n"
    )
    print(result)
    with open(out_path, "w") as out:
        out.write(result)


if __name__ == "__main__":
    main()
