#!/usr/bin/env python3
"""
Convert a biallelic-SNP VCF (already filtered and thinned) into a GENEPOP file
for NeEstimator's temporal method.

Samples are written as two "Pop" blocks in time order -- Past first, then
Present -- because NeEstimator reads samples in file order and assigns the
generation set (e.g. 0 5) to them in that order.

Genotype coding: REF = 01, ALT = 02, missing = 0000 (2-digit alleles).

Usage:
  bcftools query -l filtered.vcf.gz > samples.txt   # (only needed for checking)
  python3 vcf2genepop.py --vcf filtered.vcf.gz --popmap popmap_unrelated.txt \
          --out lepc_temporal.gen [--order Past,Present]

Requires bcftools on PATH (reads the VCF through `bcftools query`).
"""
import argparse
import subprocess
import sys

CODE = {"0": "01", "1": "02"}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--vcf", required=True)
    ap.add_argument("--popmap", required=True, help="<sample> <group>, no header")
    ap.add_argument("--out", required=True)
    ap.add_argument("--order", default="Past,Present",
                    help="group order = sampling order (generation 0 first)")
    a = ap.parse_args()

    order = a.order.split(",")
    group = {}
    with open(a.popmap) as fh:
        for line in fh:
            f = line.split()
            if len(f) >= 2:
                group[f[0]] = f[1]

    vcf_samples = subprocess.run(["bcftools", "query", "-l", a.vcf],
                                 check=True, capture_output=True,
                                 text=True).stdout.split()
    missing = [s for s in vcf_samples if s not in group]
    if missing:
        sys.exit(f"ERROR: samples in VCF but not in popmap: {missing}")
    bad = sorted({group[s] for s in vcf_samples} - set(order))
    if bad:
        sys.exit(f"ERROR: popmap groups {bad} not in --order {order}")

    loci = []
    geno = {s: [] for s in vcf_samples}
    q = subprocess.Popen(["bcftools", "query", "-f", "%CHROM:%POS[\t%GT]\n", a.vcf],
                         stdout=subprocess.PIPE, text=True)
    for line in q.stdout:
        f = line.rstrip("\n").split("\t")
        loci.append(f[0])
        for s, gt in zip(vcf_samples, f[1:]):
            al = gt.replace("|", "/").split("/")
            if len(al) == 2 and al[0] in CODE and al[1] in CODE:
                geno[s].append(CODE[al[0]] + CODE[al[1]])
            else:
                geno[s].append("0000")
    if q.wait() != 0:
        sys.exit("ERROR: bcftools query failed")

    with open(a.out, "w") as out:
        out.write(f"LEPC temporal: {len(loci)} SNPs, groups {','.join(order)}\n")
        for loc in loci:
            out.write(loc + "\n")
        for g in order:
            members = [s for s in vcf_samples if group[s] == g]
            out.write("Pop\n")
            for s in members:
                out.write(f"{s} , " + " ".join(geno[s]) + "\n")
            print(f"{g}: {len(members)} individuals", file=sys.stderr)
    print(f"Wrote {a.out}: {len(loci)} loci", file=sys.stderr)


if __name__ == "__main__":
    main()
