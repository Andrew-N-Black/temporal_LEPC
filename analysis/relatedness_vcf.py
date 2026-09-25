#!/usr/bin/env python3
"""
relatedness_vcf.py — pairwise relatedness among Past and Present LEPC samples
from called genotypes in a VCF, using allele-frequency-free statistics.

Statistics (Manichaikul et al. 2010; Waples et al. 2019, Mol Ecol 28:35):
  KING-robust kinship = (HETHET - 2*IBS0) / (Het_i + Het_j)
  R0                  = IBS0 / HETHET
  R1                  = HETHET / (IBS0 + IBS1)
where, over sites genotyped in both samples,
  HETHET = both heterozygous; IBS0 = opposite homozygotes;
  IBS1   = one heterozygous, one homozygous; Het_i = heterozygous sites in i.

None of these need allele frequencies, which matters here: 10 birds per
period is too few to estimate frequencies well, and Past and Present
frequencies differ. KING-robust is also robust to population structure.

Kinship degree cut-offs (Manichaikul et al. 2010):
  > 0.354 duplicate/identical, 0.177-0.354 1st degree, 0.0884-0.177 2nd,
  0.0442-0.0884 3rd, < 0.0442 unrelated.
R0 near 0 with 1st-degree kinship indicates parent-offspring; R0 clearly
above 0 indicates full siblings. The R1-vs-KING and R0-vs-KING plots follow
Waples et al. 2019.

Caveat: called genotypes from low-coverage or damaged (historical) samples
under-call heterozygotes, which pushes KING negative and inflates IBS0 for
those samples. Check the per-sample heterozygosity table: if Past samples
have clearly lower heterozygosity than Present ones, trust the
genotype-likelihood (NgsRelate) results for Past pairs instead.

Outputs (<outprefix>_*):
  pairs.tsv          all pairs: n_sites, KING, R0, R1, degree, pair type
  samples.tsv        per-sample heterozygosity and missingness
  R1_vs_KING.png, R0_vs_KING.png
"""

import argparse
import os
import sys
import time

import numpy as np
import pandas as pd
import allel

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt


def log(msg):
    print(f"[{time.strftime('%H:%M:%S')}] {msg}", flush=True)


def parse_args():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--vcf", required=True, help="Biallelic SNP VCF (.vcf or .vcf.gz)")
    ap.add_argument("--popmap", required=True, help="<sample> <Past|Present>, no header")
    ap.add_argument("--outprefix", default="relatedness")
    ap.add_argument("--exclude_chroms", default="",
                    help="Comma-separated contigs to skip (e.g. Z and W sex chromosomes)")
    ap.add_argument("--chunk_length", type=int, default=200000)
    return ap.parse_args()


def read_popmap(path):
    lookup = {"past": "Past", "old": "Past", "present": "Present", "new": "Present"}
    samples, pops = [], []
    with open(path) as fh:
        for line in fh:
            f = line.strip().split()
            if not f:
                continue
            lab = f[1].strip().lower() if len(f) > 1 else ""
            if lab not in lookup:
                sys.exit(f"ERROR: bad popmap line: {line!r}")
            samples.append(f[0].strip())
            pops.append(lookup[lab])
    if len(set(samples)) != len(samples):
        sys.exit("ERROR: duplicate sample IDs in popmap")
    return samples, pops


def degree(k):
    if not np.isfinite(k):
        return "NA"
    if k > 0.354:
        return "duplicate/identical"
    if k > 0.177:
        return "1st degree"
    if k > 0.0884:
        return "2nd degree"
    if k > 0.0442:
        return "3rd degree"
    return "unrelated"


def main():
    args = parse_args()
    for f in (args.vcf, args.popmap):
        if not os.path.exists(f):
            sys.exit(f"ERROR: input file not found: {f}")
    out_dir = os.path.dirname(args.outprefix)
    if out_dir:
        os.makedirs(out_dir, exist_ok=True)
    exclude = {c.strip() for c in args.exclude_chroms.split(",") if c.strip()}

    pm_samples, pm_pops = read_popmap(args.popmap)
    pop_of = dict(zip(pm_samples, pm_pops))
    vcf_samples = [str(s) for s in allel.read_vcf_headers(args.vcf).samples]
    use = [s for s in pm_samples if s in vcf_samples]
    missing = [s for s in pm_samples if s not in vcf_samples]
    if missing:
        log(f"WARNING: popmap samples not in VCF: {', '.join(missing)}")
    if len(use) < 2:
        sys.exit("ERROR: fewer than 2 popmap samples found in the VCF")

    _, samples_out, _, it = allel.iter_vcf_chunks(
        args.vcf, fields=["variants/CHROM", "calldata/GT"],
        samples=use, chunk_length=args.chunk_length)
    names = [s.decode() if isinstance(s, bytes) else str(s) for s in samples_out]
    n = len(names)
    log(f"{n} samples; excluding contigs: {', '.join(sorted(exclude)) or '(none)'}")

    HH = np.zeros((n, n))      # both het
    OPP = np.zeros((n, n))     # opposite homozygotes (IBS0)
    HHOM = np.zeros((n, n))    # one het, one hom (IBS1)
    HC = np.zeros((n, n))      # [i, j]: i het and j called
    CC = np.zeros((n, n))      # both called
    het_tot = np.zeros(n); called_tot = np.zeros(n)
    n_read = n_used = 0

    for chunk, nrec, _, _ in it:
        n_read += nrec
        chrom = chunk["variants/CHROM"].astype(str)
        gt = chunk["calldata/GT"]                              # (sites, samples, 2)
        keep = ~np.isin(chrom, list(exclude)) if exclude else np.ones(nrec, bool)
        gt = gt[keep]
        a0, a1 = gt[..., 0], gt[..., 1]
        called = (a0 >= 0) & (a1 >= 0) & (a0 <= 1) & (a1 <= 1)
        dos = np.where(called, a0 + a1, -1)
        # keep sites polymorphic among called genotypes
        alt = np.where(called, dos, 0).sum(1); nc = called.sum(1)
        poly = (alt > 0) & (alt < 2 * nc)
        if not poly.any():
            continue
        dos, called = dos[poly], called[poly]
        n_used += dos.shape[0]

        C = called.astype(np.float32)
        H = (dos == 1).astype(np.float32)
        A = (dos == 0).astype(np.float32)
        B = (dos == 2).astype(np.float32)
        HOM = A + B
        HH += H.T @ H
        OPP += A.T @ B + B.T @ A
        HHOM += H.T @ HOM + HOM.T @ H
        HC += H.T @ C
        CC += C.T @ C
        het_tot += H.sum(0); called_tot += C.sum(0)
        if n_read % (args.chunk_length * 10) == 0:
            log(f"  read {n_read:,} records")

    log(f"Read {n_read:,} records; {n_used:,} polymorphic autosomal sites used")
    if n_used == 0:
        sys.exit("ERROR: no usable sites")

    # Per-sample summary
    samp = pd.DataFrame({
        "sample": names, "period": [pop_of[s] for s in names],
        "sites_called": called_tot.astype(int),
        "missing_frac": 1 - called_tot / n_used,
        "het_rate_polymorphic_sites": het_tot / np.maximum(called_tot, 1)})
    samp.to_csv(f"{args.outprefix}_samples.tsv", sep="\t", index=False, float_format="%.5f")
    log("Per-sample heterozygosity (at polymorphic sites) by period:")
    for per, grp in samp.groupby("period"):
        log(f"  {per}: mean {grp.het_rate_polymorphic_sites.mean():.4f} "
            f"(range {grp.het_rate_polymorphic_sites.min():.4f}-"
            f"{grp.het_rate_polymorphic_sites.max():.4f}), "
            f"mean missing {grp.missing_frac.mean():.3f}")

    rows = []
    for i in range(n):
        for j in range(i + 1, n):
            hethet, ibs0, ibs1 = HH[i, j], OPP[i, j], HHOM[i, j]
            den = HC[i, j] + HC[j, i]
            king = (hethet - 2 * ibs0) / den if den > 0 else np.nan
            r0 = ibs0 / hethet if hethet > 0 else np.nan
            r1 = hethet / (ibs0 + ibs1) if (ibs0 + ibs1) > 0 else np.nan
            p1, p2 = pop_of[names[i]], pop_of[names[j]]
            ptype = f"{p1}-{p2}" if p1 == p2 else "Past-Present"
            rows.append({"id1": names[i], "id2": names[j], "pair_type": ptype,
                         "n_sites": int(CC[i, j]), "HETHET": int(hethet),
                         "IBS0": int(ibs0), "IBS1": int(ibs1),
                         "KING": king, "R0": r0, "R1": r1, "degree": degree(king)})
    pairs = pd.DataFrame(rows).sort_values("KING", ascending=False)
    pairs.to_csv(f"{args.outprefix}_pairs.tsv", sep="\t", index=False, float_format="%.5f")

    log("Pairs by inferred degree:")
    print(pd.crosstab(pairs.pair_type, pairs.degree).to_string(), flush=True)
    rel = pairs[pairs.degree != "unrelated"]
    if len(rel):
        log("Related pairs (KING > 0.0442):")
        print(rel[["id1", "id2", "pair_type", "KING", "R0", "R1", "degree"]]
              .to_string(index=False, float_format=lambda x: f"{x:.4f}"), flush=True)
    log("KING by pair type (median, min, max):")
    for pt, grp in pairs.groupby("pair_type"):
        log(f"  {pt}: {grp.KING.median():.4f}, {grp.KING.min():.4f}, {grp.KING.max():.4f}")

    colors = {"Past-Past": "#b5651d", "Present-Present": "#2e6f9e", "Past-Present": "#7f7f7f"}
    for ycol in ("R1", "R0"):
        fig, ax = plt.subplots(figsize=(6, 5))
        for pt, grp in pairs.groupby("pair_type"):
            ax.scatter(grp[ycol], grp.KING, s=22, alpha=0.8,
                       color=colors.get(pt, "black"), label=pt, edgecolors="none")
        for cut in (0.0442, 0.0884, 0.177, 0.354):
            ax.axhline(cut, color="grey", lw=0.6, ls=":")
        ax.set_xlabel(ycol); ax.set_ylabel("KING-robust kinship")
        ax.set_title(f"{ycol} vs KING (Waples et al. 2019)")
        ax.legend(frameon=False, fontsize=8)
        fig.tight_layout()
        fig.savefig(f"{args.outprefix}_{ycol}_vs_KING.png", dpi=300)
        plt.close(fig)

    log(f"Done. Outputs written with prefix: {args.outprefix}")


if __name__ == "__main__":
    main()
