#!/usr/bin/env python3
"""
sex_from_vcf.py — confirm Z-linked scaffolds and infer bird sex from a VCF.

Birds are ZW (female) / ZZ (male). Females carry one Z, so Z-linked sites
should be (nearly) never heterozygous in females. This script:

  1. Computes each sample's heterozygosity on candidate Z scaffolds relative
     to autosomes (ratio_Z_auto). Females ~0; males clearly > 0 (Z diversity
     is lower than autosomal, so males are often ~0.3-0.8, not 1).
     Calls: female < --female_max, male > --male_min, otherwise ambiguous.
  2. Scores EVERY scaffold >= --min_len for female hemizygosity:
       score = mean(female het / autosomal het) / mean(male het / autosomal het)
     Z scaffolds score ~0, autosomes ~1. This catches Z scaffolds the chicken
     synteny missed and flags false positives. Needs >= 1 female and 1 male.

Genotype errors (e.g. in historical samples) raise female Z heterozygosity
above zero, so check the per-sample table rather than trusting calls blindly.

Outputs (<outprefix>_*):
  sample_sex.tsv            per-sample heterozygosity and sex call
  scaffold_hemizygosity.tsv per-scaffold female/male heterozygosity score
  sample_sex.png            ratio_Z_auto per sample
"""

import argparse
import faulthandler
import os
import re
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
    ap.add_argument("--vcf", required=True)
    ap.add_argument("--z_scaffolds", required=True,
                    help="Comma-separated candidate Z scaffolds, or a file with one per line")
    ap.add_argument("--popmap", default=None, help="Optional <sample> <Past|Present>")
    ap.add_argument("--exclude", default="",
                    help="Comma-separated scaffolds to ignore entirely (e.g. candidate W)")
    ap.add_argument("--min_len", type=int, default=1000000,
                    help="Min scaffold length for per-scaffold scores [1000000]")
    ap.add_argument("--min_sites", type=int, default=1000,
                    help="Min polymorphic sites for a per-scaffold score [1000]")
    ap.add_argument("--female_max", type=float, default=0.15)
    ap.add_argument("--male_min", type=float, default=0.25)
    ap.add_argument("--outprefix", default="sexcheck")
    ap.add_argument("--chunk_length", type=int, default=200000)
    return ap.parse_args()


def read_list(arg):
    if not arg:
        return set()
    if os.path.exists(arg):
        with open(arg) as fh:
            return {l.strip().split()[0] for l in fh if l.strip()}
    return {x.strip() for x in arg.split(",") if x.strip()}


def main():
    faulthandler.enable(all_threads=True)   # print a traceback on SIGBUS/SIGSEGV
    args = parse_args()
    if not os.path.exists(args.vcf):
        sys.exit(f"ERROR: VCF not found: {args.vcf}")
    out_dir = os.path.dirname(args.outprefix)
    if out_dir:
        os.makedirs(out_dir, exist_ok=True)

    zset = read_list(args.z_scaffolds)
    exclude = read_list(args.exclude)
    if not zset:
        sys.exit("ERROR: no candidate Z scaffolds given")

    headers = allel.read_vcf_headers(args.vcf)
    clen = {}
    for h in headers.headers:
        m = re.match(r"##contig=<ID=([^,>]+).*?length=(\d+)", h)
        if m:
            clen[m.group(1)] = int(m.group(2))

    period = {}
    if args.popmap:
        with open(args.popmap) as fh:
            for l in fh:
                f = l.strip().split()
                if len(f) >= 2:
                    period[f[0]] = f[1]

    _, samples, _, it = allel.iter_vcf_chunks(
        args.vcf, fields=["variants/CHROM", "calldata/GT"], chunk_length=args.chunk_length)
    names = [s.decode() if isinstance(s, bytes) else str(s) for s in samples]
    if args.popmap:
        names_keep = [i for i, s in enumerate(names) if s in period]
    else:
        names_keep = list(range(len(names)))
    names = [names[i] for i in names_keep]
    ns = len(names)

    scaf_codes, scaf_names = {}, []
    het_by, called_by, sites_by = [], [], []   # lists indexed by scaffold code
    n_read = 0
    for chunk, nrec, _, _ in it:
        n_read += nrec
        chrom = chunk["variants/CHROM"].astype(str)
        gt = chunk["calldata/GT"][:, names_keep, :]
        a0, a1 = gt[..., 0], gt[..., 1]
        called = (a0 >= 0) & (a1 >= 0)
        het = called & (a0 != a1)
        alt = np.where(called, (a0 > 0).astype(int) + (a1 > 0).astype(int), 0).sum(1)
        poly = (alt > 0) & (alt < 2 * called.sum(1))
        for name in np.unique(chrom):
            if name in exclude:
                continue
            sel = (chrom == name) & poly
            if not sel.any():
                continue
            if name not in scaf_codes:
                scaf_codes[name] = len(scaf_names)
                scaf_names.append(name)
                het_by.append(np.zeros(ns)); called_by.append(np.zeros(ns)); sites_by.append(0)
            k = scaf_codes[name]
            het_by[k] += het[sel].sum(0)
            called_by[k] += called[sel].sum(0)
            sites_by[k] += int(sel.sum())
        if n_read % (args.chunk_length * 10) == 0:
            log(f"  read {n_read:,} records")
    log(f"Read {n_read:,} records across {len(scaf_names)} scaffolds")

    log("Summarizing per-sample heterozygosity...")
    het_by = np.array(het_by); called_by = np.array(called_by); sites_by = np.array(sites_by)
    np.savez(f"{args.outprefix}_counts.npz", scaffolds=np.array(scaf_names), samples=np.array(names),
             het=het_by, called=called_by, sites=sites_by)
    is_z = np.array([s in zset for s in scaf_names])
    missing_z = zset - set(scaf_names)
    if missing_z:
        log(f"NOTE: candidate Z scaffolds with no polymorphic sites in VCF: {', '.join(sorted(missing_z))}")
    if not is_z.any():
        sys.exit("ERROR: none of the candidate Z scaffolds occur in the VCF")

    het_auto = het_by[~is_z].sum(0) / np.maximum(called_by[~is_z].sum(0), 1)
    het_z = het_by[is_z].sum(0) / np.maximum(called_by[is_z].sum(0), 1)
    ratio = het_z / het_auto
    sex = np.where(ratio < args.female_max, "female",
                   np.where(ratio > args.male_min, "male", "ambiguous"))

    samp = pd.DataFrame({"sample": names,
                         "period": [period.get(s, "NA") for s in names],
                         "het_auto": het_auto, "het_Z": het_z,
                         "ratio_Z_auto": ratio, "sex_call": sex})
    samp.to_csv(f"{args.outprefix}_sample_sex.tsv", sep="\t", index=False, float_format="%.5f")
    log(f"Z sites used: {sites_by[is_z].sum():,} on {is_z.sum()} scaffolds")
    print(samp.to_string(index=False, float_format=lambda x: f"{x:.4f}"), flush=True)
    if args.popmap:
        print(pd.crosstab(samp.period, samp.sex_call).to_string(), flush=True)

    # Per-scaffold hemizygosity score
    log("Scoring scaffolds...")
    fem = sex == "female"; mal = sex == "male"
    rows = []
    for k, s in enumerate(scaf_names):
        L = clen.get(s, np.nan)
        score = np.nan
        if fem.any() and mal.any() and sites_by[k] >= args.min_sites and (np.isnan(L) or L >= args.min_len):
            h = het_by[k] / np.maximum(called_by[k], 1) / het_auto
            m_mean = h[mal].mean()
            score = h[fem].mean() / m_mean if m_mean > 0 else np.nan
        rows.append({"scaffold": s, "length": L, "n_sites": int(sites_by[k]),
                     "candidate_Z": bool(is_z[k]), "female_male_het_score": score})
    sc = pd.DataFrame(rows)
    sc["Z_like"] = sc["female_male_het_score"] < 0.3
    sc = sc.sort_values(["Z_like", "length"], ascending=[False, False])
    sc.to_csv(f"{args.outprefix}_scaffold_hemizygosity.tsv", sep="\t", index=False,
              float_format="%.4f")

    if not (fem.any() and mal.any()):
        log("NOTE: need at least one female and one male to score scaffolds; "
            "only the per-sample table is informative.")
    else:
        agree = sc[sc.female_male_het_score.notna()]
        new_z = agree[agree.Z_like & ~agree.candidate_Z]
        not_z = agree[~agree.Z_like & agree.candidate_Z]
        log(f"Scored {len(agree)} scaffolds: {int(agree.Z_like.sum())} Z-like")
        if len(new_z):
            log("Z-like scaffolds NOT in the synteny list (check these):")
            print(new_z.to_string(index=False), flush=True)
        if len(not_z):
            log("Synteny Z candidates that do NOT look hemizygous (check these):")
            print(not_z.to_string(index=False), flush=True)
        confirmed = agree[agree.Z_like].scaffold.tolist()
        with open(f"{args.outprefix}_Z_confirmed.txt", "w") as fh:
            fh.write("\n".join(confirmed) + ("\n" if confirmed else ""))

    log("Plotting...")
    fig, ax = plt.subplots(figsize=(7, 3.8))
    order = np.argsort(ratio)
    cols = {"female": "#b5651d", "male": "#2e6f9e", "ambiguous": "#7f7f7f"}
    ax.bar(range(ns), ratio[order], color=[cols[x] for x in sex[order]])
    ax.set_xticks(range(ns))
    ax.set_xticklabels([names[i] for i in order], rotation=90, fontsize=6)
    ax.axhline(args.female_max, ls=":", color="grey", lw=0.8)
    ax.axhline(args.male_min, ls=":", color="grey", lw=0.8)
    ax.set_ylabel("Z / autosomal heterozygosity")
    ax.set_title("Sex inference from Z heterozygosity")
    fig.tight_layout()
    fig.savefig(f"{args.outprefix}_sample_sex.png", dpi=300)
    plt.close(fig)
    log(f"Done. Outputs written with prefix: {args.outprefix}")


if __name__ == "__main__":
    main()
