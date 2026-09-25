#!/usr/bin/env python3
"""
fst_temporal.py — Past vs Present LEPC FST, genome-wide and in sliding windows.

Estimators (both computed with scikit-allel):
  * Hudson FST (Hudson et al. 1992), combined as a ratio of averages
    (sum of per-site numerators / sum of per-site denominators), following
    the recommendations of Bhatia et al. 2013 (Genome Res 23:1514). This is
    the primary estimate: it is unbiased with respect to sample size and
    sample-size differences, and not dominated by rare variants.
  * Weir & Cockerham (1984) theta, also as a ratio of averages (sum a /
    sum(a+b+c)), for comparability with vcftools / earlier literature.

Genome-wide values get a block-jackknife SE and 95% CI (linked SNPs are not
independent). Windows use the same ratio-of-averages within each window.

Per-site filter: >= --mincalled genotyped individuals in EACH period and
polymorphic across all samples. No MAF filter (ratio of averages is robust
to rare variants, and dropping them biases FST upward).

Input: biallelic SNP VCF (e.g. output.subset.biallelic.snps.vcf.gz from the
temporal pipeline) and popmap <sample> <Past|Present> (old/new accepted).

Outputs (<outprefix>_*):
  fst_global.txt            genome-wide Hudson and WC FST with jackknife CIs
  fst_windows.tsv.gz        sliding-window FST (Hudson + WC), SNP counts, ZFST
  fst_windows_manhattan.png window Hudson FST along the genome
  fst_sites.tsv.gz          per-site components (only with --write_sites)

Note on interpretation: with two temporal samples, FST mostly measures drift.
Under pure drift from the Past sample, expected Hudson FST ~ t / (4 Ne).
High-FST windows are the tail of a drift distribution, not evidence of
selection by themselves.
"""

import argparse
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
    ap.add_argument("--vcf", required=True, help="Biallelic SNP VCF (.vcf or .vcf.gz)")
    ap.add_argument("--popmap", required=True, help="<sample> <Past|Present>, no header")
    ap.add_argument("--outprefix", default="fst_temporal")
    ap.add_argument("--window", type=int, default=50000, help="Window size, bp [50000]")
    ap.add_argument("--step", type=int, default=10000, help="Window step, bp [10000]")
    ap.add_argument("--min_snps", type=int, default=50,
                    help="Min SNPs for a window FST to be reported [50]")
    ap.add_argument("--mincalled", type=int, default=5,
                    help="Min genotyped individuals per period per site [5]")
    ap.add_argument("--jk_blocksize", type=int, default=5000000,
                    help="Block size (bp) for the genome-wide jackknife [5000000]")
    ap.add_argument("--min_plot_contig", type=int, default=1000000,
                    help="Only plot contigs at least this long, bp [1000000]")
    ap.add_argument("--generations", type=float, default=None,
                    help="Generations between samples; if given, report Ne = t/(4*FST)")
    ap.add_argument("--chunk_length", type=int, default=200000,
                    help="VCF records read per chunk [200000]")
    ap.add_argument("--write_sites", action="store_true",
                    help="Also write per-site FST components (large file)")
    return ap.parse_args()


def read_popmap(path):
    lookup = {"past": "past", "old": "past", "present": "present", "new": "present"}
    samples, pops = [], []
    with open(path) as fh:
        for line in fh:
            fields = line.strip().split()
            if not fields:
                continue
            if len(fields) < 2:
                sys.exit(f"ERROR: popmap line has <2 columns: {line!r}")
            lab = fields[1].strip().lower()
            if lab not in lookup:
                sys.exit(f"ERROR: unrecognized popmap label '{fields[1]}' (expected Past/Present)")
            samples.append(fields[0].strip())
            pops.append(lookup[lab])
    if len(set(samples)) != len(samples):
        sys.exit("ERROR: duplicate sample IDs in popmap")
    return samples, pops


def contig_lengths(headers):
    lens = {}
    pat = re.compile(r"##contig=<ID=([^,>]+).*?length=(\d+)")
    for h in headers.headers:
        m = pat.match(h)
        if m:
            lens[m.group(1)] = int(m.group(2))
    return lens


def ratio_jackknife(num_blocks, den_blocks):
    """Delete-one-block jackknife for a ratio of sums."""
    N, D = num_blocks.sum(), den_blocks.sum()
    est = N / D
    B = len(num_blocks)
    if B < 2:
        return est, np.nan, np.nan, np.nan
    loo = (N - num_blocks) / (D - den_blocks)
    se = np.sqrt((B - 1) / B * np.sum((loo - loo.mean()) ** 2))
    return est, se, est - 1.96 * se, est + 1.96 * se


def main():
    args = parse_args()
    if args.step <= 0 or args.window <= 0 or args.step > args.window:
        sys.exit("ERROR: need 0 < --step <= --window")

    out_dir = os.path.dirname(args.outprefix)
    if out_dir:
        os.makedirs(out_dir, exist_ok=True)
    for f in (args.vcf, args.popmap):
        if not os.path.exists(f):
            sys.exit(f"ERROR: input file not found: {f}")

    pm_samples, pm_pops = read_popmap(args.popmap)
    headers = allel.read_vcf_headers(args.vcf)
    vcf_samples = list(headers.samples)
    use = [s for s in pm_samples if s in vcf_samples]
    missing = [s for s in pm_samples if s not in vcf_samples]
    if missing:
        log(f"WARNING: popmap samples not in VCF: {', '.join(missing)}")
    pop_of = dict(zip(pm_samples, pm_pops))
    if not any(pop_of[s] == "past" for s in use) or not any(pop_of[s] == "present" for s in use):
        sys.exit("ERROR: need at least one Past and one Present sample in the VCF")
    clen = contig_lengths(headers)

    log(f"VCF: {args.vcf}")
    log(f"Samples used: {len(use)} "
        f"({sum(pop_of[s] == 'past' for s in use)} past, "
        f"{sum(pop_of[s] == 'present' for s in use)} present)")

    fields, samples_out, _, it = allel.iter_vcf_chunks(
        args.vcf, fields=["variants/CHROM", "variants/POS", "calldata/GT"],
        samples=use, chunk_length=args.chunk_length)
    samples_out = [s.decode() if isinstance(s, bytes) else str(s) for s in samples_out]
    idx_past = [i for i, s in enumerate(samples_out) if pop_of[s] == "past"]
    idx_pres = [i for i, s in enumerate(samples_out) if pop_of[s] == "present"]
    min_alleles = 2 * args.mincalled

    # Per-site accumulators (contig stored as integer code)
    contig_codes, contig_names = {}, []
    ch_l, pos_l, hn_l, hd_l, wn_l, wd_l = [], [], [], [], [], []
    n_read = 0
    for chunk, n, _, _ in it:
        n_read += n
        chrom = chunk["variants/CHROM"].astype(str)
        pos = chunk["variants/POS"].astype(np.int64)
        g = allel.GenotypeArray(chunk["calldata/GT"])

        ac1 = g.count_alleles(subpop=idx_past, max_allele=1)
        ac2 = g.count_alleles(subpop=idx_pres, max_allele=1)
        tot = ac1 + ac2
        keep = ((ac1.sum(axis=1) >= min_alleles) & (ac2.sum(axis=1) >= min_alleles) &
                (tot[:, 0] > 0) & (tot[:, 1] > 0))
        if not keep.any():
            continue

        h_num, h_den = allel.hudson_fst(ac1[keep], ac2[keep])
        a, b, c = allel.weir_cockerham_fst(g[keep], subpops=[idx_past, idx_pres], max_allele=1)
        w_num = np.nansum(a, axis=1)
        w_den = np.nansum(a + b + c, axis=1)

        for name in np.unique(chrom[keep]):
            if name not in contig_codes:
                contig_codes[name] = len(contig_names)
                contig_names.append(name)
        ch_l.append(np.array([contig_codes[x] for x in chrom[keep]], dtype=np.int32))
        pos_l.append(pos[keep])
        hn_l.append(h_num); hd_l.append(h_den)
        wn_l.append(w_num); wd_l.append(w_den)
        if len(ch_l) % 10 == 0:
            log(f"  read {n_read:,} records")

    if not ch_l:
        sys.exit("ERROR: no sites passed filters")
    ch = np.concatenate(ch_l); pos = np.concatenate(pos_l)
    hn = np.concatenate(hn_l); hd = np.concatenate(hd_l)
    wn = np.concatenate(wn_l); wd = np.concatenate(wd_l)
    del ch_l, pos_l, hn_l, hd_l, wn_l, wd_l

    ok = np.isfinite(hn) & np.isfinite(hd) & np.isfinite(wn) & np.isfinite(wd) & (hd > 0) & (wd > 0)
    ch, pos, hn, hd, wn, wd = ch[ok], pos[ok], hn[ok], hd[ok], wn[ok], wd[ok]
    log(f"Read {n_read:,} records; {len(pos):,} sites passed filters")

    # ---------------- Genome-wide estimates with block jackknife -------------
    block_key = ch.astype(np.int64) * 10**7 + pos // args.jk_blocksize
    _, blk = np.unique(block_key, return_inverse=True)
    nb = blk.max() + 1
    h_est, h_se, h_lo, h_hi = ratio_jackknife(np.bincount(blk, hn, nb), np.bincount(blk, hd, nb))
    w_est, w_se, w_lo, w_hi = ratio_jackknife(np.bincount(blk, wn, nb), np.bincount(blk, wd, nb))

    with open(f"{args.outprefix}_fst_global.txt", "w") as fh:
        def out(s=""):
            print(s); fh.write(s + "\n")
        out("Genome-wide FST, Past vs Present LEPC (ratio of averages)")
        out("==========================================================")
        out(f"VCF:              {args.vcf}")
        out(f"Samples:          {len(idx_past)} past, {len(idx_pres)} present")
        out(f"Sites used:       {len(pos):,} (>= {args.mincalled} called per period, polymorphic)")
        out(f"Jackknife blocks: {nb} x {args.jk_blocksize:,} bp")
        out("")
        out(f"Hudson FST (Bhatia et al. 2013): {h_est:.6f}  SE {h_se:.6f}  95% CI {h_lo:.6f} - {h_hi:.6f}")
        out(f"Weir & Cockerham theta:          {w_est:.6f}  SE {w_se:.6f}  95% CI {w_lo:.6f} - {w_hi:.6f}")
        if args.generations:
            t = args.generations
            fmt = lambda f: f"{t / (4 * f):.1f}" if f > 0 else "Inf"
            out("")
            out(f"Drift-only Ne implied by Hudson FST, Ne = t/(4*FST), t = {t:g}:")
            out(f"  Ne = {fmt(h_est)}  (95% CI {fmt(h_hi)} - {fmt(h_lo)})")
            out("  (approximate; assumes all divergence is drift from the Past sample)")

    # ---------------- Sliding windows ------------------------------------------
    log(f"Computing {args.window:,} bp windows, {args.step:,} bp step...")
    rows = []
    for code, name in enumerate(contig_names):
        sel = ch == code
        p = pos[sel]
        order = np.argsort(p, kind="stable")
        p = p[order]
        cs = {k: np.concatenate([[0.0], np.cumsum(v[sel][order])])
              for k, v in (("hn", hn), ("hd", hd), ("wn", wn), ("wd", wd))}
        L = clen.get(name, int(p.max()))
        starts = np.arange(1, max(L - args.window + 2, 2), args.step, dtype=np.int64)
        ends = starts + args.window - 1
        i0 = np.searchsorted(p, starts, side="left")
        i1 = np.searchsorted(p, ends, side="right")
        nsnp = i1 - i0
        with np.errstate(invalid="ignore", divide="ignore"):
            hud = (cs["hn"][i1] - cs["hn"][i0]) / (cs["hd"][i1] - cs["hd"][i0])
            wc = (cs["wn"][i1] - cs["wn"][i0]) / (cs["wd"][i1] - cs["wd"][i0])
        low = nsnp < args.min_snps
        hud[low] = np.nan; wc[low] = np.nan
        rows.append(pd.DataFrame({
            "chrom": name, "start": starts, "end": np.minimum(ends, L),
            "mid": (starts + np.minimum(ends, L)) // 2, "n_snps": nsnp,
            "hudson_fst": hud, "wc_fst": wc}))
    win = pd.concat(rows, ignore_index=True)

    good = win["hudson_fst"].notna()
    mu, sd = win.loc[good, "hudson_fst"].mean(), win.loc[good, "hudson_fst"].std()
    win["zfst_hudson"] = (win["hudson_fst"] - mu) / sd
    q999 = win.loc[good, "hudson_fst"].quantile(0.999)
    win["top_0.1pct"] = good & (win["hudson_fst"] >= q999)

    win.to_csv(f"{args.outprefix}_fst_windows.tsv.gz", sep="\t", index=False,
               float_format="%.6g", compression="gzip")
    log(f"Windows: {len(win):,} total, {good.sum():,} with >= {args.min_snps} SNPs; "
        f"mean Hudson window FST {mu:.5f}, 99.9th pct {q999:.5f}")

    # ---------------- Per-site output (optional) ---------------------------------
    if args.write_sites:
        log("Writing per-site components...")
        sites = pd.DataFrame({"chrom": np.array(contig_names)[ch], "pos": pos,
                              "hudson_num": hn, "hudson_den": hd,
                              "wc_num": wn, "wc_den": wd})
        sites.to_csv(f"{args.outprefix}_fst_sites.tsv.gz", sep="\t", index=False,
                     float_format="%.6g", compression="gzip")

    # ---------------- Manhattan-style plot -----------------------------------------
    plot_contigs = [c for c in contig_names
                    if clen.get(c, win.loc[win.chrom == c, "end"].max()) >= args.min_plot_contig]
    pw = win[win.chrom.isin(plot_contigs) & good].copy()
    if len(pw):
        offsets, ticks, labels, off = {}, [], [], 0
        for c in plot_contigs:
            L = clen.get(c, int(win.loc[win.chrom == c, "end"].max()))
            offsets[c] = off
            ticks.append(off + L / 2); labels.append(c)
            off += L
        pw["x"] = pw["mid"] + pw["chrom"].map(offsets)
        color_idx = pw["chrom"].map({c: i % 2 for i, c in enumerate(plot_contigs)})
        fig, ax = plt.subplots(figsize=(14, 4.5))
        ax.scatter(pw["x"], pw["hudson_fst"], s=1.5, linewidths=0,
                   c=np.where(color_idx == 0, "#3b6ea8", "#9aa7b8"), rasterized=True)
        ax.axhline(h_est, color="black", lw=0.8, ls="--", label=f"genome-wide {h_est:.4f}")
        ax.axhline(q999, color="#c0392b", lw=0.8, ls=":", label=f"99.9th pct {q999:.4f}")
        ax.set_xlim(0, off)
        if len(plot_contigs) <= 60:
            ax.set_xticks(ticks)
            ax.set_xticklabels(labels, rotation=90, fontsize=6)
        else:
            ax.set_xticks([])
            ax.set_xlabel(f"{len(plot_contigs)} contigs >= {args.min_plot_contig:,} bp")
        ax.set_ylabel("Hudson FST (window)")
        ax.set_title(f"Past vs Present LEPC: {args.window // 1000} kb windows, "
                     f"{args.step // 1000} kb step")
        ax.legend(loc="upper right", fontsize=8, frameon=False)
        fig.tight_layout()
        fig.savefig(f"{args.outprefix}_fst_windows_manhattan.png", dpi=300)
        plt.close(fig)

    log(f"Done. Outputs written with prefix: {args.outprefix}")


if __name__ == "__main__":
    main()
