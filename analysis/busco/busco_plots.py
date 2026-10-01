#!/usr/bin/env python3
# =============================================================================
# 04_busco_plots.py — BUSCO completeness, karyotype and synteny plots
#
# Driven by 04_busco_plots.sh. Uses BUSCO-Plot-Py (buscoplotpy):
#   https://github.com/lorenzo-arcioni/BUSCO-Plot-Py
#
# Reads, per assembly:
#   <busco_dir>/<SPECIES>_<SAMPLE>_<HAP>/short_summary.specific.*.json
#   <busco_dir>/<SPECIES>_<SAMPLE>_<HAP>/run_*/full_table.tsv
#   <final_dir>/<SPECIES>_<SAMPLE>_<HAP>.pseudo_chr.fasta.fai
#
# Writes:
#   completeness/busco_barplot_<SPECIES>_completeness.png  one panel per species
#   completeness/busco_barplot_ALL_completeness.png        all assemblies together
#   karyotype/<SPECIES>_<SAMPLE>_<HAP>.png     one per assembly
#   synteny/<A>__vs__<B>.png                   one per requested pair
# Every PNG is accompanied by an .svg of the same name (vector, for figures
# that will be scaled or edited); pass --no-svg to suppress.
#
# NOTES ON buscoplotpy's EXPECTATIONS (established by reading its source,
# because they are not in the README):
#   * organism_busco_barplot() needs 'group', 'organism' and 'version' columns
#     that load_json_summary() does NOT create — we add them here.
#   * the karyotype DataFrame needs 'chr', 'end', 'sequence', 'organism', and
#     synteny additionally reads karyotype['color'][0] by LABEL, so the index
#     must be a clean RangeIndex starting at 0.
#   * karyoplot() matches karyotype['chr'] against fulltable['sequence'].
#   * synteny keeps only status == 'Complete' rows and joins on 'busco_id'.
# =============================================================================

import argparse
import glob
import json
import os
import re
import sys

import matplotlib
matplotlib.use("Agg")          # headless compute node — must precede pyplot
import matplotlib.pyplot as plt
import pandas as pd

from buscoplotpy.utils.load_busco_fulltable import load_busco_fulltable
from buscoplotpy.utils.load_json_summary import load_json_summary
from buscoplotpy.graphics.organism_busco_barplot import organism_busco_barplot
from buscoplotpy.graphics.karyoplot import karyoplot
from buscoplotpy.graphics.synteny import horizontal_synteny_plot, vertical_synteny_plot

import contextlib


@contextlib.contextmanager
def also_save_svg(enabled=True):
    """
    Write an .svg beside every .png the plotting functions produce.

    All three buscoplotpy entry points call plt.savefig() with a hard-coded
    .png path and then plt.close() before returning, so the figure is gone by
    the time control comes back here. Wrapping plt.savefig for the duration of
    the call is the only way to catch the figure while it is still current,
    and it works identically for the barplot, karyoplot and synteny functions
    without depending on their internals beyond "they call plt.savefig".
    """
    if not enabled:
        yield
        return
    original = plt.savefig

    def patched(fname, *args, **kwargs):
        result = original(fname, *args, **kwargs)
        try:
            if isinstance(fname, (str, os.PathLike)):
                path = str(fname)
                if path.lower().endswith(".png"):
                    original(path[:-4] + ".svg", *args, **kwargs)
        except Exception as exc:                      # noqa: BLE001
            print(f"    ! SVG companion failed for {fname}: {exc}", file=sys.stderr)
        return result

    plt.savefig = patched
    try:
        yield
    finally:
        plt.savefig = original


SPECIES_COLORS = {"LEPC": "#8f5317", "GRPC": "#2f6f4e", "STGR": "#3a5a8c",
                  "GALGAL": "#6b6b6b"}   # Gallus gallus reference, for contrast

# Link palette for synteny. Default link colour in buscoplotpy is a flat grey,
# which makes a 40-chromosome comparison unreadable; colouring each link by its
# source chromosome is what lets you actually see a translocation.
LINK_PALETTE = [
    "#4c72b0", "#dd8452", "#55a868", "#c44e52", "#8172b3", "#937860",
    "#da8bc3", "#8c8c8c", "#ccb974", "#64b5cd",
]
# GALGAL last so the chicken reference reads as the outgroup/yardstick in the
# combined barplot rather than as one of the study species.
SPECIES_ORDER = ["LEPC", "GRPC", "STGR", "GALGAL"]


# ---------------------------------------------------------------- discovery --
def parse_name(dirname):
    """
    <SPECIES>_<SAMPLE>_<HAP> -> (species, sample, hap), or None.

    HAP is hapN for the phased grouse assemblies, or 'ref' for an unphased
    reference genome included for contrast (e.g. GALGAL_GRCg7b_ref). Allowing
    'ref' is what lets the chicken reference flow through this script as just
    another assembly.
    """
    m = re.fullmatch(r"([A-Z]+)_([A-Za-z0-9.]+)_(hap\d+|ref)", dirname)
    return m.groups() if m else None


def find_assemblies(busco_dir):
    out = []
    for path in sorted(glob.glob(os.path.join(busco_dir, "*"))):
        if not os.path.isdir(path):
            continue
        parsed = parse_name(os.path.basename(path))
        if parsed is None:
            continue
        species, sample, hap = parsed
        summary = sorted(glob.glob(os.path.join(path, "short_summary.specific.*.json")))
        tables = sorted(glob.glob(os.path.join(path, "run_*", "full_table.tsv")))

        # Lineage actually used, taken from BUSCO's own naming rather than from
        # the JSON's lineage_dataset.name field, which has been observed to
        # report a parent lineage (e.g. "vertebrata_odb10" on a run whose
        # marker count, n=8338, is unambiguously aves_odb10). The run directory
        # and the summary filename are both written from the dataset BUSCO
        # really loaded, so they are the reliable source.
        lineage = None
        if tables:
            run_dir = os.path.basename(os.path.dirname(tables[0]))
            if run_dir.startswith("run_"):
                lineage = run_dir[4:]
        if lineage is None and summary:
            m = re.match(r"short_summary\.specific\.([^.]+)\.",
                         os.path.basename(summary[0]))
            if m:
                lineage = m.group(1)

        out.append({
            "name": os.path.basename(path),
            "species": species, "sample": sample, "hap": hap,
            "dir": path,
            "summary_json": summary[0] if summary else None,
            "full_table": tables[0] if tables else None,
            "lineage": lineage,
        })
    return out


# ------------------------------------------------------------------ loaders --
def read_fulltable(path, species, sample, hap):
    """
    load_busco_fulltable() hard-codes skiprows=2, which assumes exactly two
    comment lines before the '# Busco id' header. That holds for BUSCO 5.4.7
    but is brittle, so verify it and fall back to an equivalent parse keyed on
    the real header line rather than letting a shifted file fail obscurely.
    """
    try:
        with open(path) as fh:
            head = [next(fh, "") for _ in range(6)]
    except OSError as exc:
        print(f"    ! cannot read {path}: {exc}", file=sys.stderr)
        return None

    header_idx = next((i for i, l in enumerate(head) if l.startswith("# Busco id")), None)
    if header_idx is None:
        print(f"    ! no '# Busco id' header in {path}", file=sys.stderr)
        return None

    if header_idx == 2:
        ft = load_busco_fulltable(path, group=species, organism=sample,
                                  genome_version=hap)
    else:
        print(f"    note: header on line {header_idx + 1}, not 3 — using fallback parse")
        raw = pd.read_csv(path, skiprows=header_idx, sep="\t")
        ft = pd.DataFrame({
            "busco_id": raw["# Busco id"], "status": raw["Status"],
            "sequence": raw["Sequence"], "gene_start": raw["Gene Start"],
            "gene_end": raw["Gene End"], "strand": raw["Strand"],
            "score": raw["Score"], "length": raw["Length"],
            "group": species, "organism": sample, "genome_version": hap,
        })
        ft["sequence"] = ft["sequence"].map(
            lambda x: x.split(":")[0] if pd.notna(x) else None)
    return ft


def read_summary(path, species, sample, hap, lineage=None):
    """
    load_json_summary() indexes a fixed set of JSON keys and raises KeyError if
    the run used a different gene predictor. Fall back to the keys the barplot
    actually needs. Then add group/organism/version, which the barplot requires
    and the loader does not provide.
    """
    try:
        df = load_json_summary(path)
    except (KeyError, TypeError) as exc:
        print(f"    note: load_json_summary failed ({exc}) — using minimal parse")
        d = json.load(open(path))
        r, l = d["results"], d.get("lineage_dataset", {})
        df = pd.DataFrame({
            "dataset_name": l.get("name", "unknown"),
            "one_line_summary": r["one_line_summary"],
            "complete": r["Complete"], "single copy": r["Single copy"],
            "multi copy": r["Multi copy"], "fragmented": r["Fragmented"],
            "missing": r["Missing"], "n_markers": r["n_markers"],
        }, index=[0])

    # The barplot titles itself with dataset_name. Prefer the lineage taken
    # from BUSCO's run directory / summary filename, and say so loudly when the
    # JSON disagrees rather than silently mislabelling every figure.
    if lineage:
        reported = str(df["dataset_name"].iloc[0])
        if reported != lineage:
            print(f"    note: JSON reports lineage '{reported}' but the BUSCO run is "
                  f"'{lineage}' (n={df['n_markers'].iloc[0]}); using '{lineage}'")
        df["dataset_name"] = lineage

    df["group"] = species
    df["organism"] = sample
    df["version"] = hap          # barplot builds labels as organism + '_' + version
    return df


# Canonical chromosome order, applied identically to every assembly so that
# plots can be compared row by row: numbered autosomes ascending, then Z, W,
# MT, then any unplaced scaffolds.
CHR_TAIL_ORDER = {"chr_Z": 0, "chr_W": 1, "chr_MT": 2}


def natural_chr_key(name):
    """chr_1 < chr_2 < ... < chr_39 < chr_Z < chr_W < chr_MT < scaffold_*."""
    m = re.fullmatch(r"chr_(\d+)", name)
    if m:
        return (0, int(m.group(1)), "")
    if name.startswith("chr_"):
        return (1, CHR_TAIL_ORDER.get(name, 3), name)
    return (2, 0, name)


def build_karyotype(fai_path, organism_label, color, chrom_only=True, limit=None):
    """Karyotype frame in the shape buscoplotpy expects."""
    fai = pd.read_csv(fai_path, sep="\t", header=None,
                      names=["sequence", "end", "offset", "linebases", "linewidth"])
    if chrom_only:
        fai = fai[fai["sequence"].str.startswith("chr_")]
    if fai.empty:
        return None
    fai = fai.sort_values("sequence", key=lambda s: s.map(natural_chr_key))
    if limit:
        fai = fai.nlargest(limit, "end").sort_values(
            "sequence", key=lambda s: s.map(natural_chr_key))

    return pd.DataFrame({
        "chr": fai["sequence"].values,
        "sequence": fai["sequence"].values,
        "start": 1,
        "end": fai["end"].values,
        "organism": organism_label,
        "color": color,
    }).reset_index(drop=True)     # synteny reads karyotype['color'][0] by label


# -------------------------------------------------------------------- plots --
def plot_barplots(assemblies, summaries, outdir, svg=True):
    os.makedirs(outdir, exist_ok=True)
    made = []
    for species in SPECIES_ORDER:
        rows = [summaries[a["name"]] for a in assemblies
                if a["species"] == species and a["name"] in summaries]
        if not rows:
            continue
        df = pd.concat(rows, ignore_index=True)
        with also_save_svg(svg):
            organism_busco_barplot(df=df, group_name=species, out_path=outdir + os.sep,
                                   filename=f"busco_barplot_{species}", dpi=200)
        plt.close("all")
        # organism_busco_barplot appends "_completeness.png" to `filename`.
        made.append(f"busco_barplot_{species}_completeness.png ({len(df)} assemblies)")

    rows = [summaries[a["name"]] for a in assemblies if a["name"] in summaries]
    if rows:
        df = pd.concat(rows, ignore_index=True)
        # The barplot labels bars as organism + '_' + version, so fold the
        # species into 'organism' to keep the combined panel unambiguous.
        df = df.copy()
        df["organism"] = df["group"] + "_" + df["organism"]
        with also_save_svg(svg):
            organism_busco_barplot(df=df, group_name="all species", out_path=outdir + os.sep,
                                   filename="busco_barplot_ALL", dpi=200)
        plt.close("all")
        made.append(f"busco_barplot_ALL_completeness.png ({len(df)} assemblies)")
    return made


def plot_karyotypes(assemblies, fulltables, final_dir, outdir, chrs_limit, svg=True):
    os.makedirs(outdir, exist_ok=True)
    made, skipped = [], []
    for a in assemblies:
        ft = fulltables.get(a["name"])
        fai = os.path.join(final_dir, f"{a['name']}.pseudo_chr.fasta.fai")
        if ft is None or not os.path.exists(fai):
            skipped.append(f"{a['name']} (missing {'full_table' if ft is None else 'fai'})")
            continue
        kt = build_karyotype(fai, f'{a["species"]} {a["sample"]} {a["hap"]}',
                             SPECIES_COLORS.get(a["species"], "#666666"))
        if kt is None:
            skipped.append(f"{a['name']} (no chr_* sequences in fai)")
            continue
        out = os.path.join(outdir, f"{a['name']}.png")
        # karyoplot lowercases and may reindex the frame it is handed, so pass a copy.
        # karyoplot renders karyotype['organism'][0] + ' ' + title, so `title`
        # must complement the organism label rather than repeat it.
        with also_save_svg(svg):
            # selected_sequences is passed ALWAYS, and this is load-bearing.
            # When an assembly has more chromosomes than chrs_limit, karyoplot
            # otherwise re-sorts the karyotype by length descending AND keeps
            # only the chrs_limit chromosomes with the most BUSCO hits — so
            # assemblies with different chromosome counts came out in different
            # row orders, with some chromosomes silently dropped. Supplying
            # selected_sequences takes the branch that filters without
            # reordering, so the canonical order built above survives and every
            # assembly is plotted row-for-row comparable.
            karyoplot(karyotype=kt.copy(), fulltable=ft.copy(), output_file=out,
                      title=f"- BUSCO positions ({a['lineage'] or 'unknown lineage'})",
                      selected_sequences=list(kt["sequence"]),
                      chrs_limit=chrs_limit, dpi=200)
        plt.close("all")
        made.append(os.path.basename(out))
    return made, skipped


def plot_synteny(pairs, assemblies, fulltables, final_dir, outdir, orientation,
                 chrs_limit, svg=True):
    os.makedirs(outdir, exist_ok=True)
    by_name = {a["name"]: a for a in assemblies}
    fn = horizontal_synteny_plot if orientation == "horizontal" else vertical_synteny_plot
    made, skipped = [], []

    for n1, n2 in pairs:
        if n1 not in by_name or n2 not in by_name:
            skipped.append(f"{n1} vs {n2} (assembly not found)")
            continue
        ft1, ft2 = fulltables.get(n1), fulltables.get(n2)
        fai1 = os.path.join(final_dir, f"{n1}.pseudo_chr.fasta.fai")
        fai2 = os.path.join(final_dir, f"{n2}.pseudo_chr.fasta.fai")
        if ft1 is None or ft2 is None or not (os.path.exists(fai1) and os.path.exists(fai2)):
            skipped.append(f"{n1} vs {n2} (missing full_table or fai)")
            continue

        kt1 = build_karyotype(fai1, n1, SPECIES_COLORS.get(by_name[n1]["species"], "#8f5317"),
                              limit=chrs_limit)
        kt2 = build_karyotype(fai2, n2, SPECIES_COLORS.get(by_name[n2]["species"], "#3a5a8c"),
                              limit=chrs_limit)
        if kt1 is None or kt2 is None:
            skipped.append(f"{n1} vs {n2} (empty karyotype)")
            continue

        shared = len(set(ft1[ft1["status"] == "Complete"]["busco_id"]) &
                     set(ft2[ft2["status"] == "Complete"]["busco_id"]))
        if shared == 0:
            skipped.append(f"{n1} vs {n2} (no shared Complete BUSCOs)")
            continue

        # Links are keyed on the LEFT assembly's chromosome (sequence_x).
        link_colors = {c: LINK_PALETTE[i % len(LINK_PALETTE)]
                       for i, c in enumerate(kt1["chr"])}

        out = os.path.join(outdir, f"{n1}__vs__{n2}.png")
        # The plot renders karyotype_1['organism'][0] + ' - ' +
        # karyotype_2['organism'][0] + ' ' + title, so `title` must not repeat
        # the two names.
        with also_save_svg(svg):
            fn(ft_1=ft1.copy(), ft_2=ft2.copy(),
               karyotype_1=kt1.copy(), karyotype_2=kt2.copy(),
               title="- shared Complete BUSCOs", link_colors=link_colors,
               output_path=out, dpi=200)
        plt.close("all")
        made.append(f"{os.path.basename(out)} ({shared} shared Complete BUSCOs)")
    return made, skipped


def default_pairs(assemblies):
    """
    One representative per species (lowest sample id, hap1 — or the 'ref'
    assembly for a reference genome), compared pairwise. With the chicken
    reference present this automatically yields each grouse species against
    chicken as well as the three grouse-vs-grouse comparisons.
    Deterministic so reruns produce the same figures.
    """
    reps = {}
    for a in sorted(assemblies, key=lambda x: (x["species"], x["sample"], x["hap"])):
        if a["hap"] in ("hap1", "ref"):
            reps.setdefault(a["species"], a["name"])
    present = [s for s in SPECIES_ORDER if s in reps]
    return [(reps[present[i]], reps[present[j]])
            for i in range(len(present)) for j in range(i + 1, len(present))]


def hap_pairs(assemblies):
    """hap1 vs hap2 within each individual."""
    haps = {}
    for a in assemblies:
        haps.setdefault((a["species"], a["sample"]), {})[a["hap"]] = a["name"]
    return [(v["hap1"], v["hap2"]) for _, v in sorted(haps.items())
            if "hap1" in v and "hap2" in v]


# --------------------------------------------------------------------- main --
def main():
    p = argparse.ArgumentParser(description="BUSCO completeness, karyotype and synteny plots.")
    p.add_argument("--busco-dir", required=True, help="qc/busco — one subdir per assembly")
    p.add_argument("--final-dir", required=True, help="final/ — holds *.pseudo_chr.fasta.fai")
    p.add_argument("--out-dir", required=True)
    p.add_argument("--chrs-limit", type=int, default=30,
                   help="max chromosomes drawn per assembly (default 30)")
    p.add_argument("--synteny-orientation", choices=["horizontal", "vertical"],
                   default="horizontal")
    p.add_argument("--synteny-pairs", default="",
                   help="comma-separated A:B pairs of assembly names; "
                        "default is one representative per species, pairwise")
    p.add_argument("--synteny-haps", action="store_true",
                   help="also plot hap1 vs hap2 for every individual")
    p.add_argument("--skip", default="", help="comma-separated: barplot,karyotype,synteny")
    p.add_argument("--no-svg", action="store_true",
                   help="write PNG only; by default an SVG is written beside every PNG")
    args = p.parse_args()

    skip = {s.strip() for s in args.skip.split(",") if s.strip()}
    svg = not args.no_svg

    assemblies = find_assemblies(args.busco_dir)
    if not assemblies:
        sys.exit(f"ERROR: no <SPECIES>_<SAMPLE>_<HAP> directories under {args.busco_dir}")

    print(f">>> {len(assemblies)} assemblies found")
    for s in SPECIES_ORDER:
        n = sum(1 for a in assemblies if a["species"] == s)
        if n:
            print(f"      {s}: {n}")

    print(">>> Loading BUSCO output")
    summaries, fulltables = {}, {}
    for a in assemblies:
        if a["summary_json"]:
            summaries[a["name"]] = read_summary(a["summary_json"], a["species"],
                                                a["sample"], a["hap"], a["lineage"])
        else:
            print(f"    ! {a['name']}: no short_summary JSON")
        if a["full_table"]:
            ft = read_fulltable(a["full_table"], a["species"], a["sample"], a["hap"])
            if ft is not None:
                fulltables[a["name"]] = ft
        else:
            print(f"    ! {a['name']}: no full_table.tsv")
    print(f"    summaries: {len(summaries)}/{len(assemblies)}   "
          f"full tables: {len(fulltables)}/{len(assemblies)}")

    if "barplot" not in skip:
        print(">>> Completeness barplots")
        for line in plot_barplots(assemblies, summaries,
                                  os.path.join(args.out_dir, "completeness"), svg=svg):
            print("    " + line)

    if "karyotype" not in skip:
        print(">>> Karyotype plots")
        made, skipped = plot_karyotypes(assemblies, fulltables, args.final_dir,
                                        os.path.join(args.out_dir, "karyotype"),
                                        args.chrs_limit, svg=svg)
        print(f"    {len(made)} written")
        for s in skipped:
            print(f"    skipped: {s}")

    if "synteny" not in skip:
        print(">>> Synteny plots")
        if args.synteny_pairs:
            pairs = [tuple(x.split(":", 1)) for x in args.synteny_pairs.split(",") if ":" in x]
        else:
            pairs = default_pairs(assemblies)
        if args.synteny_haps:
            pairs = pairs + hap_pairs(assemblies)
        print(f"    {len(pairs)} pair(s)")
        made, skipped = plot_synteny(pairs, assemblies, fulltables, args.final_dir,
                                     os.path.join(args.out_dir, "synteny"),
                                     args.synteny_orientation, args.chrs_limit, svg=svg)
        for m in made:
            print("    " + m)
        for s in skipped:
            print(f"    skipped: {s}")

    print(f">>> Done. Output under {args.out_dir}")


if __name__ == "__main__":
    main()
