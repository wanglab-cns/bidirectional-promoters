"""
Step 1 (Aim 1a-1c): genome-wide catalogue of head-to-head (bidirectional) gene pairs from a GENCODE GTF,
optional CpG-island annotation, and a one-way (unidirectional) control set.

A pair = gene A on the minus strand and gene B on the plus strand whose transcription start sites (TSSs)
point away from each other. Distance d = TSS_plus - TSS_minus, minimised over all TSS combinations of the
two genes (every annotated transcript start counts, or only Ensembl-canonical ones with --canonical-only).

  d >= 0      divergent: 5' ends face away with d bp between them
  d < 0       overlapping: the two 5' ends overlap by |d| bp

Tiers written to the output:
  core        0 <= d <= --core     (default 1,000 bp; the usual literature definition)
  extended    --core < d <= --max  (default 5,000 bp; catches pairs with a wider shared promoter)
  overlapping -(--overlap) <= d < 0 (default 1,000 bp)

Optional extras:
  --cpg FILE   UCSC cpgIslandExt.txt(.gz): adds a cpg_island column (island overlapping the pair region)
  --controls   also writes one-way control promoters: genes in no pair whose Ensembl-canonical TSS has no
               opposite-strand TSS within --max bp on either side; region = canonical TSS +/- --flank
  --watch LIST comma-separated gene names: lists every pair that involves one of them in the summary

CpG islands: UCSC hg38 cpgIslandExt, https://hgdownload.soe.ucsc.edu/goldenPath/hg38/database/cpgIslandExt.txt.gz

Standard library only (no pandas needed). Usage (from the repo folder):
  python workflow/scripts/01_bidirectional_pairs.py \
      --gtf data/rawdata/gencode_v39/gencode.v39.annotation.gtf.gz --out data/results/pairs_GENCODE_v39 \
      --cpg data/rawdata/ucsc_hg38_cpg/cpgIslandExt.txt.gz --canonical-only --label canonical
"""
import argparse
import bisect
import collections
import csv
import gzip
import os
import re

ATTR_RE = re.compile(r'(\S+) "([^"]*)"')
DEFAULT_TYPES = "protein_coding,lncRNA,miRNA"


def _open(path):
    return gzip.open(path, "rt") if path.endswith(".gz") else open(path)


def parse_gtf(path, gene_types, canonical_only):
    """Return {gene_id: dict(name, type, chrom, strand, tss=set(), canonical_tss)} from transcript lines."""
    genes = {}
    with _open(path) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 9 or f[2] != "transcript":
                continue
            attrs = collections.defaultdict(list)
            for k, v in ATTR_RE.findall(f[8]):
                attrs[k].append(v)
            gtype = attrs["gene_type"][0]
            if gtype not in gene_types:
                continue
            canonical = "Ensembl_canonical" in attrs["tag"]
            if canonical_only and not canonical:
                continue
            chrom, start, end, strand = f[0], int(f[3]), int(f[4]), f[6]
            tss = start if strand == "+" else end  # 1-based, GTF convention
            gid = attrs["gene_id"][0]
            g = genes.setdefault(gid, {"name": attrs["gene_name"][0], "type": gtype, "chrom": chrom,
                                       "strand": strand, "tss": set(), "canonical_tss": None})
            g["tss"].add(tss)
            if canonical:
                g["canonical_tss"] = tss
    return genes


def tss_index(genes):
    """chrom -> strand -> (sorted positions, parallel gene ids)."""
    idx = collections.defaultdict(lambda: {"+": [], "-": []})
    for gid, g in genes.items():
        for t in g["tss"]:
            idx[g["chrom"]][g["strand"]].append((t, gid))
    out = {}
    for chrom, s in idx.items():
        out[chrom] = {}
        for strand, lst in s.items():
            lst.sort()
            out[chrom][strand] = ([p for p, _ in lst], [g for _, g in lst])
    return out


def find_pairs(genes, idx, core, max_d, overlap):
    """For every minus-strand TSS, look for plus-strand TSSs in [tss - overlap, tss + max_d]."""
    # A pair's tier is the best one any TSS combination reaches: core (0..core) beats overlapping (<0)
    # beats extended (core..max); within a tier the smallest |d| is kept.
    rank = lambda d: 0 if 0 <= d <= core else (1 if d < 0 else 2)
    best = {}  # (minus_gid, plus_gid) -> (d, tss_minus, tss_plus)
    for chrom, s in idx.items():
        plus_pos, plus_gid = s["+"]
        for tm, gm in zip(*s["-"]):
            lo = bisect.bisect_left(plus_pos, tm - overlap)
            hi = bisect.bisect_right(plus_pos, tm + max_d)
            for tp, gp in zip(plus_pos[lo:hi], plus_gid[lo:hi]):
                d = tp - tm
                key = (gm, gp)
                if key not in best or (rank(d), abs(d)) < (rank(best[key][0]), abs(best[key][0])):
                    best[key] = (d, tm, tp)
    rows = []
    for (gm, gp), (d, tm, tp) in best.items():
        tier = ("core", "overlapping", "extended")[rank(d)]
        a, b = genes[gm], genes[gp]
        rows.append({
            "chrom": a["chrom"], "tier": tier, "distance_bp": d,
            "minus_gene": a["name"], "minus_gene_id": gm, "minus_type": a["type"], "minus_tss": tm,
            "plus_gene": b["name"], "plus_gene_id": gp, "plus_type": b["type"], "plus_tss": tp,
            "pair_types": "|".join(sorted([a["type"], b["type"]])),
        })
    rows.sort(key=lambda r: (r["chrom"], min(r["minus_tss"], r["plus_tss"])))
    return rows


def load_cpg(path):
    """UCSC cpgIslandExt: chrom -> sorted list of (start0, end, name)."""
    isl = collections.defaultdict(list)
    with _open(path) as fh:
        for line in fh:
            f = line.rstrip("\n").split("\t")
            isl[f[1]].append((int(f[2]), int(f[3]), f"{f[1]}:{int(f[2]) + 1}-{f[3]}"))
    for v in isl.values():
        v.sort()
    return isl


def cpg_hit(isl, chrom, s0, e):
    """Name of a CpG island overlapping [s0, e) (0-based half-open), or ''."""
    lst = isl.get(chrom, [])
    i = bisect.bisect_left(lst, (e,))
    for j in range(i - 1, max(-1, i - 4), -1):  # islands are short and non-overlapping; look back a few
        if lst[j][1] > s0:
            return lst[j][2]
    return ""


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--gtf", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--types", default=DEFAULT_TYPES, help="comma-separated GENCODE gene_type values")
    ap.add_argument("--core", type=int, default=1000)
    ap.add_argument("--max", dest="max_d", type=int, default=5000)
    ap.add_argument("--overlap", type=int, default=1000)
    ap.add_argument("--flank", type=int, default=500, help="bp added outside the TSSs in BED regions")
    ap.add_argument("--canonical-only", action="store_true", help="use only Ensembl_canonical transcript TSSs")
    ap.add_argument("--cpg", help="UCSC cpgIslandExt.txt(.gz) for CpG-island annotation")
    ap.add_argument("--controls", action="store_true", help="also write one-way control promoters")
    ap.add_argument("--label", default="", help="suffix for output file names")
    ap.add_argument("--watch", default="", help="comma-separated gene names to list pairs for in the summary")
    a = ap.parse_args()

    os.makedirs(a.out, exist_ok=True)
    sfx = ("_" + a.label) if a.label else ""
    genes = parse_gtf(a.gtf, set(a.types.split(",")), a.canonical_only)
    idx = tss_index(genes)
    rows = find_pairs(genes, idx, a.core, a.max_d, a.overlap)
    isl = load_cpg(a.cpg) if a.cpg else None

    for r in rows:
        left, right = min(r["minus_tss"], r["plus_tss"]), max(r["minus_tss"], r["plus_tss"])
        r["region_start0"], r["region_end"] = max(0, left - 1 - a.flank), right + a.flank
        if isl is not None:
            r["cpg_island"] = cpg_hit(isl, r["chrom"], r["region_start0"], r["region_end"])

    with open(os.path.join(a.out, f"bidirectional_pairs{sfx}.tsv"), "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0].keys()), delimiter="\t")
        w.writeheader()
        w.writerows(rows)
    with open(os.path.join(a.out, f"bidirectional_promoters{sfx}.bed"), "w") as fh:
        for r in rows:
            fh.write(f'{r["chrom"]}\t{r["region_start0"]}\t{r["region_end"]}\t'
                     f'{r["minus_gene"]}|{r["plus_gene"]}\t{r["distance_bp"]}\t.\t{r["tier"]}\n')

    # One-way controls
    controls = []
    if a.controls:
        paired = {r[k] for r in rows for k in ("minus_gene_id", "plus_gene_id")}
        for gid, g in genes.items():
            t = g["canonical_tss"]
            if gid in paired or t is None:
                continue
            opp_pos, _ = idx[g["chrom"]]["-" if g["strand"] == "+" else "+"]
            lo = bisect.bisect_left(opp_pos, t - a.max_d)
            hi = bisect.bisect_right(opp_pos, t + a.max_d)
            if hi > lo:
                continue
            s0, e = max(0, t - 1 - a.flank), t + a.flank
            c = {"chrom": g["chrom"], "gene": g["name"], "gene_id": gid, "type": g["type"],
                 "strand": g["strand"], "tss": t, "region_start0": s0, "region_end": e}
            if isl is not None:
                c["cpg_island"] = cpg_hit(isl, g["chrom"], s0, e)
            controls.append(c)
        controls.sort(key=lambda c: (c["chrom"], c["tss"]))
        with open(os.path.join(a.out, f"unidirectional_controls{sfx}.tsv"), "w", newline="") as fh:
            w = csv.DictWriter(fh, fieldnames=list(controls[0].keys()), delimiter="\t")
            w.writeheader()
            w.writerows(controls)
        with open(os.path.join(a.out, f"unidirectional_controls{sfx}.bed"), "w") as fh:
            for c in controls:
                fh.write(f'{c["chrom"]}\t{c["region_start0"]}\t{c["region_end"]}\t{c["gene"]}\t0\t{c["strand"]}\n')

    # Summary
    tiers = collections.Counter(r["tier"] for r in rows)
    combos = collections.Counter((r["tier"], r["pair_types"]) for r in rows)
    pc = {gid for gid, g in genes.items() if g["type"] == "protein_coding"}
    pc_core = {r[k] for r in rows if r["tier"] == "core" for k in ("minus_gene_id", "plus_gene_id")} & pc
    pc_pc_core = sum(1 for r in rows if r["tier"] == "core" and r["pair_types"] == "protein_coding|protein_coding")
    lines = [f"GTF: {os.path.basename(a.gtf)}; gene types: {a.types}; canonical_only={a.canonical_only}",
             f"genes considered: {len(genes)}",
             f"pairs by tier: {dict(sorted(tiers.items()))}",
             f"protein-coding/protein-coding core pairs: {pc_pc_core}",
             f"protein-coding genes in >=1 core pair: {len(pc_core)} of {len(pc)} "
             f"({100 * len(pc_core) / max(1, len(pc)):.1f}%)"]
    if isl is not None:
        for tier in ("core", "extended", "overlapping"):
            sub = [r for r in rows if r["tier"] == tier]
            if sub:
                n = sum(1 for r in sub if r["cpg_island"])
                lines.append(f"CpG island overlap, {tier} pairs: {n} of {len(sub)} ({100 * n / len(sub):.1f}%)")
        if controls:
            n = sum(1 for c in controls if c["cpg_island"])
            lines.append(f"CpG island overlap, one-way controls: {n} of {len(controls)} ({100 * n / len(controls):.1f}%)")
            pcc = [c for c in controls if c["type"] == "protein_coding"]
            n = sum(1 for c in pcc if c["cpg_island"])
            lines.append(f"CpG island overlap, protein-coding one-way controls: {n} of {len(pcc)} "
                         f"({100 * n / max(1, len(pcc)):.1f}%)")
    lines.append("pairs by tier and gene types:")
    lines += [f"  {t:12s} {c:28s} {n}" for (t, c), n in sorted(combos.items())]
    known = [("BRCA1", "NBR2"), ("MLH1", "EPM2AIP1"), ("TP53", "WRAP53"), ("PRKN", "PACRG"),
             ("PTEN", "KLLN"), ("SDHD", "TIMM8B"), ("DPH3", "OXNAD1")]
    found = {frozenset((r["minus_gene"], r["plus_gene"])): r for r in rows}
    lines.append("known shared-promoter pairs recovered (literature examples):")
    for g1, g2 in known:
        r = found.get(frozenset((g1, g2)))
        lines.append(f"  {g1}/{g2}: " + (f'{r["distance_bp"]} bp [{r["tier"]}]' if r else "NOT FOUND"))
    watch = {g for g in a.watch.split(",") if g}
    if watch:
        lines.append(f"pairs involving {' / '.join(sorted(watch))}:")
        lines += [f'  {r["minus_gene"]} (-, {r["minus_type"]}) <-> {r["plus_gene"]} (+, {r["plus_type"]}): '
                  f'{r["distance_bp"]} bp [{r["tier"]}]' + (f' CpG {r["cpg_island"]}' if isl is not None else "")
                  for r in rows if {r["minus_gene"], r["plus_gene"]} & watch]
    with open(os.path.join(a.out, f"summary{sfx}.txt"), "w") as fh:
        fh.write("\n".join(lines) + "\n")
    print("\n".join(lines))


if __name__ == "__main__":
    main()
