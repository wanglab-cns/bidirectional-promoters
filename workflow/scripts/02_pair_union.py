"""
Step 2: one table of bidirectional gene pairs that keeps every pair found on GENCODE v50 or on GENCODE v39
(the annotation the CBTN data uses), with flags that say which data types each pair can be studied with.

Pairs are built with the exact rules of 01_bidirectional_pairs.py (its functions are imported): Ensembl-canonical
TSSs, gene types protein_coding / lncRNA / miRNA, tiers core (0-1 kb), extended (1-5 kb) and overlapping (up to
1 kb overlap). The two lists are then joined into one row per pair:

  matching    by Ensembl gene ID without version; a gene whose ID changed but kept its name (same chromosome,
              new ID not used before) is matched by name; pairs left unmatched are then matched by the pair of
              gene names on the same chromosome. Matching is one-to-one.
  tier        v39 tier when the pair is on v39 (the CBTN annotation), otherwise the v50 tier; tier_v50 and
              tier_v39 keep both
  reason      why a pair is on one release only, first that applies: gene not in the other release, gene type
              changed (no longer protein_coding / lncRNA / miRNA), strand or chromosome changed, canonical TSS
              moved (distance now outside -1 kb to 5 kb). Pairs on both releases with a different tier get
              "tier changed: canonical TSS moved".

Flags (genes are looked up by their v39 gene ID, version included):
  both_genes_in_cbtn_rna  both genes are rows in the CBTN RSEM files (pairs for RNA analyses); pairs without it
                          are methylation / WGS only
  expressed_flag          both genes TPM >= 1 in >= --expr-frac of the RSEM samples
  n_epic_probes           EPIC probes (v39 probe annotation) inside the shared promoter region
  cpg_island              UCSC CpG island overlapping the shared promoter region
  mirna_on_htg_panel      a miRNA partner matches a mature miRNA on the HTG panel. Rule: same mir/let prefix and
                          same number + letters (MIR133A1 -> mir-133a), and the panel name carries no copy number
                          or the gene's copy number (MIR24-1 matches miR-24-3p, not miR-24-2-5p)
  possible_v39_match      v50-only pairs with a gene missing from v39: the nearest v39 gene (any type) on the same
                          chromosome and strand whose gene body overlaps the missing gene. "same promoter" when its
                          canonical TSS is within --match-window bp of the missing gene's TSS, else "overlapping
                          gene". rna_via_possible_match: RNA data exists for both genes once same-promoter matches
                          are used; v39_pair_at_match: the v39 pair (pair_id) formed by the matched gene(s), if any.

The shared promoter region is both canonical TSSs +/- --flank bp, from v39 when the pair is on v39, else from v50.

Needs pandas and numpy, and Rscript on PATH to read the miRNA panel names from the OpenPedCan .rds file.
Usage (from the repo folder):
  pixi run python workflow/scripts/02_pair_union.py \
      --gtf50 data/rawdata/gencode_v50/gencode.v50.annotation.gtf.gz \
      --gtf39 data/rawdata/gencode_v39/gencode.v39.annotation.gtf.gz \
      --cpg data/rawdata/ucsc_hg38_cpg/cpgIslandExt.txt.gz \
      --probes data/rawdata/openpedcan_v15/infinium.gencode.v39.probe.annotations.tsv.gz \
      --mirna data/rawdata/openpedcan_v15/mirna-expression-counts.rds \
      --rsem-dir data/rawdata/cbtn_2026/files --out data/results/pairs_union_v50_v39
"""
import argparse
import bisect
import collections
import glob
import gzip
import importlib.util
import os
import re
import subprocess

import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
_spec = importlib.util.spec_from_file_location("pairs01", os.path.join(HERE, "01_bidirectional_pairs.py"))
pairs01 = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(pairs01)

TYPES = set(pairs01.DEFAULT_TYPES.split(","))
ATTR_RE = pairs01.ATTR_RE
TIERS = ("core", "extended", "overlapping")


def nov(gid):
    """Gene ID without version, keeping the PAR_Y suffix: ENSG00000182378.15_PAR_Y -> ENSG00000182378_PAR_Y."""
    return re.sub(r"\.\d+", "", gid)


def gene_table(path):
    """All genes of a GTF: {id without version: dict(id, name, type, chrom, strand, start, end, ctss, ctx)}."""
    genes = {}
    with gzip.open(path, "rt") as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            if f[2] == "gene":
                a = dict(ATTR_RE.findall(f[8]))
                genes.setdefault(nov(a["gene_id"]), {}).update(
                    id=a["gene_id"], name=a["gene_name"], type=a["gene_type"], chrom=f[0], strand=f[6],
                    start=int(f[3]), end=int(f[4]))
            elif f[2] == "transcript" and 'tag "Ensembl_canonical"' in f[8]:
                a = dict(ATTR_RE.findall(f[8]))
                genes.setdefault(nov(a["gene_id"]), {}).update(
                    ctx=nov(a["transcript_id"]), ctss=int(f[3]) if f[6] == "+" else int(f[4]))
    return genes


def build_pairs(gtf, a):
    """Canonical-TSS pairs from one GTF with the rules of step 1."""
    genes = pairs01.parse_gtf(gtf, TYPES, True)
    rows = pairs01.find_pairs(genes, pairs01.tss_index(genes), a.core, a.max_d, a.overlap)
    p = pd.DataFrame(rows)
    for s in ("minus", "plus"):
        p[f"{s}_nov"] = p[f"{s}_gene_id"].map(nov)
    p["key"] = [frozenset(k) for k in zip(p.minus_nov, p.plus_nov)]
    p["nkey"] = [(c, frozenset(k)) for c, k in zip(p.chrom, zip(p.minus_gene, p.plus_gene))]
    return p


def name_index(genes):
    idx = collections.defaultdict(list)
    for k, x in genes.items():
        if not x["name"].startswith("ENSG"):
            idx[(x["name"], x["chrom"])].append(k)
    return idx


def map_gene(gid, src, dst, dst_names):
    """Same gene in the other release: same ID, else a unique same-name gene whose ID is new (not in src)."""
    if gid in dst:
        return gid
    x = src[gid]
    hits = dst_names.get((x["name"], x["chrom"]), [])
    return hits[0] if len(hits) == 1 and hits[0] not in src else None


def reason(gm, gp, src, dst, gmap, a):
    """Why a pair on one release (genes gm, gp: IDs without version) is not a pair on the other one."""
    found = []
    for side, gid in (("minus", gm), ("plus", gp)):
        t = gmap.get(gid)
        if t is None:
            found.append((0, "gene not in other release", f"{side} {src[gid]['name']} ({src[gid]['type']})"))
        elif dst[t]["type"] not in TYPES:
            found.append((1, "gene type changed", f"{side} {src[gid]['type']} -> {dst[t]['type']}"))
        elif "ctss" not in dst[t]:
            found.append((2, "no canonical transcript", side))
        elif dst[t]["chrom"] != src[gid]["chrom"] or dst[t]["strand"] != src[gid]["strand"]:
            found.append((3, "strand or chromosome changed", side))
    if found:
        found.sort()
        return found[0][1], "; ".join(x[2] for x in found)
    d = dst[gmap[gp]]["ctss"] - dst[gmap[gm]]["ctss"]
    moved = [s for s, g in (("minus", gm), ("plus", gp)) if src[g]["ctss"] != dst[gmap[g]]["ctss"]]
    return "canonical TSS moved", f"distance in other release {d} bp; TSS moved: {', '.join(moved)}"


def panel_names(rds):
    """Row names of the HTG miRNA count matrix, read with Rscript."""
    out = subprocess.run(["Rscript", "-e", "writeLines(rownames(readRDS(commandArgs(TRUE)[1])))", rds],
                         capture_output=True, text=True, check=True).stdout
    return [x.strip() for x in out.splitlines() if x.strip()]


def mirna_matcher(names):
    pan = []
    for x in names:
        m = re.match(r"^(miR|let)-(\d+[a-z]*)(?:-(\d+))?(?:-([35]p))?$", x)
        if m:
            pan.append((m.group(1).lower(), m.group(2), m.group(3), x))

    def hits(sym):
        m = re.match(r"^MIRLET(\d+)([A-Z]*)(\d*)(?:-(\d+))?$", sym)
        pre = "let"
        if not m:
            m = re.match(r"^MIR(\d+)([A-Z]*)(\d*)(?:-(\d+))?$", sym)
            pre = "mir"
        if not m:
            return []
        stem, copy = m.group(1) + m.group(2).lower(), m.group(3) or m.group(4) or None
        return [x for p_, s_, c_, x in pan if p_ == pre and s_ == stem and (c_ is None or c_ == copy)]
    return hits, len(pan)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--gtf50", required=True)
    ap.add_argument("--gtf39", required=True)
    ap.add_argument("--cpg", required=True, help="UCSC cpgIslandExt.txt(.gz), hg38")
    ap.add_argument("--probes", required=True, help="EPIC probe annotation (Chromosome, Location, Probe_ID)")
    ap.add_argument("--mirna", required=True, help="HTG miRNA count matrix (.rds), row names = miRNAs")
    ap.add_argument("--rsem-dir", required=True, help="folder searched for */gene-expressions/*.rsem.genes.results.gz")
    ap.add_argument("--out", required=True)
    ap.add_argument("--core", type=int, default=1000)
    ap.add_argument("--max", dest="max_d", type=int, default=5000)
    ap.add_argument("--overlap", type=int, default=1000)
    ap.add_argument("--flank", type=int, default=500)
    ap.add_argument("--expr-frac", type=float, default=0.10, help="share of samples with TPM >= 1 (expressed_flag)")
    ap.add_argument("--match-window", type=int, default=1000, help="max TSS offset for a same-promoter match")
    a = ap.parse_args()
    os.makedirs(a.out, exist_ok=True)

    # ---------- pairs and genes of both releases
    p50, p39 = build_pairs(a.gtf50, a), build_pairs(a.gtf39, a)
    g50, g39 = gene_table(a.gtf50), gene_table(a.gtf39)
    n50, n39 = name_index(g50), name_index(g39)
    to39 = {g: map_gene(g, g50, g39, n39) for g in set(p50.minus_nov) | set(p50.plus_nov)}
    to50 = {g: map_gene(g, g39, g50, n50) for g in set(p39.minus_nov) | set(p39.plus_nov)}

    # ---------- one-to-one matching: gene IDs first, then gene names among pairs still unmatched
    k39 = dict(zip(p39.key, p39.index))
    nk39 = dict(zip(p39.nkey, p39.index))
    match, used = {}, set()
    for i, gm, gp in zip(p50.index, p50.minus_nov, p50.plus_nov):
        k = frozenset((to39[gm], to39[gp])) if to39[gm] and to39[gp] else None
        if k in k39 and k39[k] not in used:
            match[i] = k39[k]
            used.add(k39[k])
    by_name = set()
    for i, nk in zip(p50.index, p50.nkey):
        if i not in match and nk in nk39 and nk39[nk] not in used:
            match[i] = nk39[nk]
            used.add(nk39[nk])
            by_name.add(i)

    # ---------- union rows
    rows = []

    def gene_cols(side, id50, id39):
        x50, x39 = g50.get(id50) if id50 else None, g39.get(id39) if id39 else None
        ref = x39 or x50
        return {f"{side}_gene": ref["name"], f"{side}_gene_id_v50": x50["id"] if x50 else "NA",
                f"{side}_gene_id_v39": x39["id"] if x39 else "NA", f"{side}_type_v50": x50["type"] if x50 else "NA",
                f"{side}_type_v39": x39["type"] if x39 else "NA", f"{side}_strand": ref["strand"],
                f"{side}_tss_v50": x50.get("ctss", "NA") if x50 else "NA",
                f"{side}_tss_v39": x39.get("ctss", "NA") if x39 else "NA"}

    def add(r50, r39, i50=None):
        r = r39 if r39 is not None else r50
        row = {"chrom": r["chrom"], "tier": r["tier"], "in_v50": r50 is not None, "in_v39": r39 is not None,
               "tier_v50": r50["tier"] if r50 is not None else "NA",
               "distance_v50": r50["distance_bp"] if r50 is not None else "NA",
               "tier_v39": r39["tier"] if r39 is not None else "NA",
               "distance_v39": r39["distance_bp"] if r39 is not None else "NA"}
        for side in ("minus", "plus"):
            id50 = r50[f"{side}_nov"] if r50 is not None else to50.get(r39[f"{side}_nov"])
            id39 = r39[f"{side}_nov"] if r39 is not None else to39.get(r50[f"{side}_nov"])
            row.update(gene_cols(side, id50, id39))
        row["pair_types"] = r["pair_types"]
        left, right = min(r["minus_tss"], r["plus_tss"]), max(r["minus_tss"], r["plus_tss"])
        row["region_start0"], row["region_end"] = max(0, left - 1 - a.flank), right + a.flank
        row["region_from"] = "v39" if r39 is not None else "v50"
        row["matched_by_name"] = i50 in by_name or bool(
            r50 is not None and r39 is not None and any(r50[f"{s}_nov"] != r39[f"{s}_nov"] for s in ("minus", "plus")))
        if r50 is not None and r39 is None:
            row["reason"], row["reason_detail"] = reason(r50["minus_nov"], r50["plus_nov"], g50, g39, to39, a)
        elif r39 is not None and r50 is None:
            row["reason"], row["reason_detail"] = reason(r39["minus_nov"], r39["plus_nov"], g39, g50, to50, a)
        elif r50["tier"] != r39["tier"]:
            row["reason"], row["reason_detail"] = "tier changed: canonical TSS moved", ""
        else:
            row["reason"], row["reason_detail"] = "", ""
        rows.append(row)

    for i, r in p50.iterrows():
        add(r, p39.loc[match[i]] if i in match else None, i)
    for j, r in p39.iterrows():
        if j not in used:
            add(None, r)
    u = pd.DataFrame(rows)
    u["_chr"] = u.chrom.map(lambda c: int(c[3:]) if c[3:].isdigit() else {"X": 23, "Y": 24, "M": 25}.get(c[3:], 26))
    u = u.sort_values(["_chr", "region_start0"]).drop(columns="_chr").reset_index(drop=True)
    u.insert(0, "pair_id", [f"BP{i:05d}" for i in range(1, len(u) + 1)])

    # ---------- CpG islands and EPIC probes in the shared promoter region
    isl = pairs01.load_cpg(a.cpg)
    u["cpg_island"] = [pairs01.cpg_hit(isl, c, s, e) for c, s, e in zip(u.chrom, u.region_start0, u.region_end)]
    probe = pd.read_csv(a.probes, sep="\t", usecols=["Chromosome", "Location", "Probe_ID"], dtype=str)
    probe = probe.drop_duplicates("Probe_ID")
    probe["Location"] = pd.to_numeric(probe.Location, errors="coerce")
    pos = {c: np.sort(g.Location.dropna().values) for c, g in probe.groupby("Chromosome")}
    u["n_epic_probes"] = [0 if c not in pos else
                          int(np.searchsorted(pos[c], e, "right") - np.searchsorted(pos[c], s + 1, "left"))
                          for c, s, e in zip(u.chrom, u.region_start0, u.region_end)]

    # ---------- CBTN RNA (RSEM gene results; gene_id looks like ENSG00000000003.15_TSPAN6)
    files = sorted(glob.glob(os.path.join(a.rsem_dir, "**", "gene-expressions", "*.rsem.genes.results.gz"),
                             recursive=True))
    ids, tpm = None, []
    for f in files:
        d = pd.read_csv(f, sep="\t", usecols=["gene_id", "TPM"])
        if ids is None:
            ids = d.gene_id.values
        elif not np.array_equal(ids, d.gene_id.values):
            raise SystemExit(f"gene rows differ from the first RSEM file: {f}")
        tpm.append(d.TPM.values)
    rsem_id = [re.sub(r"_(PAR_Y_)?.*$", lambda m: "_PAR_Y" if m.group(1) else "", x) for x in ids]
    frac = dict(zip(rsem_id, (np.column_stack(tpm) >= 1).mean(axis=1)))  # genes x samples -> share with TPM >= 1
    for s in ("minus", "plus"):
        u[f"{s}_frac_tpm_ge1"] = [round(frac[x], 3) if x in frac else "NA" for x in u[f"{s}_gene_id_v39"]]
    u["both_genes_in_cbtn_rna"] = u.minus_gene_id_v39.isin(frac) & u.plus_gene_id_v39.isin(frac)
    u["expressed_flag"] = [bool(f1 != "NA" and f2 != "NA" and f1 >= a.expr_frac and f2 >= a.expr_frac)
                           for f1, f2 in zip(u.minus_frac_tpm_ge1, u.plus_frac_tpm_ge1)]

    # ---------- miRNA partners on the HTG panel
    hits, n_panel = mirna_matcher(panel_names(a.mirna))
    names = []
    for _, r in u.iterrows():
        h = []
        for s in ("minus", "plus"):
            if "miRNA" in (r[f"{s}_type_v39"], r[f"{s}_type_v50"]):
                h += hits(r[f"{s}_gene"])
        names.append(";".join(h))
    u["htg_panel_names"] = names
    u["mirna_on_htg_panel"] = u.htg_panel_names != ""

    # ---------- possible v39 match for genes of v50-only pairs that are missing from v39
    by_cs = collections.defaultdict(list)
    for k, x in g39.items():
        if "ctss" in x:
            by_cs[(x["chrom"], x["strand"])].append(k)
    arr = {cs: (np.array([g39[k]["start"] for k in ks]), np.array([g39[k]["end"] for k in ks]), ks)
           for cs, ks in by_cs.items()}
    pair_of = {frozenset((nov(m), nov(p))): pid for pid, m, p, on39 in
               zip(u.pair_id, u.minus_gene_id_v39, u.plus_gene_id_v39, u.in_v39) if on39}
    pm, pmd, rna_pm, v39pair = [], [], [], []
    for _, r in u.iterrows():
        if r.in_v39 or "gene not in other release" not in r.reason:
            pm.append(""), pmd.append(""), rna_pm.append(False), v39pair.append("")
            continue
        level, detail, ids39, rna_ok = [], [], [], True
        for s in ("minus", "plus"):
            if r[f"{s}_gene_id_v39"] != "NA":
                ids39.append(nov(r[f"{s}_gene_id_v39"]))
                rna_ok = rna_ok and r[f"{s}_gene_id_v39"] in frac
                continue
            x = g50[nov(r[f"{s}_gene_id_v50"])]
            st, en, ks = arr.get((x["chrom"], x["strand"]), (np.array([]), np.array([]), []))
            cand = [ks[i] for i in np.where((st <= x["end"]) & (en >= x["start"]))[0]]
            if not cand:
                level.append(0), detail.append(f"{s}: none"), ids39.append(None)
                rna_ok = False
                continue
            best = min(cand, key=lambda k: abs(g39[k]["ctss"] - x["ctss"]))
            off = g39[best]["ctss"] - x["ctss"]
            same = abs(off) <= a.match_window
            level.append(2 if same else 1)
            detail.append(f"{s}: {g39[best]['id']} ({g39[best]['name']}, {g39[best]['type']}), TSS offset {off} bp")
            ids39.append(best if same else None)
            rna_ok = rna_ok and same and g39[best]["id"] in frac
        lv = min(level)
        pm.append({0: "", 1: "overlapping gene", 2: "same promoter"}[lv])
        pmd.append("; ".join(detail))
        rna_pm.append(bool(lv == 2 and rna_ok))
        v39pair.append(pair_of.get(frozenset(ids39), "") if None not in ids39 else "")
    u["possible_v39_match"], u["possible_v39_match_detail"] = pm, pmd
    u["rna_via_possible_match"], u["v39_pair_at_match"] = rna_pm, v39pair

    u.to_csv(os.path.join(a.out, "pairs_union_v50_v39.tsv"), sep="\t", index=False)

    # ---------- summary
    L = [f"GTFs: {os.path.basename(a.gtf50)}, {os.path.basename(a.gtf39)}; canonical TSS; gene types {pairs01.DEFAULT_TYPES}",
         "pairs per release (check against step 1): "
         + "; ".join(f"{v} " + ", ".join(f"{t} {int((p.tier == t).sum())}" for t in TIERS)
                     for v, p in (("v50", p50), ("v39", p39))),
         f"RSEM files: {len(files)}; expressed_flag = both genes TPM >= 1 in >= {a.expr_frac:.0%} of them",
         f"HTG panel: {n_panel} mature miRNA names; EPIC probes: {len(probe)}",
         f"union rows: {len(u)} (tier = v39 tier when on v39, else v50 tier)", ""]
    hdr = f"{'':40s}" + "".join(f"{t:>13s}" for t in TIERS) + f"{'all':>9s}"
    L.append(hdr)

    def line(label, mask):
        c = u[mask].tier.value_counts()
        L.append(f"{label:40s}" + "".join(f"{int(c.get(t, 0)):>13d}" for t in TIERS) + f"{int(mask.sum()):>9d}")

    t = pd.Series(True, index=u.index)
    line("total pairs", t)
    line("  on both releases", u.in_v50 & u.in_v39)
    line("  v50 only", u.in_v50 & ~u.in_v39)
    line("  v39 only", ~u.in_v50 & u.in_v39)
    line("both genes in CBTN RNA (RNA set)", u.both_genes_in_cbtn_rna)
    line("  of which expressed_flag", u.both_genes_in_cbtn_rna & u.expressed_flag)
    line("methylation / WGS only (no RNA)", ~u.both_genes_in_cbtn_rna)
    line("  v50-only pairs with RNA", u.in_v50 & ~u.in_v39 & u.both_genes_in_cbtn_rna)
    line(">=1 EPIC probe in promoter", u.n_epic_probes >= 1)
    line("CpG island in promoter", u.cpg_island != "")
    line("pairs with a miRNA", u.pair_types.str.contains("miRNA"))
    line("  miRNA on HTG panel", u.mirna_on_htg_panel)
    line("possible v39 match: same promoter", u.possible_v39_match == "same promoter")
    line("possible v39 match: overlapping gene", u.possible_v39_match == "overlapping gene")
    line("  RNA via same-promoter match", u.rna_via_possible_match)
    line("  matched gene(s) form a v39 pair", u.v39_pair_at_match != "")
    L += ["", "reasons, pairs on one release only (all tiers):"]
    for lab, m in (("v50 only", u.in_v50 & ~u.in_v39), ("v39 only", ~u.in_v50 & u.in_v39)):
        L.append(f"  {lab}: " + ", ".join(f"{k} {v}" for k, v in u[m].reason.value_counts().items()))
    L.append(f"  on both releases, tier changed: {int((u.reason == 'tier changed: canonical TSS moved').sum())}")
    L.append(f"  matched by gene name (Ensembl ID changed): {int(u.matched_by_name.sum())}")
    with open(os.path.join(a.out, "summary.txt"), "w") as fh:
        fh.write("\n".join(L) + "\n")
    print("\n".join(L))


if __name__ == "__main__":
    main()
