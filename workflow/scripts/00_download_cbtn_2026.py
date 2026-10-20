"""
00_download_cbtn_2026.py
Download the newest CBTN files (CAVATICA project cavatica/pbta-cbttc) that curation needs for our cohort.

Cohort: patients with miRNA-Seq in OpenPedCan v15 in the six tumor types (miRNA has no newer release).
Files: RSEM gene results, consensus somatic MAF, copy-number segments (CNVkit, ControlFREEC, GATK),
EPIC methylation IDATs, plus the non-tumor methylation samples (normal brain reference).

The manifest "id" column does not match file IDs in the project, so files are matched by their path
(folder + name) and size. Files are fetched through signed URLs and written in binary mode: the
sevenbridges download() helper changes line endings on Windows and corrupts binary files.

Needs a CAVATICA token in ~/.sevenbridges/credentials (profile "cavatica"); the token is never printed.

Usage: python 00_download_cbtn_2026.py <rawdata_dir>
  reads  <rawdata_dir>/openpedcan_v15/histologies.tsv
  writes <rawdata_dir>/cbtn_2026/manifests_<date>/, <rawdata_dir>/cbtn_2026/files/<folder>/<name>,
         <rawdata_dir>/cbtn_2026/download_log.tsv
Safe to re-run: files already present with the right size are skipped.
"""
import csv
import datetime
import os
import shutil
import sys
import time
import urllib.request
from concurrent.futures import ThreadPoolExecutor, as_completed

import pandas as pd
import sevenbridges as sbg

PROJECT = "cavatica/pbta-cbttc"
SIX_TYPES = ["Low-grade glioma", "Ependymoma", "Medulloblastoma", "Ganglioglioma",
             "High-grade glioma", "Adamantinomatous Craniopharyngioma"]
KEEP_SUFFIXES = (
    ".rsem.genes.results.gz",                      # RNA-seq gene counts and TPM
    ".consensus_somatic.norm.annot.public.maf",    # somatic SNVs/indels, consensus of callers
    ".call.seg",                                   # CNVkit segments
    ".controlfreec.CNVs.p.value.txt",              # ControlFREEC calls
    ".tumor.gatk_cnv.called.seg",                  # GATK CNV segments
    ".idat", ".idat.gz",                           # EPIC methylation arrays
)

api = sbg.Api(config=sbg.Config(profile="cavatica"),
              error_handlers=[sbg.http.error_handlers.rate_limit_sleeper,
                              sbg.http.error_handlers.maintenance_sleeper])


def fetch(file, dest, tries=3):
    """Download one CAVATICA file byte-exact through its signed URL; skip if already complete."""
    if os.path.exists(dest) and os.path.getsize(dest) == file.size:
        return "present"
    os.makedirs(os.path.dirname(dest), exist_ok=True)
    for attempt in range(1, tries + 1):
        try:
            url = file.download_info().url
            with urllib.request.urlopen(url, timeout=120) as r, open(dest + ".part", "wb") as out:
                shutil.copyfileobj(r, out, length=1 << 20)
            if os.path.getsize(dest + ".part") != file.size:
                raise IOError(f"size {os.path.getsize(dest + '.part')} != {file.size}")
            os.replace(dest + ".part", dest)
            return "downloaded"
        except Exception as e:  # retry network hiccups, then report
            if attempt == tries:
                return f"failed: {type(e).__name__}: {str(e)[:120]}"
            time.sleep(5 * attempt)


def resolve_folder(path, cache):
    """Return the CAVATICA folder object for a path like 'CBTN-X01/simple-variants'."""
    if path in cache:
        return cache[path]
    parent_path, _, name = path.rpartition("/")
    if parent_path:
        hits = list(api.files.query(parent=resolve_folder(parent_path, cache), names=[name]))
    else:
        hits = list(api.files.query(project=PROJECT, names=[name]))
    if not hits:
        raise FileNotFoundError(f"folder not found in {PROJECT}: {path}")
    cache[path] = hits[0]
    return hits[0]


def main(rawdata):
    out = os.path.join(rawdata, "cbtn_2026")
    stamp = datetime.date.today().isoformat()

    # 1. Current manifests from the project's "manifest" folder (byte-exact)
    man_dir = os.path.join(out, f"manifests_{stamp}")
    folders = {}
    for f in api.files.query(parent=resolve_folder("manifest", folders)).all():
        print("manifest", f.name, fetch(f, os.path.join(man_dir, f.name)))

    # 2. Cohort and file selection
    h = pd.read_csv(os.path.join(rawdata, "openpedcan_v15", "histologies.tsv"), sep="\t", low_memory=False)
    cohort = set(h[(h.experimental_strategy == "miRNA-Seq") & (h.cancer_group.isin(SIX_TYPES))]
                 .Kids_First_Participant_ID)
    rows = []
    for name in sorted(os.listdir(man_dir)):
        with open(os.path.join(man_dir, name), newline="", encoding="utf-8") as fh:
            for r in csv.DictReader(fh, delimiter="\t"):
                in_cohort = r.get("Kids First Participant ID") in cohort
                normal_meth = name.startswith("source-data-CBTN-Methylation") and r.get("sample_type") == "Non-Tumor"
                if (in_cohort or normal_meth) and r["name"].endswith(KEEP_SUFFIXES):
                    rows.append({"name": r["name"], "size": int(r["size"]),
                                 "participant": r.get("Kids First Participant ID", ""),
                                 "biospecimen": r.get("Kids First Biospecimen ID", ""),
                                 "sample_type": r.get("sample_type", ""), "manifest": name})
    want = pd.DataFrame(rows).drop_duplicates("name")
    print(f"cohort patients {len(cohort)} | files selected {len(want)} | {want['size'].sum() / 1e9:.1f} GB")

    # 3. Match each selected path to a project file (by folder + name, checked by size)
    jobs, log = [], []
    for folder, grp in want.groupby(want["name"].str.rpartition("/")[0]):
        parent = resolve_folder(folder, folders)
        names = grp["name"].str.rpartition("/")[2].tolist()
        found = {}
        for i in range(0, len(names), 25):  # more names per query makes the URL too long
            for f in api.files.query(parent=parent, names=names[i:i + 25]).all():
                found[f.name] = f
        for _, r in grp.iterrows():
            base = r["name"].rpartition("/")[2]
            f = found.get(base)
            if f is None or f.size != r["size"]:
                log.append({**r, "status": "not found in project" if f is None else f"size differs ({f.size})"})
            else:
                jobs.append((f, os.path.join(out, "files", r["name"]), r))

    # 4. Download in parallel
    done = 0
    with ThreadPoolExecutor(max_workers=8) as pool:
        futures = {pool.submit(fetch, f, dest): r for f, dest, r in jobs}
        for fut in as_completed(futures):
            log.append({**futures[fut], "status": fut.result()})
            done += 1
            if done % 100 == 0 or done == len(jobs):
                print(f"{done}/{len(jobs)} files", flush=True)

    log = pd.DataFrame(log)
    log.to_csv(os.path.join(out, "download_log.tsv"), sep="\t", index=False)
    print(log["status"].str.split(":").str[0].value_counts().to_string())


if __name__ == "__main__":
    main(sys.argv[1])
