# Usage Guide

## Project Configuration

### 1. Set up the environment

```bash
pixi install
```

`pixi.toml` sets the folder paths as environment variables (`$RAWDATA`, `$PROCDATA`, `$RESULTS`, `$SCRIPTS`,
`$CONFIG`, `$METADATA`), so scripts never use hard-coded paths.

### 2. Add raw data

1. Download each dataset into its own folder under `data/rawdata/`, for example `data/rawdata/openpedcan_v15/`.
2. Add an entry for it in [`data_sources.md`](data_sources.md): source, version, access date, citation.
3. Never edit files in `data/rawdata/`. Anything derived from them goes in `data/procdata/`.

CBTN patient-level data must stay local. The `data/` folder is in `.gitignore`, and only its README files are tracked.

### 3. Configuration

Dataset versions and analysis cutoffs (for example, the maximum distance between the two transcription start sites
of a bidirectional pair) go in `config/`, not inside scripts.

## Running Your Analysis

Scripts live in `workflow/scripts/`, numbered in the order they run. Each step will also get a pixi task, so the
whole workflow can be run with `pixi run <task>`. The steps are:

1. Build the bidirectional gene pair list from GENCODE v39
2. Curate clinical data (ICB_Common standard, ONCOTREE names)
3. Curate molecular data (WGS, RNA-seq, methylation, microRNA) and build the `MultiAssayExperiment`
4. Screen promoter alterations against the expression of both partner genes
5. Compare by tumor type, subtype and survival

Commands for each step will be added here as the scripts are committed.
