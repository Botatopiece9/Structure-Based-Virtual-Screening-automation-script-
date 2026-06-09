# vs_screen.sh — Automated Virtual Screening Pipeline

A general-purpose bash pipeline for Structure-Based Virtual Screening using AutoDock Vina. Automates receptor and ligand preparation, docking across a compound library, score extraction, and ranked results reporting — in a single command.

![Pipeline Overview](assets/pipeline_overview.png)

---

## Features

- Screens any ligand library (`.sdf`, `.mol2`, `.pdbqt`) against any receptor
- Accepts a pre-written Vina config file **or** raw grid box coordinates as flags
- Automatically prepares receptor and ligands to `.pdbqt` using MGLTools
- Skips already-prepared files to support interrupted runs
- Extracts best binding affinity from each docking log
- Generates a ranked TSV report and a human-readable summary
- Timestamped output directories — re-runs never overwrite previous results
- Validates all inputs and dependencies before starting

---

## Dependencies

| Tool | Purpose | Install |
|---|---|---|
| AutoDock Vina | Molecular docking engine | [vina.scripps.edu](https://vina.scripps.edu/downloads/) |
| MGLTools | Receptor/ligand preparation (`prepare_receptor4.py`, `prepare_ligand4.py`) | [ccsb.scripps.edu/mgltools](https://ccsb.scripps.edu/mgltools/) |
| Python 3 | Required by MGLTools scripts | `sudo apt install python3` |
| bsdmainutils | Provides `column` for report formatting | `sudo apt install bsdmainutils` |

The script automatically checks that all tools are available before starting and prints install hints for anything missing.

---

## Installation

```bash
git clone https://github.com/<your-username>/vs_screen.git
cd vs_screen
chmod +x vs_screen.sh
```

---

## Usage

```
./vs_screen.sh [OPTIONS]

Options:
  -r, --receptor        Receptor .pdb or .pdbqt file        (required)
  -l, --library         Directory containing ligand files    (required)
  -c, --config          Vina grid box config file            (optional*)
  --center-x            Grid box center X coordinate         (optional*)
  --center-y            Grid box center Y coordinate         (optional*)
  --center-z            Grid box center Z coordinate         (optional*)
  --size                Grid box size in Angstroms           (default: 20)
  -o, --output          Output directory                     (default: results_<timestamp>)
  -n, --top             Number of top hits to report         (default: 10)
  -e, --exhaustiveness  Vina exhaustiveness                  (default: 8)
  -h, --help            Show this help message

* Provide either -c OR --center-x/y/z, not both.
```

---

## Examples

**Using a pre-written Vina config file:**
```bash
./vs_screen.sh \
    -r  proteins/mpro.pdb \
    -l  ligands/zinc_fragments/ \
    -c  configs/mpro_grid.conf \
    -n  20
```

**Using grid box coordinates directly:**
```bash
./vs_screen.sh \
    -r  proteins/mpro.pdb \
    -l  ligands/zinc_fragments/ \
    --center-x -26.3 \
    --center-y  11.4 \
    --center-z -19.2 \
    --size 22 \
    -e  12 \
    -n  20
```

**Using a pre-prepared receptor (skips MGLTools step):**
```bash
./vs_screen.sh \
    -r  proteins/mpro.pdbqt \
    -l  ligands/fda_approved/ \
    -c  configs/mpro_grid.conf
```

---

## Grid Box Config File Format

If using `-c`, your config file should follow standard Vina format:

```
center_x = -26.346
center_y =  11.357
center_z = -19.264
size_x   =  22.5
size_y   =  22.5
size_z   =  22.5
exhaustiveness = 8
num_modes      = 9
energy_range   = 3
```

If `exhaustiveness` is missing from the file, the script appends it automatically using the `-e` default.

---

## Output Structure

Each run produces a timestamped directory containing:

```
results_20250605_143022/
├── report.txt                  ← human-readable summary report
├── results_sorted.tsv          ← all compounds ranked by affinity
├── top_hits.tsv                ← top N compounds only
├── raw_scores.tsv              ← raw parsed scores
├── vina_config.txt             ← config used for this run
├── prepared/
│   ├── receptor.pdbqt
│   ├── compound_01.pdbqt
│   └── ...
├── docking_results/
│   ├── compound_01_docked.pdbqt
│   ├── compound_01.log
│   └── ...
└── logs/
    ├── receptor_prep.log
    └── ligand_prep.log
```

---

## Example Report Output

```
═══════════════════════════════════════════════════
  Virtual Screening Report
  2025-06-05 14:30:22
═══════════════════════════════════════════════════

  Receptor:        proteins/mpro.pdb
  Ligand library:  ligands/zinc_fragments/
  Config:          results_20250605_143022/vina_config.txt
  Exhaustiveness:  8

───────────────────────────────────────────────────
  Summary
───────────────────────────────────────────────────

  Compounds screened:   50
  Best compound:        ZINC000003986735
  Best affinity:        -9.3 kcal/mol
  Mean affinity:        -6.60 kcal/mol
  Top N reported:       10

───────────────────────────────────────────────────
  Top 10 Hits
───────────────────────────────────────────────────

Compound              Affinity_kcal_mol  RMSD_lb  RMSD_ub
ZINC000003986735      -9.3               0.0      0.0
ZINC000001481109      -8.7               0.0      0.0
ZINC000002033990      -8.1               0.0      0.0
...

═══════════════════════════════════════════════════
```

---

## Ligand Library Format

Place all ligand files in a single directory. Supported formats:

- `.sdf` — standard from ZINC15, ChEMBL, PubChem
- `.mol2` — standard from ZINC15, MOE
- `.pdbqt` — pre-prepared; copied directly, no MGLTools step

Mixed formats in the same directory are supported.

---

## Notes

- **Binding affinity** is reported in kcal/mol. More negative = stronger predicted binding.
- **RMSD lb/ub** are lower and upper bound RMSDs of each pose relative to the best mode.
- If a ligand fails preparation, it is skipped and logged — the run continues.
- If all docking runs fail, the script exits with an error and points to the log files.
- MGLTools is expected at `/opt/mgltools/bin` by default. If installed elsewhere, add it to your `PATH` before running.

---

## Background

This pipeline was built as part of a computational biology portfolio project applying Structure-Based Drug Design (SBDD) principles to automated virtual screening. The underlying docking engine is [AutoDock Vina](https://vina.scripps.edu/), one of the most widely used open-source molecular docking programs in drug discovery research.

---

## License

MIT License — free to use, modify, and distribute.
