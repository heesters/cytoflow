# Cytoflow-GCB

**A point-and-click flow cytometry analysis tool for germinal centre biology.**

Cytoflow-GCB is a fork of [Cytoflow](https://github.com/cytoflow/cytoflow) with additional operations, views and quality-of-life improvements developed in the [Heesters Lab](https://heesters.gitlab.io/) (Lab for Germinal Center Biology, Utrecht University).  It runs as a standalone desktop application — no Python installation required.

---

## Download

Pre-built packages for **macOS (Apple Silicon)** and **Windows** are on the [Releases](../../releases) page.  Download the latest release and follow the instructions for your platform.

### macOS (Intel and Apple Silicon)

1. Download `Cytoflow-GCB-macos_x86-64.zip` and unzip it.
2. Move **Cytoflow-GCB.app** to your Applications folder.
3. First launch: right-click (or Control-click) the app → **Open** → **Open** again in the dialog that appears.
4. If macOS still says the app is damaged, run this once in Terminal:
   ```
   xattr -dr com.apple.quarantine /Applications/Cytoflow-GCB.app
   ```
   Then try launching normally.

> **Why?** The app is not code-signed by Apple.  macOS quarantines anything downloaded from the internet; the steps above lift that quarantine.  The package is built for Intel (x86_64) and runs on Apple Silicon Macs via Rosetta 2, which macOS installs automatically on first use.

### Windows

1. Download `Cytoflow-GCB-win_x86-64.exe` and run it.
2. If Windows SmartScreen appears, click **More info** → **Run anyway**.

---

## What is Cytoflow?

Cytoflow treats a cytometry experiment as a **workflow**: a sequence of operations (gating, compensation, transforms, statistics) applied to imported FCS data, each building on the last.  Operations and views are configured in the left-hand panel; a plot updates live on the right.  Workflows are saved to `.flow` files and can be re-opened, shared and re-run on new data.

The design philosophy is [described in detail](https://github.com/cytoflow/cytoflow#whats-different-about-cytoflow) in the upstream repository.

---

## Additions in this fork

### Compensation and transforms

| Feature | Description |
|---|---|
| **AutoSpill compensation** | Compute the spillover matrix from single-colour controls without manual gating — robust regression + iterative refinement ([Roca et al. 2021](https://doi.org/10.1038/s41467-021-23126-8)) |
| **FCS spillover** | Apply the matrix stored in the FCS file (`$SPILLOVER`/`SPILL`) with one click, as FlowJo does |
| **Axis transforms** | Set per-channel logicle W / M / A and arcsinh cofactor — equivalent to FlowJo's "Customize Axis" |
| **arcsinh scale** | Standard for mass cytometry and spectral cytometry: `arcsinh(x / cofactor)` |
| **Channel scale defaults** | Persist a preferred scale (logicle, arcsinh, …) per channel across sessions |

### Gating

| Feature | Description |
|---|---|
| **Ellipse gate** | Draw an ellipse, then drag handles to resize, move, or rotate; remains a true ellipse on log/logicle axes |
| **Freehand gate** | Lasso a population; simplified to ≤ 100 vertices, gates like a polygon |
| **Magnetic gates** | Polygon, freehand and ellipse gates shift per tube/well to follow the population (`Adjust Per`) |
| **Boolean gate** | Combine existing gates with `and` / `or` / `not`, e.g. `CD4_pos and not Dead` |
| **GatingML / FlowJo import-export** | Exchange gates with FlowJo 10 workspaces and Gating-ML 2.0 files |
| **Gate Template Library** | Save any gate to a `.gatetemplate` file (Ctrl+T) and reload it into a new experiment (Ctrl+Shift+T) |
| **Improved gate drawing** | Polygon closes on double-click; vertices draggable; threshold/quadrant lines draggable; stray clicks no longer erase gates |

### Analysis

| Feature | Description |
|---|---|
| **Population Tree** | Nested populations with event counts, % of parent and % of total — as a statistic table and a tree plot; "Build from workflow" infers hierarchy from gate subsets or FlowJo import |
| **Proliferation** | Fit dye-dilution generations (CFSE, CellTrace Violet) — reports % divided, division, proliferation, expansion and replication indices |
| **Cell Cycle** | Fit G0/G1, S and G2/M to a DNA-content histogram (Dean–Jett–Fox model), with sub-G1 and above-G2/M fractions |
| **Statistical Test** | One-way ANOVA, Kruskal–Wallis, Welch's t-test or Mann–Whitney U; optional Tukey HSD post-hoc; compares per-replicate statistics or individual events; reports effect sizes and adjusted p-values |
| **Cross-Entropy Test** | Statistical test of whether groups differ in a t-SNE embedding ([Roca et al. 2023](https://doi.org/10.1016/j.crmeth.2023.100641)) |
| **PhenoGraph clustering** | k-nearest-neighbour graph with Jaccard weights + Leiden/Louvain community detection on any set of channels |
| **UMAP** | Dimensionality reduction; fits on a sample, embeds all events |
| **Derived parameter** | New channel from an arithmetic expression of other channels, e.g. `FITC_A / PE_A` |
| **Sample Labels** | Assign short, human-readable names to any condition; optionally stores them as a new condition |

### Visualisation

| Feature | Description |
|---|---|
| **Histogram overlay** | Overlay histograms per condition; normalise to mode (peak = 100%) or display staggered |
| **Back-gating** | Scatter and density plots highlight one gate's events on top of all events in a contrasting colour |
| **Zoom and pan** | Pinch / two-finger scroll, mouse wheel or Ctrl+scroll; zoom is kept during gate adjustment; press **0** to reset |
| **Undo / redo** | Edit → Undo/Redo (Cmd/Ctrl+Z) for all operation and view changes, including drawn gates |
| **Export report** | File → Export Report: selected workflow plots to PDF (one vector plot per page, text editable) or PowerPoint |
| **SVG / PDF export** | Scatter/density layers embedded as images; text remains editable in Illustrator, Affinity Designer and Inkscape |

### macOS improvements

- Dark-mode support (follows system appearance)
- Unified title+toolbar bar
- HiDPI / Retina opt-in
- Window position and size saved and restored across launches

---

## Running from source

Requires [Miniconda](https://docs.conda.io/en/latest/miniconda.html) (or Mamba).

```bash
git clone --recurse-submodules https://github.com/heesters/cytoflow.git
cd cytoflow
conda env create -f environment.yml
conda activate cytoflow
SETUPTOOLS_SCM_PRETEND_VERSION=1.3.2 pip install -e .
python cytoflowgui/run.py
```

On Windows, replace the `pip install` line with:
```
set SETUPTOOLS_SCM_PRETEND_VERSION=1.3.2 && pip install -e .
```

The `SETUPTOOLS_SCM_PRETEND_VERSION` variable is needed because the fork has no matching version tags; without it `setuptools_scm` may fail to determine a version.

To build the standalone `.app` / `.exe` installer, see the GitHub Actions workflow at `.github/workflows/build-apps.yml`.

---

## Based on

This project is built on [Cytoflow](https://github.com/cytoflow/cytoflow) by Brian Teague, distributed under the [GNU General Public License v2](LICENSE).  All additions in this fork are released under the same licence.
