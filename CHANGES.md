## Unreleased

- New **AutoSpill Compensation** operation (`AutoSpillOp`): compute the
  spillover matrix from single-color controls without gating positive and
  negative populations (robust regression + iterative refinement; a Python
  port of the AutoSpill R package, Roca et al. 2021), then compensate.
- "Magnetic" gates: polygon, freehand and ellipse gates can be adjusted per
  sample (`adjust_by`, "Adjust Per" in the GUI) -- the gate shifts to follow
  the population in each tube / well, up to `max_adjust` of its size.
- New **Axis Transforms** operation (`ScaleSettingsOp`): set the logicle
  W / M / A and arcsinh cofactor per channel (like FlowJo's "customize
  axis"); all later plots and gates use them.
- New **arcsinh** scale (`arcsinh(x / cofactor)`), standard for mass and
  spectral cytometry.
- Fix: an explicit logicle `W` passed to the scale could be ignored (the
  cached value didn't depend on it).
- Histogram overlays: "normalize to mode" (each histogram's peak = 100%)
  and staggered (offset) overlays, in the library and the GUI.
- Back-gating: scatter and density plots have a `highlight` setting
  ("Back-gate" in the GUI) that draws a gate's events on top of all events,
  in a contrasting color.
- New **Statistical Test** operation (`StatisticalTestOp`): one-way ANOVA,
  Kruskal-Wallis, Welch's t-test or Mann-Whitney U between the groups of a
  condition, with an optional Tukey HSD post-hoc test (p-values adjusted
  for multiple comparisons). Compares a statistic's replicates (eg the
  per-well geometric mean between Dox concentrations) or, optionally, the
  individual events, and reports effect sizes. Results are new statistics;
  the plot shows group means +/- SEM with significance brackets. The
  **Table View** can now show statistics whose index levels aren't
  experimental conditions (eg the pairwise Tukey table).
- New **Boolean Gate** operation (`BooleanGateOp`): combine existing gates
  with `and` / `or` / `not` into a new gate, eg `CD4_pos and not Dead` or
  `Y2_high and Quad == 'Quad_2'` (like FlowJo's Boolean gates and FCS
  Express's combination gates). The GUI lists the gates you can use.
- New **Derived Parameter** operation (`DerivedParameterOp`): make a new
  channel from an arithmetic expression of other channels, eg
  `FITC_A / PE_A` or `log10(B1_A + 1)`, with `+ - * / **` and `log`,
  `log10`, `log2`, `exp`, `sqrt`, `abs`, `arcsinh`, `min`, `max` (like
  FlowJo's Derive Parameters and FCS Express's Parameter Math). Events
  where the result isn't finite (eg division by zero) are dropped, with a
  warning. The GUI lists the channels you can use.
- Fix `BleedthroughLinearOp.set_from_fcs` reading the FCS spillover matrix
  transposed (so it compensated with the wrong coefficients).
- Undo / redo (Edit menu, Cmd/Ctrl+Z and Shift+Cmd/Ctrl+Z) for changes to
  operation and view settings -- including gates drawn on the plot -- and
  for adding or deleting operations.
- Zoom and pan on plots: pinch / two-finger scroll on a trackpad, mouse
  wheel or Ctrl/Cmd+scroll to zoom. The zoom is kept while gates are
  adjusted (the plot is redrawn after every change); press 0 to reset.
- New **Compensation** operation (`CompensationOp`): apply the spillover
  matrix the cytometer stored in the FCS file (`$SPILLOVER` / `SPILL` /
  `$SPILL`), like FlowJo and FCS Express do. Works with renamed channels,
  can load a matrix from a CSV file (eg exported from FlowJo), and the GUI
  shows the matrix as an editable table plus a heat map of the spillover.
- New UMAP (Uniform Manifold Approximation and Projection) dimensionality
  reduction operation, in both the library (`UMAPOp`) and the GUI. Like
  tSNE, it fits to a sample of the events and then embeds all of them;
  exposes the number of components, neighbors, minimum distance, metric and
  an optional random seed. Requires the `umap-learn` package.
- Channel axis labels: give each channel a free-text label at import (eg
  "CD4 (FITC)"), or fill them from the marker names recorded in the FCS
  file. Labels are used as axis titles on every plot. (`ImportOp.channel_labels`)
- SVG / PDF export keeps text as editable text, and embeds large scatter /
  density layers as an image, so exported figures open quickly in
  Illustrator, Affinity Designer and Inkscape.
- New **ellipse gate** (`EllipseOp`): drag out an ellipse, then drag its
  handles to resize it, drag the inside to move it, or press **r** and drag
  a corner handle to rotate it. The ellipse is defined on the scaled axes,
  so it's still an ellipse on log / logicle plots.
- New **freehand gate** (`FreehandOp`): drag a lasso around a population.
  The path is simplified to at most 100 vertices and gates exactly like a
  polygon gate.
- Easier gate drawing, especially on a laptop trackpad (eg, macOS):
  - Polygon gates close with a double-click; Backspace/Delete removes the
    last vertex; Escape starts over.
  - Existing polygon, range and 2D range gates can be adjusted by dragging
    their vertices / edges / corners, instead of re-drawing them.
  - A stray click no longer erases a range or 2D range gate.
- Threshold and quadrant gates: drag the line(s) to move them; a stray
  click no longer moves the gate (double-click jumps it to a new spot); a
  gate at 0 is now drawn; clicking outside the plot no longer erases a 
  quadrant gate.
- Smoother mouse tracking while drawing gates (60 updates/second, was 10).
- Keyboard shortcuts: Save As is now Ctrl/Cmd+Shift+S (was Ctrl/Cmd+E) and
  Save Plot is Ctrl/Cmd+Shift+E (was Ctrl/Cmd+X, which is normally Cut).
- New application icon.
- The GUI's plot now receives keyboard events (they were dropped before).
- Fix gates being drawn in the wrong place when the plot window is on a
  non-Retina screen (eg, an external monitor on a Mac).

## Release 1.2 "mCherry"
* Huge improvements to documentation.  Head over to https://cytoflow.readthedocs.io for tutorials, HOWTOs, explainers, and reference documentation.
* Implemented an "experiment browser" to help keep track of channels, conditions, and statistics in the GUI.
* Add an FCS export "view". Save calibrated or compensated data back to FCS files for analysis with other tools.
* Better support for high-DPI displays.
* Major GUI code restructure to separate the model from the view and controller.  As a result, cytoflow only creates one icon on the Mac taskbar now!
* Fix a number of bugs, including one high-profile bug where the GUI would crash on saving.
* Update dependencies to current versions.
  

## Release 1.1 "tdTomato"
Maintenance release.

* Update dependencies to current versions.
* Move CI infrastructure to GitHub Actions.
* Significant re-work of the GUI unit tests

## Release 1.0
Cytoflow's 1.0 release! This marks the state at which Cytoflow is, for all intents and purposes, "feature-complete." Major additions in this release include:

* Update all dependent libraries to their latest versions. Notably, this updates matplotlib to 3.0.2 and PyQt to 5.9.2.
* Overhauled the experimental setup dialog, including the ability to import CSV files with filenames and experimental conditions.
* Packaging differences make program startup much, much faster. Also, a proper Windows installer!

## Release 0.5 "GFP"
* Significant re-work of summary statistics generation and plotting
* Tighten down TASBE calibrated cytometry
* MANY many bug fixes, both in modules and GUI

## Release 0.4.1 
* A GUI for the calibrated cytometry tools!
* Bleedthrough compensation (both linear and piecewise-linear)
* Bead calibration
* Autofluorescence correction
* Color mapping
* Keep the same scales on plots in the same x, y facets
* Significant under-the-hood re-wiring
* Other bug-fixes

## Release 0.3.2
* Fix some Logicle crashers
* Add error bars to barplot, stats views
* Use FastLogicle -- it's faster!

## Release 0.3.0
Public beta!

* The point-and-click GUI is not only functional but useful. Serious restructuring behind the scenes, including a multi-process model to keep the GUI responsive.
* More operations (GMMs) and views (2D histogram, 1 and 2D KDEs, violin plots.)
*Lots (lots) of packaging work. One-click bundles and distribution infrastructure.

## Release 0.2.0
A developer's beta!

* Add a public API to Experiment to add conditions and channels.
* Add scaling options to views, data-driven modules
* Add 1D and 2D gaussian mixture models
* Clean up import mess
* Linear bleedthrough module, quad gate

