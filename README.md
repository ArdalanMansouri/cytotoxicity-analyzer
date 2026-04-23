# cytotox — Cytotoxicity Assay Analyzer

## Background

Cytotoxicity assays based on fluorescent markers such as Caspase 3/7 (an apoptosis indicator) and Propidium Iodide (PI, a necrosis/membrane integrity marker) generate large volumes of single-cell imaging data. Despite the widespread adoption of high-content imaging platforms, the bioinformatics tooling available to researchers for analyzing this type of data remains fragmented. Most existing workflows depend on manual gating in closed, GUI-driven software (e.g., proprietary image analysis platforms), spreadsheet-based calculations, or generalist data science libraries that require researchers to build analysis pipelines from scratch. There is no dedicated, open-source Python library that covers the full analytical workflow — from raw single-cell fluorescence data to classified cytotoxicity outcomes and publication-ready visualizations — in a transparent and reproducible way.

## Solution

`cytotox` is an open-source Python package that provides a complete, code-first pipeline for cytotoxicity assay analysis. It enables researchers to go from raw single-cell fluorescence intensity data (Caspase 3/7 and PI channels) to classified cell populations, replicate-level summary statistics, and interactive visualizations — all within a reproducible Jupyter notebook environment. Threshold-based and slope-gated classification are both supported, giving users the flexibility to match their experimental gating strategy. By building on standard scientific Python libraries (pandas, plotly, scipy, statsmodels, FlowCytometryTools), the package integrates naturally into existing data science workflows while providing biology-specific abstractions that eliminate repetitive boilerplate.

## Key Features

- **Dual-marker cytotoxicity classification** — Classifies individual cells into Viable, Apoptosis, Late Apoptosis, and Necrosis groups based on Caspase 3/7 and PI fluorescence intensity thresholds.
- **Slope-based gating** — Supports diagonal (line-based) gates defined by two coordinate points, enabling more precise separation of apoptotic and necrotic populations where rectangular gating is insufficient.
- **SD-threshold categorization** — Identifies inhibitors and inducers of a biological process (e.g., extracellular vesicle uptake) by computing mean/median ± N standard deviations from an untreated control group, with configurable N and aggregation function.
- **Per-well and per-sample summary tables** — Automatically aggregates single-cell data to replicate-level counts and percentages of live, Caspase-positive, and PI-positive cells via `compute_cytotox_table`.
- **Interactive scatter plots with faceting** — Generates Plotly-based faceted scatter plots for Caspase 3/7 vs. PI intensity, with automatic threshold lines, annotations, and support for both individual-well and merged-replicate views.
- **Gated scatter plots with color-coded populations** — Visualizes classified cell populations in color-coded faceted scatter plots, including customizable gating lines and a separator line between experimental groups.
- **Flexible replicate handling** — Supports merging or separating technical replicates across wells, with automatic facet layout calculation based on the number of samples.
- **Export support** — Figures can be exported directly to HTML (interactive) or static image formats from within the plotting functions.
- **Reproducible notebook workflows** — Example Jupyter notebooks demonstrate the full pipeline from raw `.txt` data files through cleaning, classification, and visualization.
