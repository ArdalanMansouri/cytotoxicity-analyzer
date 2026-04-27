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

## Interactive Analysis Notebooks
Due to limitations of GitHub rendering for Jupyter notebooks, the example analysis notebooks are provided as interactive `.html` files, accessable at [Cytotoxicity Plotting](https://ardalanmansouri.github.io/cytotoxicity-analyzer/) 

The original `.ipynb` files remain available in the [`notebooks/`](notebooks/) folder for those who wish to run or modify the analysis locally.

## Installation

To install the required dependencies, ensure you have Python 3.10 or later installed. It is recommended to use the provided `environment.yml` to create a Conda environment with all dependencies:

```bash
conda env create -f environment.yml
conda activate cytotox
```

Alternatively, install the package and its dependencies directly with pip:

```bash
pip install .
```

## Usage

1. Clone the repository:
   ```bash
   git clone <repository-url>
   ```
2. Navigate to the project directory:
   ```bash
   cd cytotoxicity-analyzer
   ```
3. Activate the Conda environment:
   ```bash
   conda activate cytotox
   ```
4. Open the Jupyter notebooks in the `notebooks/` directory to explore the data analysis workflows:
   - `cytotox_data_cleaning.ipynb` — Load and clean raw single-cell `.txt` files from high-content imaging.
   - `cytotox_plotting.ipynb` — Generate interactive scatter plots with threshold annotations.
   - `cytotoxicity_gated_plots.ipynb` — Classify cell populations and produce color-coded gated scatter plots.

## Contributing

Contributions are welcome! If you have suggestions for improvements or new features, feel free to open an issue or submit a pull request.

## License

This project is licensed under the MIT License. See the `LICENSE` file for details.

## Acknowledgments

This project is funded by Marie Curie Actions under the European Union’s Horizon 2020 research and innovation program for project proEVLifeCycle, grant No 860303.

Special thanks to the open-source community for providing the foundational libraries (Pandas, NumPy, Plotly, SciPy, Statsmodels, FlowCytometryTools) that make this project possible.
