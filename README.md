# ProtVis

[![](https://www.repostatus.org/badges/latest/active.svg)](https://www.repostatus.org/#active)
[![](https://img.shields.io/badge/GitHub-ProtVis-blue.svg)](https://github.com/anhuikylin/ProtVis)
[![](https://img.shields.io/badge/R-Shiny-orange.svg)](https://github.com/anhuikylin/ProtVis)
[![](https://img.shields.io/badge/license-MIT-green.svg)](https://github.com/anhuikylin/ProtVis/blob/dev/LICENSE)

<img src="https://raw.githubusercontent.com/anhuikylin/ProtVis/dev/app/www/ProtVis_ico.png" alt="ProtVis Logo" align="right" width="300"/>

**ProtVis: interactive visualization and downstream interpretation for proteomics and metaproteomics data**

ProtVis is an R/Shiny platform for reproducible proteomics analysis, from imported software outputs or raw mzML data to preprocessing, differential analysis, enrichment, GSEA, pathway interpretation, PTM/PSM inspection, multi-omics, and metaproteomics.

For detailed tutorials, parameter explanations, examples, and complete workflows, see the **ProtVis Cookbook**:

**https://anhuikylin.github.io/ProtVis-cookbook/**

<br clear="right"/>

---

## Highlights

- **Multiple input sources**: Raw mzML + FASTA, MaxQuant, DIA-NN, Spectronaut, FragPipe, Proteome Discoverer, Skyline, Mascot, OpenMS, and custom matrices.
- **Database search**: integrated Sage workflow and FragPipe headless backend.
- **Pre-processing**: filtering, transformation, imputation, normalization, missing-value visualization, dimensionality reduction, and QC.
- **Differential protein analysis**: Recommended DEP plus limma, DEqMS, proDA, and MSstats workflows.
- **Functional analysis**: enrichment analysis, directional KEGG, GSEA, pathway visualization, co-enrichment, PPI, WGCNA, Venn, and expression profiling.
- **PTM / PSM inspection**: peptide/PSM browsing, PTM-aware fragment matching, annotated spectra, and built-in phosphorylation / maize K-acetylation benchmarks.
- **Multi-omics and metaproteomics**: taxonomy, function, taxon × function, peptide-centric analysis, and compatible multi-omics inputs.
- **Reproducible project state**: `ProtVis_dataset`, checkpoints, provenance, append-only analysis runs, and resumable workflows.
- **GUI and CLI**: Shiny interface plus `run_protvis_cli()` for scripted workflows.

---

## Installation

Install the current development version from GitHub:

```r
options(repos = c(CRAN = "https://cloud.r-project.org"))
install.packages("pak")

pak::pak(c(
  "anhuikylin/ProtVisDatabase",
  "anhuikylin/ProtVis@dev"
))
```

`ProtVisDatabase` supplies the locally installed examples, backgrounds,
templates, benchmark files, and Sage executables used by ProtVis. No bundled
resource is downloaded while ProtVis is running.

For a clean reinstall/update:

```r
source(
  "https://raw.githubusercontent.com/anhuikylin/ProtVis/dev/install_ProtVis.R"
)
```

---

## Run ProtVis

```r
library(ProtVis)
run_ProtVis()
```

---

## Typical workflow

```text
Project init
   ↓
Import software output
or Raw mzML + FASTA → Search
   ↓
Pre-processing
   ↓
ProtVis_dataset
   ↓
Differential analysis
   ↓
Enrichment / GSEA / Pathway / Network
   ↓
PTM / Multi-omics / Metaproteomics
   ↓
Export + provenance
```

---

## Analysis conventions

Sage results are filtered after scoring using the requested q-value cutoff at
spectrum, peptide, protein and LFQ levels. Raw outputs are retained. Protein
quantification sums unique-peptide LFQ intensities; shared peptides remain in
the LFQ table but are excluded from individual-protein sums. Sage 0.14.7 traces
LFQ features only from peptide identifications passing its internal 1% gate.
Sage telemetry is disabled by the integrated runner.

Positive log2FC consistently means Group1 is higher than Group2. Differential
analysis requires verified log2-scale input and rejects raw, ln/log10 and
standardized values. For an imported log2 matrix without recorded scale, the
programmatic differential step accepts `input_scale = "log2"`.

See [validation results](dev/validation/analysis_repairs_20261007.md) for the
executed checks and remaining raw-data validation limit.

---

## Supported input sources

| Source | Typical input |
|---|---|
| Raw | mzML + protein FASTA |
| MaxQuant | proteinGroups / quantitative output |
| DIA-NN | report tables |
| Spectronaut | protein-group reports |
| FragPipe | protein reports |
| Proteome Discoverer | protein/peptide exports |
| Skyline | report tables |
| Mascot | exported result tables |
| OpenMS | consensus / protein quantification tables |
| Custom | protein-by-sample quantitative matrix |

All supported inputs are harmonized into the ProtVis project/data model for downstream analysis.

---

## Documentation

Full documentation and worked examples are maintained in the **ProtVis Cookbook**:

**https://anhuikylin.github.io/ProtVis-cookbook/**

The cookbook contains installation notes, input preparation, raw-data search, preprocessing, DEP, enrichment, GSEA, pathway analysis, PTM/PSM workflows, multi-omics, metaproteomics, and troubleshooting.

---

## Command-line use

```r
ProtVis::run_protvis_cli(c(
  "--input", "proteins.tsv",
  "--source", "DIA-NN",
  "--output", "results"
))
```

Use:

```r
ProtVis::protvis_cli_help()
```

for the available CLI options.

---

## Help

For usage details, first check the **ProtVis Cookbook**.  
For bugs or feature requests, please open an issue in this repository.

---

## Citation

If you use ProtVis in your research, please cite the corresponding ProtVis publication when available.

---

## Author

**Fei Liang**  
State Key Laboratory of Crop Stress Adaptation and Improvement  
Henan Joint International Laboratory for Crop Multi-Omics Research  
School of Life Sciences, Henan University  
Kaifeng 475004, China

---

## License

ProtVis is distributed under the repository license.
