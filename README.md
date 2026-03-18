# STMotif

STMotif is an R package for mining spatial-temporal motifs from numeric matrices using SAX encoding and block-based motif discovery.

## Project Structure

- `R/`: Core R functions and package implementation (main algorithms, utilities, and visualization).
- `data/`: Example datasets used by package examples and tests.
- `inst/`: Installed files and documentation resources (e.g., `doc/STMotif.html`).
- `man/`: R documentation files generated from roxygen comments.
- `vignettes/`: Usage vignettes and examples.

## Key Functions

- `STSADatasetAdjust(D, tb, sb)`: Adjust dataset size to fit block dimensions.
- `NormSAX(D, a)`: Normalize and SAX-encode the dataset.
- `SearchSTMotifs(D, DS, w, a, sb, tb, si, ka)`: Find candidate spatial-temporal motifs.
- `RankSTMotifs(stmotifs)`: Rank motifs by quality metrics.
- `display_motifsDataset(...)`, `display_motifsSTSeries(...)`: Visualization functions.

## Documentation

- API documentation is in `man/` and generated help in `inst/doc/STMotif.html`.
- Usage vignette: `vignettes/STMotif.Rmd`.

## Quick Start

1. Install package from source using `devtools::install()`.
2. Load package and run example workflow:
   ```r
   library(STMotif)
   D <- STMotif::example_dataset
   DS <- NormSAX(D, 5)
   stmotifs <- SearchSTMotifs(D, DS, w=4, a=5, sb=4, tb=10, si=2, ka=2)
   rstmotifs <- RankSTMotifs(stmotifs)
   display_motifsDataset(D, rstmotifs[1:4], alpha=5)
   ```

