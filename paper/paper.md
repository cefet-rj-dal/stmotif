---
title: "STMotif: An R Package for Spatial-Time Motif Discovery"
authors:
  - name: Heraldo Borges
    affiliation: 1
  - name: Antonio Castro
    affiliation: 1
  - name: Rafaelli Coutinho
    affiliation: 1
  - name: Eduardo Ogasawara
    affiliation: 1
affiliations:
  - name: Federal Center for Technological Education of Rio de Janeiro (CEFET/RJ), Brazil
    index: 1
bibliography: paper.bib
repository: https://github.com/cefet-rj-dal/STMotif
documentation: https://cran.r-project.org/package=STMotif
tags: [spatial-time series, motif discovery, time series, SAX, R]
date: 17 March 2026
---

# Summary

**STMotif** is an R package for discovering, ranking, and visualizing motifs in **spatial-time series**. Unlike traditional motif discovery tools that analyze each time series independently, STMotif targets patterns that are constrained simultaneously in space and time and may emerge only when neighboring series are analyzed jointly. The package implements the **Combined Series Approach (CSA)** introduced by @borges_spatial-time_2020 and provides end-to-end support for normalization, SAX encoding, motif search, ranking, and visual inspection of the resulting motifs.

# Statement of need

Many real-world phenomena can be represented as time series, and some of them are naturally associated with spatial location, producing spatial-time series datasets [@shumway_time_2017; @atkinson_gis_2000]. Examples include sensor networks, environmental monitoring, ocean currents, seismic data, and other geo-referenced temporal observations [@fu_review_2011]. In these scenarios, relevant patterns may not be strongly expressed in a single series, but rather distributed across multiple neighboring series.

Motif discovery aims to identify previously unknown and recurrent subsequences in time series [@lin_finding_2002; @mueen_time_2014]. Existing research has produced important advances in motif discovery algorithms and representations [@torkamani_survey_2017; @yeh_time_2018], but these methods mainly focus on single-series analysis. This creates a gap for applications in which the pattern of interest is both temporally localized and spatially distributed.

STMotif addresses this gap by providing a software implementation tailored to **motifs constrained in space and time**. The package is intended for researchers and practitioners who need to identify recurrent spatiotemporal patterns, rank them according to informative criteria, and inspect them visually in the original data. The main audience includes users working with georeferenced temporal data such as sensor networks, environmental monitoring, and related spatial analytics applications.

# State of the field

Symbolic representations such as SAX are widely used to reduce the complexity of motif discovery over real-valued time series [@lin_experiencing_2007]. They enable efficient indexing and comparison of subsequences while preserving relevant shape information. Traditional motif discovery methods, however, generally assume that the input is a single series or a collection of series analyzed separately.

The CSA method proposed by @borges_spatial-time_2020 extends this perspective by partitioning a spatial-time dataset into spatial and temporal blocks, combining the observations within each block, and then searching for motifs that satisfy both occurrence and spatial-distribution constraints. STMotif operationalizes this method in software and, to our knowledge, is the first R package focused on motif discovery under joint spatial and temporal constraints. This is therefore not a repackaging of an existing workflow: the contribution is the software realization of a distinct spatial-time motif mining method together with its ranking and visualization pipeline.

# Software design

STMotif is organized as a pipeline that separates symbolic preprocessing, constrained motif search, ranking, and visualization. This separation matters because users may need to inspect encoded data, tune discovery thresholds, or reuse intermediate results while exploring spatial-time datasets.

The motif discovery process begins with a spatial-time dataset in which each time series is associated with a spatial position. Since direct motif discovery over real-valued subsequences is inefficient, the package first applies **z-score normalization** followed by **SAX indexing** [@keogh_exact_2005; @lin_experiencing_2007]. SAX transforms numeric subsequences into symbolic words drawn from an alphabet of size $a$, producing a representation more suitable for frequent-pattern discovery.

The indexed dataset is then partitioned into blocks according to two user-defined constraints: spatial block size ($sb$) and temporal block size ($tb$). Inside each block, the package combines the subsequences into a single time series and examines all subsequences of size $w$. A candidate subsequence is considered a spatial-time motif if it satisfies both a minimum number of occurrences ($\sigma$) and a minimum number of distinct spatial-time series with occurrences ($\kappa$) inside the block [@borges_spatial-time_2020].

After discovery, STMotif ranks motifs according to three complementary criteria: number of occurrences, spatial-temporal proximity of occurrences, and entropy. The ranking phase favors motifs that are frequent, spatially and temporally coherent, and information-rich. This design choice is important because motif mining can return many candidates, and the ranking stage helps users focus on the most informative patterns without losing access to the complete result set.

![Figure 1: Main functionalities of STMotif, including normalization and SAX indexing, motif search, ranking, and visualization.](figures/STMotif_package.png)

The package exposes this workflow through a compact set of user-facing functions:

- `NormSAX()` applies normalization and SAX encoding to the input dataset.
- `SearchSTMotifs()` discovers motifs from the original and encoded datasets using the parameters $w$, $a$, $sb$, $tb$, $\sigma$, and $\kappa$.
- `RankSTMotifs()` ranks the discovered motifs according to occurrence count, proximity, and entropy.
- `CSAMiningProcess()` executes the complete workflow in a single call.
- `display_motifsSTSeries()` visualizes motifs over selected spatial-time series.
- `display_motifsDataset()` displays motif occurrences over the whole dataset.

This design separates preprocessing, mining, ranking, and visualization, while also allowing the full pipeline to be executed compactly when desired. In practice, this supports both exploratory use by domain researchers and scripted, reproducible use in larger analyses.

# Research impact statement

STMotif operationalizes the CSA method published by @borges_spatial-time_2020 as reusable research software, lowering the barrier for motif discovery in spatial-time datasets. Beyond motif identification, the package emphasizes interpretability through dedicated visualization functions, which help users inspect the spatial-temporal distribution of the discovered patterns and connect symbolic motifs back to the original data.

This capability is relevant for scientific and applied domains in which analysts need to detect recurring local phenomena that are not visible when each series is processed in isolation. By packaging discovery, ranking, and visualization into a single R tool distributed through CRAN, STMotif supports reproducible experimentation, method reuse, and practical adoption of spatial-time motif mining. The associated CSA publication provides direct evidence of research use and scientific relevance.

# AI usage disclosure

Generative AI tools were used only for language support during manuscript adaptation to the JOSS format. The authors reviewed and verified the resulting text, and the technical content, software description, methodological claims, and conclusions remain the responsibility of the authors.

# Acknowledgements

The authors thank CAPES (finance code 001), CNPq, and FAPERJ for partially funding this research.

# References
