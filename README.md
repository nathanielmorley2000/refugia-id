# Overview

This data and code repository is supplemental to the manuscript "Ice age
refugium shows potential for geohistorical data to guide modern
conservation efforts" currently under review at *Communications Earth &
Environment*. All data files were pulled from the Neotoma database using
scripts contained herein. Geochronologic controls were pulled on January
12, 2025, and pollen data were last pulled on January 1, 2026.

All code is self-contained and do not require additional materials from
outside of this repository.


# Software Requirements

This code was written using the programming language `R` v. 4.3.3. The
following R packages are required to execute the code: 
* `Bchron` 
* `dplyr` 
* `gganimate` 
* `ggplot2` 
* `gifski` 
* `neotoma2` 
* `raster` 
* `rnaturalearth` 
* `rnaturalearthdata` 
* `sf` 
* `stringr`
* `tidyr` 
* `viridis`

The full references for these packages can be found in Supplementary
Table 1 of the manuscript.


# Repository Guide

Three files are included in the root of the repository to help with organization and documentation: 
* README.md - This document, a brief overview of the repository and the material it contains. 
* .gitignore - A list of files and file types to be excluded from the GitHub repository (mostly
temporary files). 
* refugia-id.Rproj - R project file to organize repository.

Three `R` scripts are also included in the root of the repository to execute various parts of the analysis: 
* time-calibration.R - A script for pulling geochrnologies and pollen data from Neotoma database,
performing the necessary time calibrations, and organizing the data in a way that is useful for the analysis. 
* animations.R - A script to create the animated heat maps and stills necessary seen in Fig. 4 and
Supplementary Figs. 1 and 2. 
* regional-analytics.R - A script for the remaining analyses, including the temporal and spatial Monte Carlo
analysis, time series shown in Figs. 3-4 and Supplementary Figs. 1-2, and correlations between minimum and maximum sampling protocols.

The *TempFiles* folder contains temporary files that should be reserved for analysis in later scripts.

The *IndividualSummaries* folder contains .csv files for each of the three taxa over the past 20 ka, organized according to the minimum
and maximum sample protocols (see Methods of manuscript for details).

The *Results* folder contains results needed for figures. These include the Picea stills that are used in 
Fig. 4, the Populus stills needed for Supplementary Figs 1 and 2, the time series used in Figs. 3-4 and Supplementary Figs. 1-2,
and the animated heat maps used in Supplementary videos 1-3. Note that any statistical outputs will be entered in your console when running
regional-analytics.R and were directly input into the manuscript.


# Contact

If there are any issues with this code, please contact Nathaniel Morley
via email at [nmorley\@ualberta.ca](mailto:nmorley@ualberta.ca).
