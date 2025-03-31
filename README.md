# PFAS Fingerprinting in Massachusetts Surface Waters

This small project explores spatial patterns of PFAS (Per- and Polyfluoroalkyl Substances) detected in surface water across Massachusetts using public data from MassDEP.

We:
- Filtered the dataset to real samples (not QC or duplicates)
- Estimated `<MDL` (below detection) values using ½ the MDL, only when a compound was detected elsewhere
- Averaged PFAS concentrations across sampling sites per waterbody
- Performed hierarchical clustering on PFAS fingerprints
- Visualized the spatial distribution of clusters with an interactive leaflet map

**Data Source**: [MassDEP PFAS Surface Water & Fish Tissue Study](https://www.mass.gov/info-details/pfas-in-surface-water-and-fish-tissue#:~:text=PFAS%20in%20Massachusetts%20Rivers%20(2020)&text=PFAS%20were%20detected%20in%20all%2027%20of%20the%20rivers%20sampled,between%200.3%20and%20399%20ppt.)

**Packages used**: `dplyr`, `tidyr`, `leaflet`, `readxl`, `factoextra`, `ggplot2`
