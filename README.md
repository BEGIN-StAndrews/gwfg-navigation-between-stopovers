# README: Navigation Between Stopovers by Greater White-Fronted Geese - Code & Workflow

## Overview:
This repository contains the analysis code and workflow for the study "Navigation Between Stopovers by Greater White-Fronted Geese: Comparing Compass Mechanisms and Efficiency Benchmarks". It covers the full pipeline from GPS preprocessing and quality control, through stopover detection and segment construction, to modelling and evaluating alternative navigation strategies. We implement five compass-based routes: geographic loxodrome (GEO), geomagnetic loxodrome (MAG), magnetoclinic (MCL), time-compensated sun compass (SUN), and local wind-aligned (LW); and two efficiency benchmarks: great-circle (GC) and global wind-optimal (GWO); then compare simulated routes against observed tracks using multiple similarity metrics.

## Workflow (run in sequence)
- Preprocess GPS (01)
- Stopovers (02)
- Segmentation & Filtering (03-04)
- Kp & Initial Heading (05-06)
- Geomagnetic Corridors (07_1-07_2)
- Hourly Resample & Kinematics (08)
- Route Simulations (09_1-09_5)
- Route Similarity Measures (10_1-10_2)
- Evaluation Figures (11-13)
- Add Covariates & Route Labels (14_1-14_3)
- Modelling & Repeatability (15_1-15_4)
- Mapping (16_1-16_2)

## Data availability:
- **Observed GPS tracks:** Access may be restricted; final availability and conditions will be stated in the manuscript's Data Availability section after publication.  
- **External datasets:**  
  - **ERA5 hourly wind (u100, v100):** Copernicus Climate Data Store, DOI: 10.24381/cds.adbb2d47.  
  - **Geomagnetic fields:** MagGeo outputs (declination, inclination); details in the manuscript.

## Contact:
For questions or issues, please contact:
- Ali Moayedi
- University of St Andrews, UK
- Email: am636@st-andrews.ac.uk


## Citation:
If you use this code, please cite the associated paper:  

Moayedi, A., et al. "Navigation Between Stopovers by Greater White-Fronted Geese: Comparing Compass Mechanisms and Efficiency Benchmarks", Movement Ecology (2026). DOI: [Pending]  

