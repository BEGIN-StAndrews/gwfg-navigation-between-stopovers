## Navigation Between Stopovers by Greater White-fronted Geese

*Code and workflow for route-simulation and navigation-mechanism analysis*

---

## Overview

This repository contains the analysis code and workflow for a study of route choice and orientation mechanisms in migratory greater white-fronted geese. The workflow links GPS preprocessing, stopover detection, route-segment construction, route simulation, route-similarity analysis, statistical modelling, and mapping outputs.

The associated manuscript is provisionally titled:

> **What can analysis of routes between stopovers reveal about orientation mechanisms in a long-distance migratory bird?**

The analysis compares observed routes between stopovers with simulated routes generated from compass-based orientation mechanisms and efficiency benchmarks. Compass-based routes include geographic loxodrome (GEO), geomagnetic loxodrome (MAG), magnetoclinic (MCL), time-compensated sun compass (SUN), and local wind-aligned routes (LW). Efficiency benchmarks include great-circle routes (GC) and global wind-optimal routes (GWO).

## Workflow

The scripts are organised to be run in sequence:

1. **GPS preprocessing** (`01`)
2. **Stopover detection** (`02`)
3. **Segmentation and filtering** (`03`–`04`)
4. **Kp index and initial-heading extraction** (`05`–`06`)
5. **Geomagnetic corridors** (`07_1`–`07_2`)
6. **Hourly resampling and movement kinematics** (`08`)
7. **Route simulations** (`09_1`–`09_5`)
8. **Route-similarity measures** (`10_1`–`10_2`)
9. **Evaluation figures** (`11`–`13`)
10. **Covariates and route labels** (`14_1`–`14_3`)
11. **Modelling and repeatability analysis** (`15_1`–`15_4`)
12. **Mapping outputs** (`16_1`–`16_2`; `17`)

## Data availability

- **Observed GPS tracks:** Access may be restricted. Final availability and access conditions will be stated in the manuscript's Data Availability section after publication.
- **ERA5 hourly wind data:** Copernicus Climate Data Store, DOI: `10.24381/cds.adbb2d47`.
- **Geomagnetic fields:** Generated using MagGeo-derived geomagnetic variables; methodological details are provided in the associated manuscript.

## Citation

If you use this code or workflow, please cite the associated manuscript:

Moayedi, A., Long, J. A., Kölzsch, A., Kruckenberg, H., Benítez-Páez, F., & Demšar, U. *What can analysis of routes between stopovers reveal about orientation mechanisms in a long-distance migratory bird?* DOI: Pending.

## Contact

Ali Moayedi  
University of St Andrews, UK  
am636@st-andrews.ac.uk
