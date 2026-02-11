# EmpStudyProject

## Overview
The EmpStudyProject is designed for robust change point detection in time series data, specifically focusing on temperature anomalies. The project implements various statistical methods to identify significant changes in trends over time.

## Files
- **EmpStudy.r**: Contains the main R code for robust change point detection. It includes functions for input validation, candidate change point generation, test statistics computation, and threshold selection for controlling false discovery rate.
  
- **data/temperature_anomalies.RData**: Contains raw temperature anomaly data used in the analysis.
  
- **data/temperature_anomalies_adj.RData**: Contains adjusted temperature anomaly data used in the analysis.
  
- **MethodCode/PELTtrendARp.R**: Implements the PELT (Pruned Exact Linear Time) algorithm for detecting change points in time series data.

## Usage
1. Load the necessary R packages.
2. Source the `EmpStudy.r` file to access the functions.
3. Load the data files as needed for analysis.
4. Call the `robust_change_point_detection` function with the appropriate parameters to perform change point detection.

## Dependencies
- R (version 4.0 or higher)
- Required R packages: `strucchange`, `zoo`, `WeightedPortTest`

## Setup Instructions
1. Clone the repository to your local machine.
2. Ensure that all dependencies are installed.
3. Load the data files and run the analysis as described in the usage section.