# Longitudinal Data Analysis: PBC Case Study

This repository contains the final project for Course STAT 5060, focusing on advanced statistical methods for longitudinal clinical data. The analysis integrates Latent Class Mixed Models (GMM), Missing Not At Random (MNAR) mechanisms, and Multilevel Growth Curve Models.

## 📊 Live Demo

**Interactive Clinical Dashboard**: [https://conniemessi.shinyapps.io/final_project/](https://conniemessi.shinyapps.io/final_project/)

*Note: If the link is inaccessible, you can run the app locally:*
```r
# Ensure packages are installed: shiny, ggplot2, dplyr, tidyr, bslib
shiny::runApp("app.R")
```

## 📂 Project Structure

- **`4_growth.R`**: Multilevel Growth Curve Models (Chapter 4). Includes random intercept/slope models, quadratic growth, and stratified analysis by GMM class.
- **`5_GMM.R`**: Finite Mixture Models (Chapter 5). Identifies latent patient subgroups (latent classes) with distinct disease trajectories.
- **`6_missing.R`**: Missing Data Analysis (Chapter 6). Compares Naive, Heckman Selection, Pattern Mixture, Multiple Imputation, and Common Factor models for handling MNAR data.
- **`app.R`**: Source code for the Shiny dashboard.
- **`pbc_longitudinal.csv`**: The dataset used for analysis (derived from `survival::pbcseq`).
- **`report_comprehensive.pdf`**: The final comprehensive report integrating all analyses.
- **`plots/`, `plots_growth/`, `plots_gmm/`**: Generated visualization outputs.