# Moselle_CEFIPRA_projet
This repository consists of all Codes and Data sets for reproducibility for the analysis of Invasion of Round Goby impact in the Moselle River. 
# Moselle_CEFIPRA_projet

This repository contains the data and code used in the study:

**"GNN-based invasion spread and ecological risk assessment in the Moselle River"**

---

## Overview

This project investigates the spread dynamics and ecological impact of invasive species in the Moselle River using an integrated modeling framework combining:

- Hierarchical Modelling of Species Communities (HMSC)
- Graph Neural Networks (GNN)
- Graph Attention Networks (GAT)
- Ecological risk assessment (Vulnerability → Risk → Probability of decline)

The workflow links species-environment relationships, spatial network structure, and community interactions to quantify invasion risk and ecological consequences.

---

## Repository Structure
Moselle_CEFIPRA_projet/
├── README.md
├── data/
│ ├── raw/ # Raw datasets (if included)
│ └── processed/ # Cleaned datasets used in analyses
├── code/
│ ├── R/ # HMSC modeling and statistical analyses
│ └── Python/ # GNN, GAT, and risk modeling scripts
├── figures/ # Figures generated from the analysis (optional)
├── results/ # Intermediate outputs (optional)

---

## Methods Summary

### 1. Data Preparation
- Fish community data aggregated at station × campaign × year level  
- Environmental covariates standardized prior to modeling  

### 2. HMSC Modeling (R)
- Estimation of species–environment relationships (β)
- Trait–environment interactions (Γ)
- Residual species associations (Ω)

### 3. Network-Based Modeling (Python)
- GNN used to estimate invasion spread probability  
- GAT used to quantify:
  - Node centrality (cᵢ)
  - Cascade impact (Iᵢ)

### 4. Impact Assessment
- Intrinsic and contextual vulnerability combined  
- Structured risk index calculated  
- Probability of species decline (P_decline) estimated  

---

## Reproducibility

The scripts are organized to follow the full analysis pipeline:

1. Data preprocessing  
2. HMSC model fitting  
3. GNN/GAT modeling  
4. Risk and decline computation  
5. Figure generation  

All processed datasets required to reproduce the analyses are included.

---

## Requirements

### R
- Hmsc
- coda
- tidyverse  
- other standard statistical packages  

### Python
- numpy  
- pandas  
- torch  
- torch_geometric  
- matplotlib / seaborn
- other standard python libraries

---

## Data Availability Note

This repository includes all data and code necessary to reproduce the analyses presented in the manuscript.

Processed datasets used in the analyses are provided. Raw data are not publicly shared due to data-sharing constraints.

---

## Status

This repository accompanies a manuscript currently under review. Data and code will be made fully public upon acceptance.

---

## Contact

For questions regarding the code or data, please contact the authors.
