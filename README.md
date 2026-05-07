# 🌿 Macroevolutionary Patterns of Extinction Risk in Iberian Vascular Plants

This repository contains the complete codebase, datasets, and supplementary materials for the study:

**“Macroevolutionary patterns of extinction risk in Iberian vascular plants”**  


---

## 📁 Repository Structure

The repository is organized into **five main methodological sections**, each containing one or more scripts depending on the complexity of the task. Scripts pull input data from **`Data/Raw`** and generate outputs into **`Data/Processed`**, **`Results`**, or **`Figures`**.

---

### 1. Species List and Extinction Risk Data  
- `01_conservation_data.R`: This script processes and formats data frames for further analyses from two plant extinction risk assessments corresponding to Peninsular Spain and Eastern Andalusia.

---

### 2. Phylogenetic Analyses and EDGE calculations 

This section contains the scripts used to process phylogenetic data, estimate diversification dynamics and phylogenetic structure, and calculate EDGE2 metrics for the floras of Peninsular Spain and Eastern Andalusia. The resulting datasets constitute the main inputs for subsequent comparative and macroevolutionary analyses.

- `02_phylo_manipulation.R`: Processes the global plant phylogeny to generate region-specific genus-level phylogenies through pruning and collapsing procedures. The script also calculates phylogenetic signal in proportion of threatened species within genera and species' probability of extinction, estimates corrected taxon ages under three extinction scenarios, and produces the final data frames used in downstream analyses.
- `02.1_julia_ClaDS.jl`: Estimates genus-level diversification rates from the global plant phylogeny using the ClaDS framework implemented in Julia.
- `02.2_EDGE_metrics.R`: Applies the EDGE2 framework on extinction risk assessments and the prunned phylogenies to estimate extinction probabilities and evolutionary distinctiveness metrics for plant species in Peninsular Spain and Eastern Andalusia.

---

### 3. Statistical Analyses  
**Scripts:** `03_GLM_proportion.R` & `06_Quantile_regs.R`  
Fits GLMs (quasibinomial) to assess relationships between threat proportions and phylogenetic predictors. Also includes quantile regressions to capture non-linear and threshold effects of genus age on extinction risk.

---

### 4. Phylogenetic Diversity Loss  
**Script:** `04_phylo_div_loss_curves.R`  
Simulates the sequential extinction of genera based on threat status and calculates PD loss curves. Results are compared against null expectations from randomized extinction sequences.

---

### 5. Sensitivity Analysis  
**Script:**  `07_sensitivity_randomization.R`  
Tests the robustness of results by randomly reassigning DD and NE species as threatened or non-threatened, depending on the proportions, across 100 replicates.

---

## 🧰 Supporting Scripts

- `functions.R`: Collection of custom functions used across scripts.  
- `05_plots_figures.R`: Script used to generate figures and plots (output in `Figures/`).

---

## 💻 Environment

Main analyses were performed using:  
- **R** with packages: `ape`, `phytools`, `picante`, `DescTools`, `quantreg`, etc.  
- **Julia** with the `PANDA` package for ClaDS implementation.  


---

## 📫 Contact

For questions or collaborations, feel free to contact:

**Carlos Calderón del Cid** – caldecid@gmail.com  
**Bruno Vilela** – bvilela.bv@gmail.com







