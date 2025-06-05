# 🌿 Macroevolutionary Patterns of Extinction Risk in Iberian Vascular Plants

This repository contains the complete codebase, datasets, and supplementary materials for the study:

**“Macroevolutionary patterns of extinction risk in Iberian vascular plants”**  


---

## 📁 Repository Structure

The repository is organized into **five main methodological sections**, each containing one or more scripts depending on the complexity of the task. Scripts pull input data from **`Data/Raw`** and generate outputs into **`Data/Processed`**, **`Results`**, or **`Figures`**.

---

### 1. Species List and Extinction Risk Data  
**Script:** `01_conservation_data.R`  
Processes and formats extinction risk data for Peninsular Spain and Eastern Andalusia.

---

### 2. Phylogenetic Analyses  
**Scripts:** `02_phylo_manipulation.R` & `02.1_julia_ClaDS.jl`  
Generates region-specific genus-level phylogenies by pruning and collapsing a global phylogeny. Calculates phylogenetic signal in threat proportions and derives four genus-level predictors: species richness, corrected taxon age, speciation rate (via ClaDS), and evolutionary distinctiveness.

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







