# orexin-ligand-virtual-screening
A five-step computational pipeline for virtual screening of orexin  receptor ligands, combining PubChem similarity-based library construction,  physicochemical filtering, ADMET prediction, CNS-weighted composite  scoring, and dual-receptor molecular docking against OX1R and OX2R  crystal structures.
---

## Prerequisites

### Python Environment
All Python notebooks were developed using **Python 3.11** in a conda 
environment. Create and activate the environment before running:

```bash
conda create -n orexin_env python=3.11
conda activate orexin_env
```

### Python Packages
Install required packages:

```bash
pip install pubchempy pandas numpy requests admet-ai
conda install -c conda-forge rdkit openbabel
```

### R Environment
Step 5 requires **R version 4.5.0** or later. Install required packages 
by running this in R or RStudio:

```r
install.packages(c(
  "tidyverse", "ggpubr", "factoextra", "corrplot",
  "ggrepel", "patchwork", "rstatix"
))
```

### External Tools (Step 4 only)
Step 4 requires two command-line tools installed separately:

**Open Babel v3.1.0**
```bash
conda install -c conda-forge openbabel
```
Verify: `obabel --version`

**AutoDock Vina v1.2.7**  
Download the executable for your operating system from:  
https://github.com/ccsb-scripps/AutoDock-Vina/releases

- **Windows:** download `vina_1.2.7_win.exe`, rename to `vina.exe`, 
  place in a folder on your system PATH
- **Mac (Apple Silicon):** download `vina_1.2.7_mac_aarch64`
```bash
  chmod +x vina_1.2.7_mac_aarch64
  mv vina_1.2.7_mac_aarch64 ~/bin/vina
```
- **Linux:** download `vina_1.2.7_linux_x86_64`
```bash
  chmod +x vina_1.2.7_linux_x86_64
  sudo mv vina_1.2.7_linux_x86_64 /usr/local/bin/vina
```
- **Google Colab:**
```python
  !wget https://github.com/ccsb-scripps/AutoDock-Vina/releases/download/v1.2.7/vina_1.2.7_linux_x86_64
  !chmod +x vina_1.2.7_linux_x86_64
  !mv vina_1.2.7_linux_x86_64 /usr/local/bin/vina
```

Verify: `vina --version`

---

## Running the Pipeline

Run the notebooks in order. Each step reads from the previous step's 
output CSV file.

### Step 1 — Ligand Library Construction
**Notebook:** `step1_build_library.ipynb`  
**Input:** None (queries PubChem directly)  
**Output:** `orexin_ligand_library_raw.csv`  
**Environment:** `orexin_env`  
**Description:** Retrieves seed compounds (daridorexant, lemborexant, 
suvorexant, oveporexton) from PubChem and performs 2D fingerprint 
similarity searches (Tanimoto threshold 0.85, 20 results per seed) 
to build an 80-compound ligand library.

### Step 2 — Library Cleaning and ADMET Annotation
**Notebook:** `step2_clean_annotate.ipynb`  
**Input:** `orexin_ligand_library_raw.csv`  
**Output:** `orexin_ligand_library_clean.csv`  
**Environment:** `orexin_env`  
**Description:** Removes salt forms, applies Lipinski Rule-of-5 
filtering, predicts 13 ADMET endpoints using admet-ai, and computes 
CNS-specific drug-likeness scores.

### Step 3 — Screening and Ranking
**Notebook:** `step3_screen_rank.ipynb`  
**Input:** `orexin_ligand_library_clean.csv`  
**Output:** `orexin_library_ranked.csv`, `orexin_shortlist_for_docking.csv`  
**Environment:** `orexin_env`  
**Description:** Applies hard physicochemical filters, computes a 
weighted composite CNS drug-likeness score across eight properties, 
and selects the top 20 compounds for docking.

### Step 4 — Molecular Docking
**Notebook:** `step4_docking.ipynb`  
**Input:** `orexin_shortlist_for_docking.csv`  
**Output:** `final_docking_results.csv`  
**Environment:** `docking_env` (Python 3.11 with Open Babel and AutoDock Vina)  
**Prerequisites:** Open Babel and AutoDock Vina must be installed 
before running this notebook (see Prerequisites above)  
**Description:** Downloads OX1R (PDB: 4ZJC) and OX2R (PDB: 4S0V) 
crystal structures, prepares receptor and ligand files, and runs 
AutoDock Vina docking calculations for all 20 shortlisted compounds 
against both receptors.

### Step 5 — Statistical Analysis and Figures
**Script:** `step5_statistics.R`  
**Input:** `final_docking_results.csv`, `orexin_library_ranked.csv`  
**Output:** Figures 2, 3; Tables 2, 3  
**Environment:** R 4.5.0 with required packages (see Prerequisites)  
**Description:** Generates publication-quality figures including box 
plots, selectivity scatter plot, PCA, and Spearman correlation matrix. 
Performs Wilcoxon rank-sum tests and Spearman correlation tests.

---

## Data Files

| File | Description | Compounds |
|------|-------------|-----------|
| `orexin_ligand_library_raw.csv` | Raw library from PubChem | 80 |
| `orexin_ligand_library_clean.csv` | Cleaned and ADMET-annotated library | 79 |
| `orexin_library_ranked.csv` | Full ranked library with composite scores | 79 |
| `orexin_shortlist_for_docking.csv` | Top 20 candidates for docking | 20 |
| `final_docking_results.csv` | Docking scores against OX1R and OX2R | 20 |

---

## Key Findings

- 78 of 79 compounds (98.7%) satisfied Lipinski Rule-of-5 criteria
- Lemborexant scaffold analogs dominated the top-ranked shortlist
- CID 56944043 (Compound 38, Yoshida et al. 2015) achieved the most 
  balanced dual-receptor docking profile (OX1R −8.796, OX2R −7.780 
  kcal/mol) with published experimental Ki values of 2 nM (OX1R) 
  and 4 nM (OX2R)
- Significant positive correlation between hERG liability and OX1R 
  docking score (ρ = 0.58, p = 0.008)
- Significant negative correlation between CNS drug-likeness and 
  OX1R docking score (ρ = −0.54)

---

## Citation

If you use this code or data, please cite:

Erskine, H. (2025). Virtual screening of an orexin receptor ligand 
library identifies lemborexant analogs with favorable CNS drug-likeness 
and dual-receptor binding profiles. [Preprint]. bioRxiv.

---

## License

This project is released under the MIT License.

---

## Contact

H.P. Erskine  
University of Wisconsin
erskinehp03@uww.edu
