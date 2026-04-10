# Levy-Lobster-R 

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)  
[![R Version](https://img.shields.io/badge/R-≥4.0-blue.svg)](https://www.r-project.org/)  

## Overview
Levy-Lobster-R is an R-based project that applies Levy processes to ecological modeling, with a focus on lobster population dynamics. It provides scripts for statistical analysis, visualization, and reproducibility testing.

This repository is designed for researchers, students, and data scientists interested in stochastic processes, population modeling, and ecological data analysis.

## Features
- Levy process simulations implemented in R  
- Ecological modeling applied to lobster populations  
- Visualization scripts for generating plots and figures  
- Test scripts to validate functions and workflows  

## Repository Structure
Levy-Lobster-R/
├── scripts/
│   ├── fig3.R
│   ├── fig4.R
│   ├── fig6.R
│   └── fig7.R
│
├── scripts/plots/
│   ├── fig1.R
│   ├── fig2.R
│   ├── fig5.R
│
├── scripts/tables/
│   ├── table1
│   ├── table2_bg_model.R
│   ├── table2_bl_model.R
│   ├── table2_gg_model.R
│   ├── table2_gl_model.R
│   ├── table2_bf.R
│   ├── table2_gf.R
│   ├── table2_igf.R
│   ├── supplementary_bg_model.R
│   ├── supplementary_bl_model.R
│   ├── table3.R
│   ├── table4.R
│   ├── table5_independent_model.R
│   └── table5_joint_models.R
│
├── results/
│   ├── tables/
│   └── plots/
│
└── README.md


## Installation & Usage
Clone the repository and navigate into the folder:

## ```bash
git clone https://github.com/ygwang2018/Levy-Lobster-R.git
cd Levy-Lobster-R

## Install required R packages:
Before running any scripts, make sure the following R packages are installed:

- **dplyr** – data manipulation (filtering, grouping, joins)
- **ggplot2** – plotting and figure generation (used in `scripts/plots/`)
- **MASS** – additional statistical functions (Gamma GLM support)
- **stats** – base R package (provides `glm`, `optim`, distributions)
- **utils** – base R package (provides `write.csv`, etc.)
- **graphics / grDevices** – base R packages (plotting functions)

You can install the non‑base packages with:

install.packages(c(
  "dplyr",
  "ggplot2",
  "MASS"
))

## Run the main analysis:
source("scripts/main_analysis.R")

## Generate plots:
source("scripts/plot_results.R")

## Run example tests:
source("scripts/test_functions.R")

## Results
The analysis produces statistical outputs and visualizations that demonstrate the application of Levy processes to ecological data.
Example outputs can be found in the results/ directory.
(Consider adding sample plots or screenshots here for clarity.)

## Testing
To validate functions and ensure reproducibility, run:
source("scripts/test_functions.R")

## Contributing
Contributions are welcome!
1. Fork the repository
2. Create a feature branch (git checkout -b feature-name)
3. Commit changes (git commit -m "Add feature")
4. Push to the branch (git push origin feature-name)
5. Open a Pull Request
Please follow coding style guidelines and include tests where appropriate.

## License
This project is licensed under the MIT License. See the LICENSE file for details.

## Acknowledgments
• Inspiration from stochastic process modeling in ecology
• R community for packages like tidyverse, ggplot2, and dplyr
• Collaborators and contributors who supported the project
