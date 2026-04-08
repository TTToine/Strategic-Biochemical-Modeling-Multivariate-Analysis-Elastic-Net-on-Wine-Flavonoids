# Strategic Biochemical Modeling: Multivariate Analysis & Elastic Net on Wine Flavonoids

## Project Overview
This project explores the chemical profile of three wine cultivars using the famous `wine` dataset. Rather than focusing on simple classification, this study aims to reverse-engineer the biochemical ecosystem. By treating **Flavonoids** (a vital parameter for wine structure, aging potential, and antioxidant capacity) as the target variable, this project maps the metabolic drivers and competitive pathways within the wine matrix.

The ultimate goal is to bridge the gap between advanced statistical inference and industrial operations, translating regression coefficients into actionable business and agronomic levers.

## Objectives
* Overcome the limits of univariate analysis by modeling a complex, highly multicollinear biological matrix.
* Isolate true causal biochemical drivers from spurious correlations using penalized selection criteria (AIC/BIC).
* Handle structural multicollinearity natively using **Elastic Net Regularization**.
* Map the conditional dependencies of the system to visualize metabolic competition.
* Provide data-driven strategic recommendations for winery operations.

## Methodology & Tech Stack
* **Language & Environment:** R, RStudio
* **Exploratory Data Analysis (EDA):** Density plots, LOESS curves, and bivariate analysis to identify the baseline additive structure driven by the cultivar.
* **Model Selection:** Stepwise feature selection comparing pure additive, polynomial, and fully saturated spaces. **BIC (Bayesian Information Criterion)** was favored over AIC to strictly penalize complexity and extract the true chemical backbone of the system.
* **Regularization:** Implementation of **Elastic Net (L1 + L2 penalty)** to stabilize the model against the endogenous multicollinearity (e.g., the structural co-dependence between Phenols and Flavonoids).
* **Network Analysis:** Computation of Partial Correlations (Gaussian Graphical Models logic) to map conditional dependencies and mathematically isolate the metabolic competition.

## Key Business & Operational Insights
The statistical output was translated into two main strategic levers for a winery:
1. **Upstream Strategy (Agronomic Selection):** The choice of the cultivar (specifically `v1`) acts as a biological switch. It establishes an insurmountable baseline advantage in Flavonoid abundance that cannot be replicated downstream.
2. **Downstream Operations (Winery Extraction):** The negative elasticity of *Non-Flavanoids* reveals a strong metabolic competition. Extraction technologies in the cellar must be calibrated not only to push phenolic extraction but to **actively inhibit** the collateral synthesis of non-flavonoid compounds.

## Repository Structure
* `scripts/`: R scripts containing the data pipeline, EDA, Elastic Net modeling, and Network generation.
* `report/`: The comprehensive 21-page analytical report (PDF).
* `presentation/`: The final pitch deck created in LaTeX Beamer, tailored for a managerial/academic audience.
* `plots/`: High-resolution outputs of the network graphs and diagnostic plots.

## How to Run
1. Clone this repository: `git clone https://github.com/YourUsername/wine-flavonoids-multivariate.git`
2. Open the R project.
3. Ensure the required packages are installed (`glmnet`, `MASS`, `car`, `qgraph`, etc.).
4. Run `main_analysis.R` to reproduce the models and visual outputs.

---
*Author: Antoine Peluso* *MSc in Management Engineering*
