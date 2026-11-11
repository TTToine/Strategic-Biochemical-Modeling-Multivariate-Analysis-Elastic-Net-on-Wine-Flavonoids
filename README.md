# Wine Flavonoids: Model Selection, Elastic Net and Partial-Correlation Networks

Multivariate analysis of the chemical profile of three wine cultivars, with **flavanoid
content** as the response variable. The project compares stepwise model selection under
AIC and BIC, stabilises estimation under multicollinearity with an **Elastic Net**, and
maps **conditional (partial) correlations** between compounds as a network.

## Motivation

Flavanoids drive wine structure, ageing potential and antioxidant capacity. The chemical
variables in the `wine` dataset are strongly collinear (e.g. total phenols and
flavanoids), so marginal correlations are a poor guide to which compounds are associated
with flavanoid content *once the others are held fixed*. The goal is to separate those
conditional associations from redundant marginal ones and to translate them into
hypotheses for vineyard and cellar operations.

## Data

The classic `wine` dataset (178 samples, 3 cultivars, 13 chemical measurements), loaded
from the R package [`gRbase`](https://cran.r-project.org/package=gRbase) with `data(wine)`.
Variables are renamed in the script (`Type`, `Alcohol`, `Malic`, …, `Flavanoids`, …, `Proline`).

## Method

1. **Exploratory analysis**: descriptive statistics by cultivar, density and box plots,
   LOESS curves of each variable against flavanoids, correlation matrices by cultivar.
2. **Model selection**: linear models for `Flavanoids` compared across three model spaces
   (additive; additive + selected quadratic terms; all pairwise interactions), each reduced
   by bidirectional stepwise search under **AIC** (k = 2) and **BIC** (k = log n).
   BIC is preferred to favour a parsimonious model; the additive BIC model is retained.
3. **Diagnostics of the selected model**: VIF, coefficient estimates with 95% CIs,
   fitted vs observed, residual diagnostics and partial-residual plots.
4. **Elastic Net** (`glmnet`, α = 0.5, 10-fold CV, `lambda.min`) as a regularised
   check on the stepwise selection under multicollinearity.
5. **Partial-correlation network**: partial correlations from the inverse correlation
   matrix (Gaussian graphical model logic), thresholded at |r| > 0.15 and drawn with `qgraph`.

## Main findings

- **Cultivar is the dominant factor**: it sets a baseline level of flavanoids that the
  other chemical variables only modulate.
- **Non-flavanoid phenols enter with a negative conditional association**, consistent with
  a competition between phenolic pathways.

These are conditional associations in observational data, not causal effects. Read as
hypotheses for operations: cultivar choice is an upstream lever, and cellar extraction
could be tuned to favour flavanoid over non-flavanoid phenols. Testing either would
require controlled experiments.

## Limitations

- 178 observations; stepwise selection is unstable and its p-values are optimistic after selection.
- The network uses Pearson partial correlations, which assume roughly Gaussian, linear relations.
- No out-of-sample validation of the selected model.

## Reproduce

Requires R (≥ 4.1) and the packages:

```r
install.packages(c("gRbase", "ggplot2", "gridExtra", "knitr", "broom", "car",
                   "ggcorrplot", "qgraph", "dplyr", "tidyr", "reshape2", "glmnet"))
```

Then run from the repository root:

```bash
Rscript Vin.R
```

Figures are written to `grafici_wine/`; model comparison tables are printed to the console.

## Repository structure

```
Vin.R       # full pipeline: EDA, model selection, Elastic Net, diagnostics, network
README.md
```

## Author

Antoine Peluso, MSc in AI Engineering, University of Florence.
