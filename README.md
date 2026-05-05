# irtpop
An R package for computing population-based item difficulty and discrimination metrics for IRT models fitted with the [mirt package](https://cran.r-project.org/web/packages/mirt/index.html). Unlike traditional IRT parameters, these metrics integrate over a target population distribution, giving a single summary of how difficult and discriminating each item is for a specific population.

The package can be installed directly from GitHub using the devtools package as follows:
```R
devtools::install_github("joakimwallmark/irtpop")
```
If you don't have devtools installed yet, you can install it via:
```R
install.packages("devtools")
```
After installing, refer to the package help files for more information on how to use the package:
```R
library(irtpop)
?population_difficulty
?population_discrimination
?item_metrics
```

## Example

```R
library(irtpop)
library(mirt)

# Fit a graded response model to polytomous data
mod <- mirt(Science, 1, itemtype = "graded", verbose = FALSE)
item_metrics(mod)

# Fit a 2PL model to binary data and compute metrics with standard errors
mod_se <- mirt(LSAT7, 1, itemtype = "2PL", verbose = FALSE, SE = TRUE)
item_metrics(mod_se, se = TRUE)
```
