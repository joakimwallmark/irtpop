# irtpop
An R package for computing population-based item difficulty and discrimination metrics for IRT models fitted with the [mirt package](https://cran.r-project.org/web/packages/mirt/index.html). Unlike traditional IRT parameters, these metrics integrate over a target population distribution, giving a single summary of how difficult and discriminating each item is for a specific population.

The package can be installed directly from GitHub using [pak](https://pak.r-lib.org/):

```R
# install.packages("pak")
pak::pak("joakimwallmark/irtpop")
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
data <- expand.table(LSAT7)
mod_se <- mirt(data, 1, itemtype = "2PL", verbose = FALSE, SE = TRUE)
item_metrics(mod_se, se = TRUE)
```

