# spotr
An R-package for estimating population indices from spatio-temporal models fitted to monitoring data. 
The package works with models fitted by `mgcv` or `brms` and computes aggregated indices of relative population change across space and time.

## Installation

The package is on CRAN and can be installed by
```{r}
install.packages("spotr")
```

Alternatively, the package can be installed from github via
```{r}
library(remotes)
install_github("jknape/spotr")
```
