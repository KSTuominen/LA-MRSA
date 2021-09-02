## Code for modelling the spread of LA-MRSA in a Swedish pig herd
This repository contains the code used in the manuscript "Modelling the spread of livestock-associated methicillin-resistant Staphylococcus aureus in a Swedish pig herd" by Tuominen et al. The model is written using R language (version 4.0.3).
### Dependencies
Running the model requires installing the SimInf package. The latest released version can be installed from [CRAN](https://cran.r-project.org/web/packages/SimInf/index.html) by
`install.packages("SimInf")`

Running the approximate Bayesian computation (ABC) currently requires a [development version of SimInf](https://github.com/stewid/siminf) which can installed through GitHub by using the remote package. At the time of writing, the version was 8.2.0.9000.
```R
library(remotes)
install_github("stewid/SimInf")
```
### Running the model
Baseline model can be run by executing model.R. This runs the model from day 0 to 3000 and generates events "on the go". The indirect transmission rates used in the model run are the median rates from the medium parameter set that were obtained through parametrization (see *Approximate Bayesian computation (ABC)* below). The model.R produces a data frame "final_result", which has the status of each node (pen) at each day of the model run.

For faster run, code that uses pre-generated events is available in model-premade-events.R. This runs the model from day 731 to 3000, where whole herd has been infected on day 730. The indirect transmission rates and result data frame are the same as in model.R.

### Approximate Bayesian computation (ABC)
The code used for parametrization is available in abc.R. Note that ABC functionality is currently available only in a development version of SimInf (see *Dependencies* above). 

The ABC uses the model with pre-generated events and pre-calculated target parameter sets to generate fit and saves the result graph under /fit. The graphs for the last generations of each target parameter sets that were produced for the manuscript are readily provided as .pdf in the subdirectories of /fit. Similarly, the last generation model fits are provided as .Rda files and can be reloaded by typing `load("fit/SET/fitted.Rda")`, where "SET" is one of the target parameter set subdirectories (low, med or high). 

The loaded fit results can then be viewed by typing `fit`.

### Authors
[![orcid](https://info.orcid.org/wp-content/uploads/2019/11/orcid_16x16.png)](https://orcid.org/0000-0002-6576-9668) Thomas Rosendal, [![orcid](https://info.orcid.org/wp-content/uploads/2019/11/orcid_16x16.png)](https://orcid.org/0000-0002-2223-9376) Krista Tuominen (Maintainer)
