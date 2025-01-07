
### Bagging ridge model
If user sought to use CellOracle, another method for GRN inference,
It uses Bagging ridge and Bayesian ridge regression models from sklearn (python). 
scMore leverages Pando-based apporach taht use reticulate to interact with python and implement these models also here. 
Thus, user should install scikit-learn and pandas in python in advance.

## MAC or Linux system

open terminal


```bash

source ~/.virtualenvs/r-reticulate/bin/activate

pip install pandas
Successfully installed pandas-<version>

pip install scikit-learn
Successfully installed scikit-learn-<version>


```


## Open R
### make sure that 'reticulate' uses the vitual environment
```r
library(reticulate)
use_virtualenv("~/.virtualenvs/r-reticulate", required = TRUE)
pd <- import("pandas")
sklearn <- import("sklearn")
