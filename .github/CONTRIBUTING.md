# Contributing to melt 

The R package **melt** is in a stable state of development, with some degree of subsequent development planned by the authors. 
We welcome your contributions to the package. 
Please submit questions, bug reports, and requests in the [issues](https://github.com/markean/melt/issues).
If you’ve found a bug, please file an issue that illustrates the bug with a minimal [reprex](https://www.tidyverse.org/help/#reprex).
If you want to make a bigger change, it's a good idea to first file an issue and then fork the repository if needed.

## Roadmap
TO DO:

- Additional options of families and link functions for `el_glm()`. 

- Additional model diagnostic tools.

- Additional multiple testing options.

- Additional optimization routines.

NOT TO DO:

- Other parallelization schemes such as [boot](https://cran.r-project.org/package=boot), [RcppParallel](https://cran.r-project.org/package=RcppParallel), etc. 
OpenMP is employed for parallel computing. Some unit tests are also written with the assumption that OpenMP is available, although the availability does not affect the test results. 

- Methods for over-identified models.


## Code style
- We generally follow the tidyverse [style guide](https://style.tidyverse.org). 
You can use the [styler](https://CRAN.R-project.org/package=styler) package to apply these styles. 
All function and variable names should use `snake_case`. 
S4 class names should use `UpperCamelCase`. 
S4 generics and methods should use `lowerCamelCase`. 
These guidelines also apply to C++ code style.

- We use [roxygen2](https://cran.r-project.org/package=roxygen2), 
with [Markdown syntax](https://cran.r-project.org/web/packages/roxygen2/vignettes/rd-formatting.html), 
for documentation.  

- We use [testthat](https://cran.r-project.org/package=testthat) for unit tests. 
Contributions with test cases included are easier to accept.  

## Code of Conduct
Please note that the melt project is released with a [Contributor Code of Conduct](CODE_OF_CONDUCT.md). 
By contributing to this project you agree to abide by its terms.
