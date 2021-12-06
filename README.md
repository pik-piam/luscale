# PIK Landuse Group Data Scaling Tools

R package **luscale**, version **2.24.4**

[![CRAN status](https://www.r-pkg.org/badges/version/luscale)](https://cran.r-project.org/package=luscale) [![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.1158584.svg)](https://doi.org/10.5281/zenodo.1158584) [![R build status](https://github.com/pik-piam/luscale/workflows/check/badge.svg)](https://github.com/pik-piam/luscale/actions) [![codecov](https://codecov.io/gh/pik-piam/luscale/branch/master/graph/badge.svg)](https://codecov.io/gh/pik-piam/luscale) [![r-universe](https://pik-piam.r-universe.dev/badges/luscale)](https://pik-piam.r-universe.dev/ui#builds)

## Purpose and Functionality

A collection of tools which allow to aggregate und disaggregate data in various ways.


## Installation

For installation of the most recent package version an additional repository has to be added in R:

```r
options(repos = c(CRAN = "@CRAN@", pik = "https://rse.pik-potsdam.de/r/packages"))
```
The additional repository can be made available permanently by adding the line above to a file called `.Rprofile` stored in the home folder of your system (`Sys.glob("~")` in R returns the home directory).

After that the most recent version of the package can be installed using `install.packages`:

```r 
install.packages("luscale")
```

Package updates can be installed using `update.packages` (make sure that the additional repository has been added before running that command):

```r 
update.packages()
```

## Questions / Problems

In case of questions / problems please contact Jan Philipp Dietrich <dietrich@pik-potsdam.de>.

## Citation

To cite package **luscale** in publications use:

Dietrich J, Bodirsky B, Bonsch M, Patrick, von Jeetze, Kreidenweiss U, Hennig R, Humpenoeder F (2021). _luscale: PIK Landuse Group Data Scaling Tools_. doi: 10.5281/zenodo.1158584 (URL: https://doi.org/10.5281/zenodo.1158584), R package version 2.24.4, <URL: https://github.com/pik-piam/luscale>.

A BibTeX entry for LaTeX users is

 ```latex
@Manual{,
  title = {luscale: PIK Landuse Group Data Scaling Tools},
  author = {Jan Philipp Dietrich and Benjamin Leon Bodirsky and Markus Bonsch and {Patrick, von Jeetze} and Ulrich Kreidenweiss and Roman Julius Hennig and Florian Humpenoeder},
  year = {2021},
  note = {R package version 2.24.4},
  doi = {10.5281/zenodo.1158584},
  url = {https://github.com/pik-piam/luscale},
}
```
