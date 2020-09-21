# PIK Landuse Group Data Scaling Tools

R package **luscale**, version **2.18.2**

[![Travis build status](https://travis-ci.com/pik-piam/luscale.svg?branch=master)](https://travis-ci.com/pik-piam/luscale) [![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.1158584.svg)](https://doi.org/10.5281/zenodo.1158584) 

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

Dietrich J, Bodirsky B, Bonsch M, Kreidenweiss U, Hennig R, Humpenoeder F (2020). _luscale: PIK Landuse Group Data
Scaling Tools_. doi: 10.5281/zenodo.1158584 (URL: https://doi.org/10.5281/zenodo.1158584), R package version 2.18.2,
<URL: https://github.com/pik-piam/luscale>.

A BibTeX entry for LaTeX users is

 ```latex
@Manual{,
  title = {luscale: PIK Landuse Group Data Scaling Tools},
  author = {Jan Philipp Dietrich and Benjamin Leon Bodirsky and Markus Bonsch and Ulrich Kreidenweiss and Roman Julius Hennig and Florian Humpenoeder},
  year = {2020},
  note = {R package version 2.18.2},
  doi = {10.5281/zenodo.1158584},
  url = {https://github.com/pik-piam/luscale},
}
```

