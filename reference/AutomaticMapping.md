# AutomaticMapping

Function automatically finds fitting mapping for provided MAgPIE object
for a given target aggregation.

## Usage

``` r
AutomaticMapping(x,mapping=NULL,from=NULL,to=NULL)
```

## Arguments

- x:

  MAgPIE object

- mapping:

  A array or data.frame containing a mapping query or mapping name

- from:

  Only required if query is not NULL. Column of the query with original
  dimnames of the incoming dataset

- to:

  Only required if query is not NULL. Column of the query with the
  target dimnames of the outcoming dataset

## Author

Benjamin Leon Bodirsky
