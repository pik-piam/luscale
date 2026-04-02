# getAggregationMatrix

Function which converts the supplied regionmapping file to a
transformation matrix which then can be used for aggregation with
[`speed_aggregate`](speed_aggregate.md).

## Usage

``` r
getAggregationMatrix(rel, from=NULL, to=NULL)
```

## Arguments

- rel:

  file name of a region mapping (".csv" as file ending, ";" used as
  separator and with columns containing the region names) or a mapping
  matrix (matrix with region names in columns)

- from:

  Name of the first column to be used (if not set the first or second
  column will be used)

- to:

  Name of the second column to be used (if not set the second or third
  column will be used)

## Value

A matrix nregions1 x nregions2 with 1s and 0s showing the mapping of
countries to regions.

## Author

Jan Philipp Dietrich

## Examples

``` r
if (FALSE) { # \dontrun{
  x <- cbind(reg=c("REG1","REG2","REG2","REG3"),country=c("C1","C2","C3","C4")) 

  getAggregationMatrix(x)
  
  getAggregationMatrix(x,from="reg",to="reg")
} # }
```
