# superAggregate

Function which applies aggregation functions over a subset of an magpie
or lpj object.

## Usage

``` r
superAggregate(data, aggr_type, level="reg", weight=NULL, na.rm=TRUE,
crop_aggr=FALSE, ...)
```

## Arguments

- data:

  An MAgPIE or LPJ object

- aggr_type:

  Aggregation Type. Can be any function for one or two objects (data and
  weight) of the same size. Currently pre-supported functions:
  "sum","mean","weighted_mean".

- level:

  Aggregation level: Either a level type (as a name) or a vector of
  3-Character region names with the same length as the data has cells.
  Allowed level types are global "glo", regional "reg", per Country
  "country" and per REMIND regions "remind_reg". "country" and
  "remind_reg" are only supported for 0.5 grid datacells with 59199
  cells and are always returned as arrays. If you use a vector of
  regions the aggregation will take place according to your regions.

- weight:

  Currently only used for weighted_mean (see
  [`weighted.mean`](https://rdrr.io/r/stats/weighted.mean.html), yet
  also applicable for individualized functions. Has to be of the same
  size as data.

- na.rm:

  If TRUE, NAs are ignored both in data and weight.

- crop_aggr:

  determines whether output should be crop-specific (FALSE) or
  aggregated over all crops (TRUE). The method used for aggregation is
  set by aggr_type (Currently works only for levels "reg" and "glo")

- ...:

  additional arguments for the aggregation method for the standard
  functions (not for self-created ones)

## Value

In the case of level="glo" or "reg", the function returns a MAgPIE
object. In the case of level="country" an array is returned. In the case
of an LPJ object being aggregated, a list of MAgPIE objects is returned,
each entry being one of the 4th dimension slices.

## See also

[`colSums`](https://rdrr.io/r/base/colSums.html)

## Author

Benjamin Bodirsky, Jan Philipp Dietrich, Florian Humpenoeder

## Examples

``` r
data(population_magpie)
superAggregate(population_magpie,"sum",level="glo")
#>        scenario
#> t              A2       B1
#>   y1995  5472.867 5472.867
#>   y2005  6470.020 6471.290
#>   y2015  7349.480 7251.820
#>   y2025  8263.170 7899.360
#>   y2035  9130.160 8356.570
#>   y2045  9896.750 8616.840
#>   y2055 10555.010 8684.690
#>   y2065 11128.330 8564.640
#>   y2075 11619.640 8294.200
#>   y2085 12024.700 7891.590
#>   y2095 12294.300 7356.890
#>   y2105 12386.320 7056.120
#>   y2115 12386.320 7056.120
#>   y2125 12386.320 7056.120
#>   y2135 12386.320 7056.120
#>   y2145 12386.320 7056.120
superAggregate(population_magpie,"mean",level="glo")
#>        scenario
#> t              A2       B1
#>   y1995  547.2867 547.2867
#>   y2005  647.0020 647.1290
#>   y2015  734.9480 725.1820
#>   y2025  826.3170 789.9360
#>   y2035  913.0160 835.6570
#>   y2045  989.6750 861.6840
#>   y2055 1055.5010 868.4690
#>   y2065 1112.8330 856.4640
#>   y2075 1161.9640 829.4200
#>   y2085 1202.4700 789.1590
#>   y2095 1229.4300 735.6890
#>   y2105 1238.6320 705.6120
#>   y2115 1238.6320 705.6120
#>   y2125 1238.6320 705.6120
#>   y2135 1238.6320 705.6120
#>   y2145 1238.6320 705.6120
superAggregate(population_magpie,"weighted_mean",level="glo",weight=population_magpie)
#>        scenario
#> t              A2        B1
#>   y1995  817.3915  817.3915
#>   y2005  942.7120  934.8523
#>   y2015 1080.3979 1042.6060
#>   y2025 1229.7711 1135.2091
#>   y2035 1377.6027 1199.1047
#>   y2045 1515.4201 1234.2838
#>   y2055 1639.4112 1241.8551
#>   y2065 1749.7419 1219.3721
#>   y2075 1840.1164 1170.1472
#>   y2085 1905.7570 1096.3516
#>   y2095 1937.0978 1000.1402
#>   y2105 1943.0084  947.2417
#>   y2115 1943.0084  947.2417
#>   y2125 1943.0084  947.2417
#>   y2135 1943.0084  947.2417
#>   y2145 1943.0084  947.2417
aggregation_function<-function(func_data,func_weight) {
  colMeans(func_data)
}
superAggregate(population_magpie,aggregation_function,level="glo",weight=population_magpie)
#>        scenario
#> t              A2       B1
#>   y1995  547.2867 547.2867
#>   y2005  647.0020 647.1290
#>   y2015  734.9480 725.1820
#>   y2025  826.3170 789.9360
#>   y2035  913.0160 835.6570
#>   y2045  989.6750 861.6840
#>   y2055 1055.5010 868.4690
#>   y2065 1112.8330 856.4640
#>   y2075 1161.9640 829.4200
#>   y2085 1202.4700 789.1590
#>   y2095 1229.4300 735.6890
#>   y2105 1238.6320 705.6120
#>   y2115 1238.6320 705.6120
#>   y2125 1238.6320 705.6120
#>   y2135 1238.6320 705.6120
#>   y2145 1238.6320 705.6120
```
