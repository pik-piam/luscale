# rename_dimnames

Renames the dimnames of an array or MAgPIE object after a query.

## Usage

``` r
rename_dimnames(data,dim=1,query=NULL,from=NULL,to=NULL)
```

## Arguments

- data:

  Array

- dim:

  The dimension to be renamed.

- query:

  If NULL, query is automatically searched for. Otherwhise an array,
  data.frame or the path of a csv with at least two columns. One column
  has to have the name of "from", the other one the name of "to". Some
  queries can be found in the svn-folder tools/queries.

- from:

  Only required if query is not NULL. Column of the query with original
  dimnames of the incoming dataset

- to:

  Only required if query is not NULL. Column of the query with the
  target dimnames of the outcoming dataset

## Value

An array with different dimnames

## Note

translate.with.query has the same functionality, is more efficient, yet
more complicated to use.

## Author

Benjamin Bodirsky, Ulrich Kreidenweis
