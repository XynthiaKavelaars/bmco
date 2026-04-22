# Round to Chosen Interval

Function to round number up or down to a chosen interval. Adapted from:
https://stackoverflow.com/a/32508105

## Usage

``` r
round_choose(x, round_to, dir = 1)
```

## Arguments

- x:

  Scalar. Number to be rounded.

- round_to:

  Scalar. Interval to be rounded to. E.g. 5, to round to the next 5th
  number.

- dir:

  Integer. "1" for rounding up; "0" for rounding down. Defaults to 1.

## Value

Scalar. Rounded number.
