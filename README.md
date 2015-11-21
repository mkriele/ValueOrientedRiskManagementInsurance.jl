# ValueOrientedRiskManagementInsurance

[![Build Status](https://travis-ci.org/mkriele/ValueOrientedRiskManagementInsurance.jl.svg?branch=master)](https://travis-ci.org/mkriele/ValueOrientedRiskManagementInsurance.jl)

**This package has been tested with Julia 0.4.x only.  Other versions of Julia
may generate errors or warnings.**

The package provides example calculations for the second edition
of the book

Kriele M. and Wolf, J. _Wertorientiertes Risikomanagement von
 Versicherungsunternehmen_, 2nd edition, Springer-Verlag, Berlin Heidelberg,
 2016 (to be published)

It is also intended to use these examples for future editions of the English
 translation of this book,
 _Value-Oriented Risk Management of Insurance Companies_.  (The examples in the
  current English edition are written in R).

 The package consists of 4 distinct illustrations:

 - An extremely simplified example of the SST (Swiss Solvency
   Test) calculation for life insurance
 - A simplified example of the S2 (Solvency 2) calculation for
   non-life insurance
 - A simplified example of the S2 calculation for life insurance
 - An extremely simplified example of an internal economic capital
   model for non-life insurance. In the book, this model is used
   to illustrate concepts in value based management.

## SST Life Calculation

SST stands for "Swiss Solvency Test", which is the Swiss regulatory capital
requirement.  The resulting monetary requirement is referred
to as the target capital `ZK`.

The files *SST_Types.jl* and *SST_Functions.jl* contain types and functions
which can be used to model a *simplified* version of the
the SST standard life model. The calculation is basically as follows. The  change of the
"risk bearing capital",
`ΔRTK`<sup>1</sup> over 1 year is approximated by a quadratic function of the
 risk factors. As the risk factors are assumed to be multinormally distributed,
 the distribution of `ΔRTK` is known.  The target capital `ZK` is the sum of the 99% expected shortfall and the market value margin.  The market value margin is calculated using a cost of capital method, where the capital is given by the 99% expected shortfall of `ΔRTK` and the cost of capital factor is assumed to be 6%.

 The example calculation provided here is based on an extremely simplified life insurance portfolio.

## Solvency 2 Life Calculation

## Solvency 2 Non-Life Calculation

## Internal Economic Capital Model

The internal economic capital model is a Monte Carlo model of an extremely simplified non-life insurance company.  Reserves are completely ignored in the model.



## Footnotes

<sup>1</sup> The abbreviation `RTK` stems from the  original German term "Risikotragendes Kapital". Observe that `RBC` is usually understood to mean risk based capital which is different from risk bearing capital.  
