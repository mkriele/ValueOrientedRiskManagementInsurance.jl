# ValueOrientedRiskManagementInsurance

[![Build Status](https://travis-ci.org/mkriele/ValueOrientedRiskManagementInsurance.jl.svg?branch=master)](https://travis-ci.org/mkriele/ValueOrientedRiskManagementInsurance.jl)

**This package requires Julia 0.4.x.**

***These files are part of branch `de_2ed`
(relating to the 2nd edition of the German version of the book)***


The package provides example calculations for the second edition of the book

> Kriele M. and Wolf, J. _Wertorientiertes Risikomanagement von  Versicherungsunternehmen_, 2nd edition, Springer-Verlag, Berlin Heidelberg,  2016 (to be published)

It is also intended to use these examples for future editions of the English  translation of this book,  _Value-Oriented Risk Management of Insurance Companies_.  (The examples in the current English edition are written in R).

 The package consists of 4 distinct parts:

 - *SSTLife*: An extremely simplified example of the SST (Swiss Solvency Test) calculation for life insurance. The Swiss Solvency Test is the Swiss regulatory capital requirement.  The resulting monetary requirement is referred to as the "target capital" `ZK`<sup>1</sup>.
 - *S2Life*: A simplified example of the S2 (Solvency 2) calculation for non-life insurance. Solvency 2 is the new regulatory capital requirement in the European Union. The resulting monetary requirement is referred to as the "Solvency capital requirement" `SCR`.
 - *S2NonLife*: A simplified example of the S2 calculation for life insurance
 - *ECModel*: An extremely simplified example of an internal economic capital model for non-life insurance. This model is used to illustrate some techniques used in value based management.

Note that we have simplified and (in part changed for our exposition) the regulatory requirements for SST and Solvency 2. Also note that the implementation of Solvency 2 may be slightly different in different EU countries. For definitive information about SST or Solvency 2, please consult the original literature and any guidance issued by the supervisory authorities in the jurisdiction of interest.

## Files in the Folder "src" and its Subfolders

The folder "src" contains types and functions which are meant to be used in specific examples

### SST Life Calculation

The files *SST__Types.jl* and *SST_Functions.jl* contain types and functions which can be used to model a simplified version of the the SST standard life model. The calculation is basically as follows. The change of the "risk bearing capital", `ΔRTK`<sup>2</sup> over 1 year is approximated by a quadratic function of the risk factors. As the risk factors are assumed to be multinormally distributed, the distribution of `ΔRTK` is known.  The target capital `ZK` is the sum of the 99% expected shortfall and the market value margin.  The market value margin is calculated using a cost of capital method, where the capital is given by the 99% expected shortfall of `ΔRTK` and the cost of capital factor is assumed to be 6%.

 The example calculation provided here is based on an extremely simplified life insurance portfolio.

### Solvency 2 Life Calculation

The files *Life__Types.jl*, *Life_Constructors.jl* and *Life_Functions.jl* implement a simplified model of a life insurance company.  A simplified Solvency 2 calculation for this model is implemented in the files *S2Life__Types.jl*, *S2Life_Constructors.jl* und *S2Life_Functions.jl*. The Solvency 2 Standard Formula for Life insurers consists of several modules for the different risk factors. For each module a capital requirement is calculated first, and then all these requirements are aggregated to the overall Solvency Capital Requirement `SCR`. Most of the individual capital requirements are calculated through deterministic stress tests, which are calibrated to correspond to the 99.5% quantile.  Only a subset of the Solvency 2 modules applicable for life insurance companies have been implemented, but it should not be difficult to extend the code to more modules.


### Solvency 2 Non-Life Calculation

The files *S2NonLife__Types.jl* and *S2Life_Functions.jl* contain a simplified implementation of the Solvency 2 capital requirement for a highly simplified non-life insurance company.

### Internal Economic Capital Model

The internal economic capital model is a highly simplified Monte Carlo model of an extremely simplified non-life insurance company.  Reserves are completely ignored in the model. The purpose of this model is to illustrate some techniques used in value based management.

## Files in the Folder "test"

The folder "test" contains the files files *x_Input.jl*, *x.jl*, and *x_test.jl* for each of the four parts *x* ∈ {*SSTLife*, *S2Life*, *S2NonLife*, *ECModel*}. The files *x_Input.jl* and *x.jl* replicate the calculation in the book and the files *x_test.jl* are used for automatic testing.  In addition, the file *Life_Input.jl* contains the description of the example insurer, which is used  in the S2 life calculation.  *runtests.jl* controls the automated tests.

## Footnotes

<sup>1</sup> The abbreviation `ZK` comes from the original German term "Zielkapital". "Target capital" is the literal English translation.

<sup>2</sup> The abbreviation `RTK` comes from the  original German term "Risikotragendes Kapital". Observe that `RBC` is usually understood to mean "risk based capital" which has a different meaning than "risk bearing capital".  
