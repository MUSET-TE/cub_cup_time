# cub_cup_time
Project to develop  time/longitudinal versions of the CUB and CUP models, in simple and hierarchical variants

Panel CUB models for repeated ordinal ratings, with a **free shelter effect**:
an `R`/`shiny` application and a dependency-light estimation engine.

A rating is modelled as a mixture of three components, each with a logistic
predictor of its own:

* a **shifted binomial**, carrying evaluation (`phi`);
* a **discrete uniform**, carrying indecision (`upsilon`);
* a **point mass** at an analyst-declared category, carrying retreat to a
  refuge (`delta`).

Because the three predictors are separate, a covariate may raise indecision
while leaving evaluation untouched, and the software reports that as such.
A positive coefficient always means *more of the quantity named in its own
equation*; complements are never tabulated.

## Install and run

```r
install.packages("shiny")            # required
install.packages("readxl")           # optional: .xlsx upload
shiny::runApp(".")
```

## Licence

All code and data in this repository are licensed under the terms described in
[`LICENSE`](./LICENSE).

## Acknowledgement
This work was supported by a grant of the Ministry of Research, Innovation and Digitization, CNCS - UEFISCDI, project number PN-IV-P2-2.1-TE-2023-1317, within PNCDI IV. 

## Cite
Use is restricted to academic, non-commercial purposes only. Citation of the following article is required: TBA upon publication
