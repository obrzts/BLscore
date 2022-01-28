
<!-- README.md is generated from README.Rmd. Please edit that file -->

# BLscore

<!-- badges: start -->
<!-- badges: end -->

BLscore provides functions to identify clusters of convergent TCRs that
are likely to recognize the same peptide. The clustering is based on
pairwise comparisons of all TCRs by aligning their CDR1-3 and
calculating the alignment scores using BLOSUM62 substitution matrix
reflecting evolutionary amino acid interchangeability. The three
alignment scores (CDR1, CDR2 and CDR3) are logistically transformed into
a single score (BL-score, BLOSUM-logistic) ranging from 0 to 1, where 1
corresponds to the highest probability of specificity match and 0 to the
lowest. Because the CDRs differ in their impact on TCR specificity, they
are weighted differently in the logistic function. The weights and
clustering thresholds were established using available data sets of TCRs
with known specificity (VDJdb and IEDB).

## Installation

You can install the development version of BLscore from
[GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("obrzts/BLscore_private")
```

## Usage

``` r
library(BLscore)
#head(sample_data_B)
```

``` r
#clusterize_TCR(sample_data_B, chains="AB", tmp_folder=".", ncores=4)
```
