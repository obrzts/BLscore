
<!-- README.md is generated from README.Rmd. Please edit that file -->

# BLscore

<!-- badges: start -->
<!-- badges: end -->

BLscore provides functions to identify clusters of convergent TCRs that
are likely to recognize the same epitope. The clustering is based on
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
devtools::install_github("obrzts/BLscore")
```

## Usage

BLscore’s main function is `clusterize_TCR` which takes a data.frame
with TCR sequence data, calculate pairwise BL-scores, uses them to
define clusters of similar TCRs and return a data.frame with cluster
ids.

The input data.frame must contain following fields:

-   junction_beta - amino acid sequence of CDR3 plus the two flanking
    conserved residues;
-   v_beta, j_beta - V and J gene with or without allele; allele
    information is not used for score calculation.

If paired chain clustering is desired junction_alpha, v_alpha and
j_alpha must be provided too.

An example data.frame is provided with the package:

``` r
library(BLscore)
head(example_TCR_df)
#>     v_beta  j_beta   junction_beta    junction_alpha      v_alpha j_alpha id
#> 1   TRBV15 TRBJ1-2 CATRRNRGNTYGYTF    CAVRQTAAGNKLTF       TRAV21  TRAJ17  1
#> 2   TRBV13 TRBJ2-2    CASRQTSGELFF     CAVKGGGADGITF      TRAV8-1  TRAJ45  2
#> 3  TRBV7-9 TRBJ1-5 CASSSSLAGDQPQHF     CATDGGGAQKLVF       TRAV17  TRAJ54  3
#> 4 TRBV12-4 TRBJ2-1 CASSLSGGSYNEQFF    CVVNPVDSSYKLIF     TRAV12-1  TRAJ12  4
#> 5 TRBV12-4 TRBJ2-7 CASSSSGVGFYEQYF     CARGETSYDKVIF     TRAV13-1  TRAJ50  5
#> 6 TRBV12-3 TRBJ2-7    CASSFGVYEQYF CAYSSGAGGTSYGKLTF TRAV38-2/DV8  TRAJ52  6
```

Note, that the default clustering thresholds were defined to optimally
detect clusters of TCRs recognizing the same epitope. If instead of full
junction only CDR3 sequence witout flanking residues is provided the
scores will be overestimated which may lead to wrong cluster assignment.

Usage example:

``` r
clusters = clusterize_TCR(example_TCR_df, chains="AB", id_col = "id", tmp_folder=".", ncores=4)
head(clusters)
#>   cluster_id   v_beta  j_beta   junction_beta    junction_alpha      v_alpha
#> 1          1   TRBV15 TRBJ1-2 CATRRNRGNTYGYTF    CAVRQTAAGNKLTF       TRAV21
#> 2          2   TRBV13 TRBJ2-2    CASRQTSGELFF     CAVKGGGADGITF      TRAV8-1
#> 3          3  TRBV7-9 TRBJ1-5 CASSSSLAGDQPQHF     CATDGGGAQKLVF       TRAV17
#> 4          4 TRBV12-4 TRBJ2-1 CASSLSGGSYNEQFF    CVVNPVDSSYKLIF     TRAV12-1
#> 5          5 TRBV12-4 TRBJ2-7 CASSSSGVGFYEQYF     CARGETSYDKVIF     TRAV13-1
#> 6          6 TRBV12-3 TRBJ2-7    CASSFGVYEQYF CAYSSGAGGTSYGKLTF TRAV38-2/DV8
#>   j_alpha id
#> 1  TRAJ17  1
#> 2  TRAJ45  2
#> 3  TRAJ54  3
#> 4  TRAJ12  4
#> 5  TRAJ50  5
#> 6  TRAJ52  6
```
