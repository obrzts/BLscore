
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

BLscore’s main function is `clusterize_TCR()` which takes a data.frame
with TCR sequence data, calculate pairwise BL-scores, uses them to
define clusters of similar TCRs and returns a data.frame with cluster
ids.

The input data.frame must contain the following fields:

-   junction_beta - amino acid sequence of CDR3 plus the two conserved
    anchor residues;
-   v_beta, j_beta - V and J gene with or without allele; allele
    information is not used for the score calculation.

If paired chain clustering is desired junction_alpha, v_alpha and
j_alpha must be provided too.

Here is an example:

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

Note, that the default thresholds for clustering were defined for full
junction sequences. Providing only CDR3 sequence without anchor residues
may lead to less accurate cluster assignment. If you have only CDR3
sequences you can add anchor residues with `add_cdr3_anchors()`
function:

``` r
# make example table without anchor residues
df <- example_TCR_df
df$cdr3_beta <- substr(df$junction_beta, 2, nchar(df$junction_beta) - 2)

# generate column with beta chain junction sequence
df <- add_cdr3_anchors(df, chains="B", species="human")
head(df)
#>    j_beta  v_beta     junction_beta    junction_alpha  v_alpha j_alpha  id
#> 1 TRBJ1-1 TRBV4-1     CASSQGQGNTEAF CAEKRYAGGTSYGKLTF TRAV13-2  TRAJ52  45
#> 2 TRBJ1-1 TRBV5-1        CASVSGNEAF     CAVNGQAGTALIF  TRAV8-6  TRAJ15 786
#> 3 TRBJ1-1 TRBV6-3      CASRKTLNTEAF   CALSGLNSGYSTLTF   TRAV16  TRAJ11  96
#> 4 TRBJ1-1 TRBV5-4 CASSLSPGQGFRNTEAF   CAGQNRGIGANQFYF   TRAV25  TRAJ49  61
#> 5 TRBJ1-1  TRBV28      CASSFQGFTEAF       CADYYGQNFVF    TRAV3  TRAJ26 766
#> 6 TRBJ1-1 TRBV4-1   CASSQDNQQGGTEAF   CALGRYNNAGNMLTF   TRAV19  TRAJ39 813
#>         cdr3_beta
#> 1     ASSQGQGNTEA
#> 2        ASVSGNEA
#> 3      ASRKTLNTEA
#> 4 ASSLSPGQGFRNTEA
#> 5      ASSFQGFTEA
#> 6   ASSQDNQQGGTEA
```

Then run clustering:

``` r
clusters = clusterize_TCR(df, chains="AB", id_col = "id", tmp_folder=".", ncores=4)
head(clusters)
#>   cluster_id  j_beta   v_beta  junction_beta    junction_alpha      v_alpha
#> 1          1 TRBJ1-2   TRBV15 CATRRNRGNTYGYF    CAVRQTAAGNKLTF       TRAV21
#> 2          2 TRBJ2-2   TRBV13    CASRQTSGELF     CAVKGGGADGITF      TRAV8-1
#> 3          3 TRBJ1-5  TRBV7-9 CASSSSLAGDQPQF     CATDGGGAQKLVF       TRAV17
#> 4          4 TRBJ2-1 TRBV12-4 CASSLSGGSYNEQF    CVVNPVDSSYKLIF     TRAV12-1
#> 5          5 TRBJ2-7 TRBV12-4 CASSSSGVGFYEQF     CARGETSYDKVIF     TRAV13-1
#> 6          6 TRBJ2-7 TRBV12-3    CASSFGVYEQF CAYSSGAGGTSYGKLTF TRAV38-2/DV8
#>   j_alpha id
#> 1  TRAJ17  1
#> 2  TRAJ45  2
#> 3  TRAJ54  3
#> 4  TRAJ12  4
#> 5  TRAJ50  5
#> 6  TRAJ52  6
```
