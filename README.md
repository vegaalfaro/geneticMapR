
<!-- README.md is generated from README.Rmd. Please edit that file -->

# geneticMapper

<!-- badges: start -->
<!-- badges: end -->

<img src="man/figures/logo.png" align="right" width="140" />

`geneticMapper` simplifies and streamlines the construction of genetic
maps through an organized, reproducible workflow. Built on top of
[MapRtools](https://github.com/jendelman/MapRtools) and
[R/qtl](https://rqtl.org/), it makes genetic mapping more accessible
while supporting flexible analysis. Optimized for F2 diploid plant
populations, many functions were generalized to potentially other
species and different **ploidy levels**.

## Value

Genetic map construction involves several practical challenges.
`geneticMapper` was designed to help with reproducible genetic map
construction and quantitative trait loci (QTL) analysis.

## Installation

You can install the development version of `geneticMapper` from
[GitHub](https://github.com/) with:

``` r
# Get pak if needed
if (!requireNamespace("pak", quietly = TRUE)) {
  install.packages("pak")
}

#Install geneticMapper
pak::pak("vegaalfaro/geneticMapper")
#> ℹ Loading metadata database
#> ✔ Loading metadata database ... done
#> 
#> 
#> → Will install 66 packages.
#> → Will update 1 package.
#> → Will download 33 CRAN packages (33.04 MB), cached: 33 (58.14 MB).
#> → Will download 1 package with unknown size.
#> + abind                       1.4-8       ⬇ (65.11 kB)
#> + backports                   1.5.0      🔧 ⬇ (122.02 kB)
#> + broom                       1.0.8      
#> + car                         3.1-3      
#> + carData                     3.0-5       ⬇ (1.83 MB)
#> + cli                         3.6.4      🔧
#> + colorspace                  2.1-1      🔧 ⬇ (2.67 MB)
#> + corrplot                    0.95        ⬇ (3.83 MB)
#> + cowplot                     1.1.3      
#> + Deriv                       4.1.6       ⬇ (150.55 kB)
#> + distributional              0.5.0      
#> + doBy                        4.6.25     
#> + dplyr                       1.1.4      🔧
#> + fansi                       1.0.6      🔧 ⬇ (383.06 kB)
#> + farver                      2.1.2      🔧 ⬇ (1.97 MB)
#> + Formula                     1.2-5       ⬇ (158.56 kB)
#> + generics                    0.1.3       ⬇ (81.91 kB)
#> + geneticMapper  0.0.0.9000 → 0.0.0.9000 👷🏻‍♂️🔧 ⬇ (GitHub: caa16f9)
#> + ggdist                      3.3.2      🔧
#> + ggplot2                     3.5.1      
#> + ggpubr                      0.6.0      
#> + ggrepel                     0.9.6      🔧
#> + ggsci                       3.2.0      
#> + ggsignif                    0.6.4      
#> + glue                        1.8.0      🔧 ⬇ (173.70 kB)
#> + gridExtra                   2.3         ⬇ (1.11 MB)
#> + gtable                      0.3.6      
#> + HMM                         1.0.1      
#> + isoband                     0.2.7      🔧 ⬇ (1.87 MB)
#> + labeling                    0.4.3       ⬇ (61.49 kB)
#> + lifecycle                   1.0.4       ⬇ (124.78 kB)
#> + lme4                        1.1-37     🔧
#> + magrittr                    2.0.3      🔧 ⬇ (233.52 kB)
#> + MatrixModels                0.5-4      
#> + microbenchmark              1.5.0      🔧 ⬇ (72.58 kB)
#> + minqa                       1.2.8      🔧 ⬇ (340.28 kB)
#> + modelr                      0.1.11     
#> + munsell                     0.5.1       ⬇ (246.54 kB)
#> + nloptr                      2.2.1      🔧
#> + numDeriv                    2016.8-1.1  ⬇ (114.03 kB)
#> + pbkrtest                    0.5.3      
#> + pillar                      1.10.2     
#> + pkgconfig                   2.0.3       ⬇ (18.45 kB)
#> + polynom                     1.4-1       ⬇ (402.59 kB)
#> + purrr                       1.0.4      🔧
#> + qtl                         1.70       🔧 ⬇ (6.38 MB)
#> + quadprog                    1.5-8      🔧 ⬇ (40.38 kB)
#> + quantreg                    6.1        🔧
#> + R6                          2.6.1      
#> + rbibutils                   2.3        🔧 ⬇ (1.29 MB)
#> + RColorBrewer                1.1-3       ⬇ (53.32 kB)
#> + Rcpp                        1.0.14     🔧 ⬇ (3.36 MB)
#> + Rdpack                      2.6.3      
#> + reformulas                  0.4.0      
#> + rlang                       1.1.5      🔧 ⬇ (1.90 MB)
#> + rstatix                     0.7.2      
#> + scales                      1.3.0      
#> + SparseM                     1.84-2     🔧 ⬇ (942.60 kB)
#> + stringi                     1.8.7      🔧
#> + stringr                     1.5.1      
#> + tibble                      3.2.1      🔧
#> + tidyr                       1.3.1      🔧 ⬇ (1.32 MB)
#> + tidyselect                  1.2.1      
#> + utf8                        1.2.4      🔧 ⬇ (206.91 kB)
#> + vctrs                       0.6.5      🔧
#> + viridisLite                 0.4.2       ⬇ (1.30 MB)
#> + withr                       3.0.2       ⬇ (222.97 kB)
#> ℹ Getting 33 pkgs (33.04 MB) and 1 pkg with unknown size, 33 (58.14 MB) cached
#> ✔ Cached copy of geneticMapper 0.0.0.9000 (source) is the latest build
#> ✔ Cached copy of Deriv 4.1.6 (aarch64-apple-darwin20) is the latest build
#> ✔ Cached copy of Formula 1.2-5 (aarch64-apple-darwin20) is the latest build
#> ✔ Cached copy of RColorBrewer 1.1-3 (aarch64-apple-darwin20) is the latest build
#> ✔ Cached copy of Rcpp 1.0.14 (aarch64-apple-darwin20) is the latest build
#> ✔ Cached copy of SparseM 1.84-2 (aarch64-apple-darwin20) is the latest build
#> ✔ Cached copy of abind 1.4-8 (aarch64-apple-darwin20) is the latest build
#> ✔ Cached copy of backports 1.5.0 (aarch64-apple-darwin20) is the latest build
#> ✔ Cached copy of carData 3.0-5 (aarch64-apple-darwin20) is the latest build
#> ✔ Cached copy of colorspace 2.1-1 (aarch64-apple-darwin20) is the latest build
#> ✔ Cached copy of corrplot 0.95 (aarch64-apple-darwin20) is the latest build
#> ✔ Cached copy of fansi 1.0.6 (aarch64-apple-darwin20) is the latest build
#> ✔ Cached copy of farver 2.1.2 (aarch64-apple-darwin20) is the latest build
#> ✔ Cached copy of generics 0.1.3 (aarch64-apple-darwin20) is the latest build
#> ✔ Cached copy of glue 1.8.0 (aarch64-apple-darwin20) is the latest build
#> ✔ Cached copy of gridExtra 2.3 (aarch64-apple-darwin20) is the latest build
#> ✔ Cached copy of isoband 0.2.7 (aarch64-apple-darwin20) is the latest build
#> ✔ Cached copy of labeling 0.4.3 (aarch64-apple-darwin20) is the latest build
#> ✔ Cached copy of lifecycle 1.0.4 (aarch64-apple-darwin20) is the latest build
#> ✔ Cached copy of magrittr 2.0.3 (aarch64-apple-darwin20) is the latest build
#> ✔ Cached copy of microbenchmark 1.5.0 (aarch64-apple-darwin20) is the latest build
#> ✔ Cached copy of minqa 1.2.8 (aarch64-apple-darwin20) is the latest build
#> ✔ Cached copy of munsell 0.5.1 (aarch64-apple-darwin20) is the latest build
#> ✔ Cached copy of numDeriv 2016.8-1.1 (aarch64-apple-darwin20) is the latest build
#> ✔ Cached copy of pkgconfig 2.0.3 (aarch64-apple-darwin20) is the latest build
#> ✔ Cached copy of polynom 1.4-1 (aarch64-apple-darwin20) is the latest build
#> ✔ Cached copy of qtl 1.70 (aarch64-apple-darwin20) is the latest build
#> ✔ Cached copy of quadprog 1.5-8 (aarch64-apple-darwin20) is the latest build
#> ✔ Cached copy of rbibutils 2.3 (aarch64-apple-darwin20) is the latest build
#> ✔ Cached copy of rlang 1.1.5 (aarch64-apple-darwin20) is the latest build
#> ✔ Cached copy of tidyr 1.3.1 (aarch64-apple-darwin20) is the latest build
#> ✔ Cached copy of utf8 1.2.4 (aarch64-apple-darwin20) is the latest build
#> ✔ Cached copy of viridisLite 0.4.2 (aarch64-apple-darwin20) is the latest build
#> ✔ Cached copy of withr 3.0.2 (aarch64-apple-darwin20) is the latest build
#> ✔ Installed geneticMapper 0.0.0.9000 (github::vegaalfaro/geneticMapper@caa16f9) (96ms)
#> ✔ Installed Deriv 4.1.6  (118ms)
#> ✔ Installed Formula 1.2-5  (124ms)
#> ✔ Installed HMM 1.0.1  (125ms)
#> ✔ Installed MatrixModels 0.5-4  (127ms)
#> ✔ Installed R6 2.6.1  (129ms)
#> ✔ Installed RColorBrewer 1.1-3  (129ms)
#> ✔ Installed Rcpp 1.0.14  (147ms)
#> ✔ Installed Rdpack 2.6.3  (38ms)
#> ✔ Installed SparseM 1.84-2  (36ms)
#> ✔ Installed abind 1.4-8  (35ms)
#> ✔ Installed backports 1.5.0  (36ms)
#> ✔ Installed broom 1.0.8  (73ms)
#> ✔ Installed carData 3.0-5  (73ms)
#> ✔ Installed car 3.1-3  (38ms)
#> ✔ Installed cli 3.6.4  (60ms)
#> ✔ Installed colorspace 2.1-1  (42ms)
#> ✔ Installed corrplot 0.95  (39ms)
#> ✔ Installed cowplot 1.1.3  (37ms)
#> ✔ Installed distributional 0.5.0  (38ms)
#> ✔ Installed doBy 4.6.25  (38ms)
#> ✔ Installed dplyr 1.1.4  (39ms)
#> ✔ Installed fansi 1.0.6  (38ms)
#> ✔ Installed farver 2.1.2  (37ms)
#> ✔ Installed generics 0.1.3  (37ms)
#> ✔ Installed ggdist 3.3.2  (40ms)
#> ✔ Installed ggplot2 3.5.1  (44ms)
#> ✔ Installed ggpubr 0.6.0  (43ms)
#> ✔ Installed ggrepel 0.9.6  (89ms)
#> ✔ Installed ggsci 3.2.0  (73ms)
#> ✔ Installed ggsignif 0.6.4  (16ms)
#> ✔ Installed glue 1.8.0  (16ms)
#> ✔ Installed gridExtra 2.3  (17ms)
#> ✔ Installed gtable 0.3.6  (15ms)
#> ✔ Installed isoband 0.2.7  (19ms)
#> ✔ Installed labeling 0.4.3  (11ms)
#> ✔ Installed lifecycle 1.0.4  (16ms)
#> ✔ Installed lme4 1.1-37  (45ms)
#> ✔ Installed magrittr 2.0.3  (16ms)
#> ✔ Installed microbenchmark 1.5.0  (14ms)
#> ✔ Installed minqa 1.2.8  (14ms)
#> ✔ Installed modelr 0.1.11  (13ms)
#> ✔ Installed munsell 0.5.1  (14ms)
#> ✔ Installed nloptr 2.2.1  (21ms)
#> ✔ Installed numDeriv 2016.8-1.1  (13ms)
#> ✔ Installed pbkrtest 0.5.3  (14ms)
#> ✔ Installed pillar 1.10.2  (18ms)
#> ✔ Installed pkgconfig 2.0.3  (11ms)
#> ✔ Installed polynom 1.4-1  (14ms)
#> ✔ Installed purrr 1.0.4  (18ms)
#> ✔ Installed qtl 1.70  (50ms)
#> ✔ Installed quadprog 1.5-8  (12ms)
#> ✔ Installed quantreg 6.1  (21ms)
#> ✔ Installed rbibutils 2.3  (24ms)
#> ✔ Installed reformulas 0.4.0  (12ms)
#> ✔ Installed rlang 1.1.5  (42ms)
#> ✔ Installed rstatix 0.7.2  (15ms)
#> ✔ Installed scales 1.3.0  (17ms)
#> ✔ Installed stringi 1.8.7  (73ms)
#> ✔ Installed stringr 1.5.1  (17ms)
#> ✔ Installed tibble 3.2.1  (21ms)
#> ✔ Installed tidyr 1.3.1  (21ms)
#> ✔ Installed tidyselect 1.2.1  (14ms)
#> ✔ Installed utf8 1.2.4  (15ms)
#> ✔ Installed vctrs 0.6.5  (25ms)
#> ✔ Installed viridisLite 0.4.2  (13ms)
#> ✔ Installed withr 3.0.2  (15ms)
#> ✔ 1 pkg + 74 deps: kept 3, upd 1, added 66 [7.3s]

# Load library
library(geneticMapper)
```

## Example

This is a basic example which shows you how to solve a common problem:

``` r
# Load the example dataset
data("simulated_geno")

# Check markers previous to recoding
print(simulated_geno)
#>         Parent1 Parent2 F2_1 F2_2 F2_3
#> Marker1       0       2    0    1    2
#> Marker2       2       0    2    0    1
#> Marker3       0       2    1    2    1
#> Marker4       2       0    2    0    0
#> Marker5       0       2    0    2    1
#> Marker6       2       0    2    0    0

# Recode the markers using the recode() function
phased <- geneticMapper::recode(simulated_geno, parent1 = "Parent1", parent2 = "Parent2")
```

``` r
# Print the output
print(phased)
#>         Parent1 Parent2 F2_1 F2_2 F2_3
#> Marker1       0       2    0    1    2
#> Marker2       0       2    0    2    1
#> Marker3       0       2    1    2    1
#> Marker4       0       2    0    2    2
#> Marker5       0       2    0    2    1
#> Marker6       0       2    0    2    2
```

What is special about using `README.Rmd` instead of just `README.md`?
You can include R chunks like so
