<!--
--- 
title: "R Package LSSPCA" 
---
author: "Giovanni Merola" 
output: html 
-->
# R Package LSSPCA
<!-- badges: start -->
[![R build status](https://github.com/merolagio/LSSPCA//workflows/R-CMD-check/badge.svg)](https://github.com/merolagio/LSSPCA//actions)
[![Project Status: WIP – Initial development is in progress, but there has not yet been a stable, usable release suitable for the public.](https://www.repostatus.org/badges/latest/wip.svg)](https://www.repostatus.org/#wip)
 [![license](https://img.shields.io/github/license/mashape/apistatus.svg)](https://choosealicense.com/licenses/mit/)
[![Last-changedate](https://img.shields.io/badge/last%20change-2021--04--10-yellowgreen.svg)](/commits/main)"
<!-- badges: end -->

  This R package is a companion to my tutorial paper published on xxx and computes sparse principal components. 
  
  LSSPCA provides only 3 functions 
  
-  *lsspca* is the basic function that computes the sparse PCs. 

-  *lsspca_blocked* computes lsspca separately on specified subsets of variables.

-  *makevexp* computes the variance explained by any set of components, also not computed with *lsspca*, of course.

It also provides 2 utilities (not methods, for various reasons)

- *print_spca* prints the sparse loadings

- *summary_spca* prints summary statistics for the SPCs produced with *lsspca* or *lsspca_blocked*

  This is a lightweight package not designed to handle large matrices. A dated (but working, made obsolete on CRAN)  package with methods for visualizing the PCs is available [here](https://github.com/merolagio/spca) (or *devtools::install_github("merolagio/spca")*). A new version with fast C++ code and PSPCA is in the making and will be released on CRAN one day. By the way, if you find this package useful, please acknowledge my work. It will make my manager happy :wink:

  The number of nonzero loadings can be controlled by changing the parameter *alpha* (which is the minimal proportion of variance explained by the Principal Components to be reproduced).
  
  Orthogonal sparse components (USPCA) are computed by choosing *spcaMethod = "u"*. Correlated components (CSPCA) can be obtained by choosing *spcaMethod = "c"*. The higher order CSPCA components may explain a bit more variance at the price of being correlated. 
  
  The variables can be selected by using exhaustive, stepwise, forward or backward selection via the argument *variableSelection* options "e", "s", "f" and "b", respectively.
  
  An example of sparse loadings is in this image
  
![](man/figures/readme_fig1.png)

  The sparse PCs are combinations of only 2, 3 and 4 variables out of 16 but are a (very) close approximation to the original PCs. See by yourself:

![](man/figures/readme_fig2.png)

  Well, the idea is to approximate the data as well as possible with sparse components. Since the PCs give the best approximation of the data, approximating the PCs with sparse components is pretty much the same thing. So, decent sparse components can be obtained by simply projecting (yes, by linear regression, PSPCA) the PCs on a subset of variables, option *spcaMethod = "p"*.
  
  Explanations, details and examples about LS SPCA can be found in the tutorial paper xxx.  
  
  More can be found in the mathematically oriented papers:

-  Merola, G. M. (2015). Least squares sparse principal component analysis: a backward elimination approach to attain large loadings. Australia & New Zealand Journal of Statistics, 57:391–429. [pdf](https://arxiv.org/abs/1406.1381)

-  Merola, G. M. and Chen, G. (2019). Projection sparse principal component analysis: An efficient least squares method. Journal of Multivariate Analysis, 173:366–382. [pdf](https://arxiv.org/abs/1612.00939)

  A discussion about sparsifying rotated principal components can be found in

-  Merola, G. M. (2020). Simpca: a framework for rotating and sparsifying principal components. Journal of Applied Statistics, 47(8):1325–1353. [pdf](https://arxiv.org/abs/1910.03266)


  Ahhh, I almost forgot: LS SPCA is much better than other SPCA methods because it maximizes the variance explained and produces orthogonal components which are combinations of key variables. Other methods do not and have serious problems with correlated variables. All explained and proved in my papers.
