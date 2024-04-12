# EGG: Estimating Genetic Gaussian Networks from GWAS Summary Data

The EGG (Estimating Genetic Gaussian Networks) package is designed to estimate multiple variable Genetic Gaussian Networks using GWAS summary data. 

It utilizes advanced statistical methods to infer the genetic relationships and dependencies between different traits or variables, providing a comprehensive view of the genetic architecture. 

This package is particularly useful for researchers and geneticists looking to understand the complex interplay between multiple genetic factors.

## Installation

You can install the development version of EGG from GitHub using the `devtools` package:

```R
if (!requireNamespace("devtools", quietly = TRUE))
    install.packages("devtools")
devtools::install_github("harryyiheyang/EGG")
```
Ensure that you have the data.table and glasso packages installed as EGG depends on these packages. They can be installed from CRAN using:
```
install.packages(c("data.table", "glasso"))
```
In particular, most matrix calculations in EGG package rely on CppMatrix Package, that convert the default function in R to C++ based functions:
```
devtools::install_github("harryyiheyang/CppMatrix")
```
## Example: Estimating the Genetic Network for 20 Traits in the EAS Population

This example demonstrates how to estimate the genetic network of 20 traits for the EAS population using publicly available GWAS summary data.

Step 1: Import Packages and Load Data
```R
library(igraph)
library(RColorBrewer)
library(EGG)
data("EAS_Error_COV")
data("EAS_Zscore")
BETA = EASZZ[,-1]
R = EASZR
head(BETA)
```
Step 2: Estimate the Genetic Network Using EGG
```R
fit1 = entropy.mcp.spearman.sampling(BETA, Rnoise = R, lamvec = c(25:52)/1125, max.eps = 0.005, max.iter = 25, rho = 0.05, mineig = 0.01, subfrac = 0.5, subthres = 0.5, subtime = 300, alpha = 0)
plot(rowMeans(fit1$cv.error), main = "CV Error of EGG")
```
Step 3: Visualize the Estimated Genetic Networks
```R
ThetaEUR = -fit1$Theta * (fit1$K > 0.95)
KEUR = fit1$K
nam = rownames(ThetaEUR)
colnames(ThetaEUR) = rownames(ThetaEUR) = nam
nam = sort(nam)
ThetaEUR = ThetaEUR[nam, nam]
gEUR <- graph_from_adjacency_matrix(ThetaEUR, weighted = TRUE, diag = FALSE, mode = "undirected")
V(gEUR)$size <- 12
E(gEUR)$color <- ifelse(E(gEUR)$weight > 0, "black", "grey80")
plot(gEUR, layout = layout.circle,
     vertex.size = V(gEUR)$size,
     vertex.label = V(gEUR)$name,
     vertex.label.cex = 2,
     vertex.label.color = "black",
     edge.color = E(gEUR)$color,
     vertex.frame.color = V(gEUR)$color,
     edge.width = 2, main = "Network estimate of EGG")
```
## License

This package is under the MIT License.

## Contact

For any questions or issues, please contact Yihe Yang at yxy1234@case.edu
