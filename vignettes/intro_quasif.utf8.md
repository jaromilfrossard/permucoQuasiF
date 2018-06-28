---
title: "An introduction to the Quasi F statistics"
date: "2018-06-28"
author: "Jaromil Frossard"
output: rmarkdown::html_vignette
vignette: >
  %\VignetteIndexEntry{"An introduction to the Quasi F statistics"}
  %\VignetteEngine{knitr::rmarkdown}
  %\VignetteEncoding{UTF-8}
---

# The quasi-F statistic

Given overparametrized data of a fixed design $[X^0 ~D^0]$, participants indicator $Z_p^0$ matrix, items indicator $Z_i^0$ matrix, and the response $y$ variable, we are interested to test the effect of $X^0$. 

We code the matrices $X= X^0C_X$, $D= D^0C_D$ with ($\textrm{contr.sum}$) coding $C_X$ and $C_D$ such that $X$ and $D$ are full rank matrices.

$X$ is the matrices coding the main effect or interaction of factor that are either between participants, between items, or within both participants and subjects. Then the full rank matrix $Z_p = Z_p^0pC_p$ is computed using $C_p$, which is constructed using a ($\textrm{contr.sum}$) coding for each level of effect between participants. The same procedure is use to create $Z_i = Z_i^0C_i$. Then we define $Z_{pi} = \left( Z_{p}^\top * Z_i^\top\right)^\top (C_p\otimes C_i)\\$ where $*$ denotes the column wise Khatri-Rao product and $\otimes$ the Kronecker.

Then we compute the matrices:

$$
Z_{p:x} = \left( Z_p^\top *X_{p:within}^\top \right)^\top\\
Z_{i:x} = \left( Z_i^\top *X_{i:within}^\top \right)^\top\\
Z_{pi:x} = \left( Z_{pi}^\top *X_{pi:within}^\top \right)^\top\\
$$
 where $X_{p:within}$, $X_{i:within}$, $X_{pi:within}$ are matrices corresponding only to the column of the effect within participants  (respectively items, and both) of the design matrix $X$, bounded with a vector of one.
 
Using this construction and a balanced design (with respect to the experimental conditions, participants and items), the matrix $[X~D~Z_{p:x}~Z_{i:x}~Z_{pi:x}]$ is a square full rank matrix. In that setting, we write the saturated model:

$$
y = X\beta +D\eta + Z_{p:x}\gamma_p + Z_{i:x}\gamma_i + Z_{pi:x}\gamma_{pi} + \epsilon
$$
where $\gamma_p$, $\gamma_i$, $\gamma_{pi}$ and $\epsilon$ are random effect and we want to test the effect of $\beta$

Then, we write the quasi - F statistics:

$$qF = \frac{Num_1+Num_2}{Den_1+Den_2}\times \frac{df_{N}}{df_{D}}$$

where, 

$$Num_1 = y^\top H(X)y/ \textrm{rank}(X)\\
Num_2 = y^\top H(Z_{pi:x})y/ \textrm{rank}(Z_{pi:x}) \\
Den_1 = y^\top H(Z_{p:x})y/ \textrm{rank}(Z_{p:x}) \\
Den_2 = y^\top H(Z_{i:x})y/ \textrm{rank}(Z_{i:x}) \\$$

where $H(\cdot)$ is a fonction that compute the hat matrix and, \textrm{rank}(\cdot)$ is the rank of the matrix and,

$$
df_N = (Num_1+Num_2)^2+(Num_1^2 / \textrm{rank}(X) + Num_2^2 / \textrm{rank}(Z_{pi:x})\\
df_D = (Den_1+Den_2)^2+(Den_1^2 / \textrm{rank}(Z_{p:x}) + Den_1^2 / \textrm{rank}(Z_{i:x})
$$
Then the statistic follows approximatively a Fisher statistc under the null hypothesis: $qF \overset{approx}{\sim} F(df_N,df_D)$.
