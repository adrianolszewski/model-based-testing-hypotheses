# Proving the equivalence between the 2-sample Wald’s z-statistic for comparing proportions with unpooled variances and the Average Marginal Effect over logistic regression with a single binary predictor
#### Adrian Olszewski, 02.03.2025

The Wald’s z-statistic for difference in 2 proportions with unpooled variances is of the following form:

``` math
\begin{equation}
z=\ \frac{\widehat{p_B}-\widehat{p_A}}{\sqrt{\frac{\widehat{p_A}\left(1-\widehat{p_A}\right)}{n_A}+\frac{\widehat{p_B}\left(1-\widehat{p_B}\right)}{n_B}}}
\hspace{2cm} (1)
\end{equation}
```
Where $\widehat{p_1}\$ stands for the estimated probability (sample proportion, %) in the 1st group, $\widehat{p_2}$ is the estimated probability (sample proportion, %) in the 2nd group, $n_1$ and $n_2$ denote respective group sizes.
I will show that this statistic is 1:1 equivalent to the average marginal effect (AME) of the logistic regression with a single binary predictor playing role of indicator for the compared samples.

To simplify calculations, I will show the equivalence of the nominator and denominator of the z statistic, i.e. the difference between two probabilities and its variance, separately.

**The equation of logistic regression**

Let’s start from the equation of the logistic regression with a single binary predictor:
``` math
\begin{equation}
logit\left(E\left(Y\middle| X\right)=logit\left(\hat{p}(Y=1|X\right)\right)=\ln{\left(\frac{\hat{p}\left(Y=1\middle| X\right)}{1-\hat{p}\left(Y=1\middle| X\right)}\right)}=\beta_0+X_1\beta_1
\hspace{2cm} (2)
\end{equation}
```
or equivalently, after applying the inverse-logit, i.e. sigmoid function (let’s also simplify X1 to X)
``` math
\begin{equation}
\hat{p}\left(Y=1\middle| X\right)=\frac{e^{\beta_0+X\beta_1}}{1+e^{\beta_0+X\beta_1}}
\hspace{2cm} (3)
\end{equation}
```
where $$Y_i$$​ are independent Bernoulli random variables with probabilities $$\hat{p_i}$$​.

**(Average) Marginal Effect**

Now, let’s introduce the marginal effect (ME) of a regression model. The ME of a given predictor variable X is the slope of the regression surface with respect to X, reflecting the rate at which Y changes at a given point with respect to X, while holding values of all other predictors constant. In other words, ME is an instantaneous rate of change, calculated as partial derivative of the regression equation with respect to X. For example, for an equation: Y = β0 + β1X1+ β2X2 + β3X1X2  corresponding ME(Y, X2) = ∂Y/ ∂X2 = β2 + β3X1
The average marginal effect (AME) calculates ME at every observed value of X and averages across the resulting effect estimates.  In other words, AME is an average effect of an infinitesimal change in X across all observations:

``` math
\begin{equation}
{AME}_{Xk}(model)=\beta_k\times\frac{1}{N}\sum_{i=1}^{N}\frac{\partial E(Y_i|{Xk}_i,\ covariates)}{\partial Xk}
\hspace{2cm} (4)
\end{equation}
```
 
For a generally defined equation:
``` math
\begin{equation}
\sigma\left(x\right)=\frac{e^x}{1+e^x}
\hspace{2cm} (5)
\end{equation}
```

the partial derivative (using the quotient rule and rewriting back in terms of σ) is of the form:
``` math
\begin{equation}
\frac{d\sigma(x)}{dx}=\sigma(x)(1-\sigma\left(x\right))
\hspace{2cm} (6)
\end{equation}
```

Therefore,
``` math
\begin{equation}
{AME}_X(model)=\beta\times\frac{1}{N}\sum_{i=1}^{N}{\hat{p}\left(Y=1\middle| X=x_i\right)\times\hat{p}\left(Y=0\middle| X=x_i\right)}
\hspace{2cm} (7)
\end{equation}
```

**AME for the binary predictor**

For a categorical predictor, however, there is no something like “infinitesimal change”. There is just switch between categories, so the AME becomes a contrast, i.e. a difference between values of the appropriate partial derivative calculated at the selected category levels, thus:

``` math
\begin{equation}
{AME}_X\left(model\right)=\hat{p}\left(Y=1\middle| X=B\right)-\hat{p}\left(Y=1\middle| X=A\right)=\ \widehat{p_B}-\widehat{p_A}
\hspace{2cm} (8)
\end{equation}
```

Which means, that the AME for such defined logistic regression corresponds to a difference in two estimated group probabilities, expressed in percentage points.

**Variance of AME for the binary predictor**

Now, I will show the equivalence between the variance in the Wald’s z statistic and the variance of the AME for the binary predictor.

``` math
\begin{equation}
var\left(AME\right)=var(\widehat{p_B}-\widehat{p_A})=\ \frac{\widehat{p_A}\left(1-\widehat{p_A}\right)}{n_A}+\frac{\widehat{p_B}\left(1-\widehat{p_B}\right)}{n_B}
\hspace{2cm} (9)
\end{equation}
```
Several replacements will save a lot of typing:
-	$$P_A=\ \hat{p}\left(Y=1\middle| X=A\right),\ P_B=\ \hat{p}\left(Y=1\middle| X=B\right)$$
-	$$P_i\times\left(1-P_i\right)=P_iQ_i$$

Let’s introduce the AME function:
``` math
\begin{equation}
AME=\ g\left(\beta_0,\beta_1\right)=P_B-P_A
\hspace{2cm} (10)
\end{equation}
```

Let’s also encode the two levels {A, B} using a single binary predictor X such that: A: X=0, B: X=1, let’s express $$P_A$$ and $$P_B$$ in terms of beta coefficients:
``` math
\begin{equation}
\left\{
\begin{aligned}
    P_A &= P(Y = 1, X = 0) = \frac{e^{\beta_0}}{1 + e^{\beta_0}} \\
    P_B &= P(Y = 1, X = 1) = \frac{e^{\beta_0 + \beta_1}}{1 + e^{\beta_0 + \beta_1}}
\end{aligned}
\right.
\hspace{2cm} (11)
\end{equation}
```
so the AME function is now expressed as:

``` math
\begin{equation}
g\left(\beta_0,\beta_1\right)=\ \frac{e^{\beta_0+\beta_1}}{1+e^{\beta_0+\beta_1}}-\frac{e^{\beta_0}}{1+e^{\beta_0}}
\hspace{2cm} (12)
\end{equation}
```

The variance for the AME is typically obtained by the _delta method_:
``` math
\begin{equation}
var\left(AME\right)=var\left(g\left(\beta_0,\beta_1\right)\right)\approx{\nabla g}^T\left(\beta_0,\beta_1\right)\times\Sigma\times\nabla g\left(\beta_0,\beta_1\right)
\hspace{2cm} (13)
\end{equation}
```

Let’s first obtain the derivatives. 
Recall, that:

``` math
\begin{equation}
\begin{aligned}
    \text{for } \sigma(x) &= \frac{e^x}{1 + e^x} \\
    \frac{d\sigma(x)}{dx} &= \sigma(x)(1 - \sigma(x)) = \frac{e^x}{1 + e^x} \times \frac{1}{1 + e^x} = \frac{e^x}{(1 + e^x)^2}
\end{aligned}
\hspace{2cm} (14)
\end{equation}
```
Therefore,

``` math
\begin{equation}
    \nabla g\left(\beta_0,\beta_1\right) =
    \left[
        \begin{matrix}
            \frac{\partial g}{\partial\beta_0} \\
            \frac{\partial g}{\partial\beta_1}
        \end{matrix}
    \right]
    =
    \left[
        \begin{matrix}
            P_B\left(1 - P_B\right) - P_A\left(1 - P_A\right) \\
            P_B\left(1 - P_B\right)
        \end{matrix}
    \right]
    =
    \left[
        \begin{matrix}
            P_B Q_B - P_A Q_A \\
            P_B Q_B
        \end{matrix}
    \right]
    =
    \left[
        \begin{matrix}
            \frac{e^{\beta_0+\beta_1}}{\left(1+e^{\beta_0+\beta_1}\right)^2} \\
            \frac{e^{\beta_0}}{\left(1+e^{\beta_0}\right)^2}
        \end{matrix}
    \right]
    \hspace{1.5cm} (15)
\end{equation}
```
Now, we need the variance-covariance matrix, i.e.
``` math
\begin{equation}
\Sigma=\left[\begin{matrix}var(\beta_0)&cov(\beta_0,\beta_1)\\cov(\beta_0,\beta_1&var(\beta_1)\\\end{matrix}\right]
\hspace{2cm} (16)
\end{equation}
```

This can be obtained by inverting the Fisher information matrix given by:
``` math
\begin{equation}
\Sigma=I^{-1}=\left(X^TWX\right)^{-1}
\hspace{2cm} (17)
\end{equation}
```

where X is the design matrix with 2 columns ($β_0$ of ones and $β_1$ indicating when X=1), with $n_A$ and $n_B$ number of rows corresponding to group A and B, respectively.

``` math
\begin{equation}
X=\left[\begin{matrix}1&A=0\\1&A=0\\\vdots&\vdots\\1&A=0\\1&B=1\\\vdots&\vdots\\1&B=1\\\end{matrix}\right]
\hspace{2cm} (18)
\end{equation}
```
and W is the diagonal matrix of weights, of the block-diagonal form:
``` math
\begin{equation}
W=diag\left(P_i\times\left(1-P_i\right)\right)=\left[\begin{matrix}P_AQ_A&0&\ldots&0&0&0&\ldots&0\\0&P_AQ_A&\cdots&0&0&0&\ldots&0\\\vdots&\vdots&\ddots&\vdots&\vdots&\vdots&\ddots&\vdots\\0&0&\ldots&P_AQ_A&0&0&\ldots&0\\0&0&\ldots&0&P_BQ_B&0&\ldots&0\\0&0&\ldots&0&0&P_BQ_B&\ldots&0\\\vdots&\vdots&\ddots&\vdots&\vdots&\vdots&\ddots&\vdots\\0&0&\ldots&0&0&0&\ldots&P_BQ_B\\\end{matrix}\right]
\hspace{2cm} (19)
\end{equation}
```
which can be simplified to:

``` math
\begin{equation}
W=diag\left(P_i\times\left(1-P_i\right)\right)=\left[\begin{matrix}P_AQ_AI_{nA}&0\\0&P_BQ_BI_{nB}\\\end{matrix}\right]
\hspace{2cm} (20)
\end{equation}
```
where $I_nA$ and $I_nB$ are respective identity matrices.

The multiplication result can be abbreviated to appropriate sums:
``` math
\begin{equation}
X^TWX=\left[\begin{matrix}\Sigma1P_iQ_i&\Sigma X P_iQ_i\\\Sigma X P_iQ_i&\Sigma X^2P_iQ_i\\\end{matrix}\right]
\hspace{2cm} (21)
\end{equation}
```

where **1** is the result of multiplying 1 x 1 (the $β_0$ vector), and **X**s refer to the other products of the $β_0$ and $β_1$ vectors. Since these vectors consist only of 0 and 1, and 0 refers to the “A” group and 1 refers to the “B” group, their products keep only the “B” part.
Therefore, the final result is:
``` math
\begin{equation}
X^TWX=\left[\begin{matrix}{n_AP}_AQ_A+{n_BP}_BQ_B&{n_BP}_BQ_B\\{n_BP}_BQ_B&{n_BP}_BQ_B\\\end{matrix}\right]
\hspace{2cm} (22)
\end{equation}
```

So the inverse can be computed (remembering that both covariance elements are equal):
``` math
\begin{equation}
\Sigma=\left(X^TWX\right)^{-1}=\frac{1}{\left({n_AP}_AQ_A+{n_BP}_BQ_B\right)\times{n_BP}_BQ_B-\left({n_BP}_BQ_B\right)^2}\left[\begin{matrix}{n_BP}_BQ_B&-{n_BP}_BQ_B\\-{n_BP}_BQ_B&{n_AP}_AQ_A+{n_BP}_BQ_B\\\end{matrix}\right]
\hspace{2cm} (23)
\end{equation}
```

After simplifying the denominator term:
``` math
\begin{equation}
\left({n_AP}_AQ_A+{n_BP}_BQ_B\right)\times{n_BP}_BQ_B-\left({n_BP}_BQ_B\right)^2={n_AP}_AQ_A\times{n_BP}_BQ_B+\left({n_BP}_BQ_B\right)^2-\left({n_BP}_BQ_B\right)^2={n_AP}_AQ_A\times{n_BP}_BQ_B 
\hspace{1cm} (24)
\end{equation}
```
finally:

``` math
\begin{equation}
\Sigma=\frac{1}{{n_AP}_AQ_A\times{n_BP}_BQ_B}\left[\begin{matrix}{n_BP}_BQ_B&-{n_BP}_BQ_B\\-{n_BP}_BQ_B&{n_AP}_AQ_A+{n_BP}_BQ_B\\\end{matrix}\right]=\left[\begin{matrix}\frac{1}{{n_AP}_AQ_A}&-\frac{1}{{n_AP}_AQ_A}\\-\frac{1}{{n_AP}_AQ_A}&\frac{{n_AP}_AQ_A+{n_BP}_BQ_B}{{n_AP}_AQ_A\times{n_BP}_BQ_B}\\\end{matrix}\right] 
\hspace{2cm} (25)
\end{equation}
```

By recalling formula #13 and matrix #15, we can express the variance of AME as:
``` math
\begin{equation}
var\left(AME\right)=\left(\frac{\partial g}{\partial\beta_0}\right)^2var\left(\beta_0\right)+\left(\frac{\partial g}{\partial\beta_1}\right)^2var\left(\beta_1\right)+2\ast\frac{\partial g}{\partial\beta_0}\frac{\partial g}{\partial\beta_1}covar\left(\beta_0,\ \beta_1\right)
\hspace{2cm} (26)
\end{equation}
```

which expands to:

``` math
\begin{equation}
var\left(AME\right)=\frac{\left(P_BQ_B-P_AQ_A\right)^2}{{n_AP}_AQ_A}+\frac{\left(P_BQ_B\right)^2\times\left({n_AP}_AQ_A+{n_BP}_BQ_B\right)}{{n_AP}_AQ_A\times{n_BP}_BQ_B}-2\ast\frac{\left(P_BQ_B-P_AQ_A\right)\times P_BQ_B}{{n_AP}_AQ_A}
\hspace{2cm} (27)
\end{equation}
```

Let's replace $P*Q$ with V to simplify (note: I occasionally use "×" to denote arithmetic multiplication to facilicate reading; A×B is just AB):

``` math
\begin{equation}
\begin{aligned}
    var(AME) &= \frac{(V_B - V_A)^2}{n_A V_A} + \frac{V_B^2 \times (n_A V_A + n_B V_B)}{n_A V_A \times n_B V_B} - 2 \times \frac{(V_B - V_A) \times V_B}{n_A V_A} \\[10pt]
    &= \frac{(V_B - V_A)(V_B - V_A - 2V_B)}{n_A V_A} + \frac{V_B \times (n_A V_A + n_B V_B)}{n_A V_A \times n_B} \\[10pt]
    &= \frac{(V_B - V_A)(-V_B - V_A)}{n_A V_A} + \frac{V_B n_A V_A + V_B^2 n_B}{n_A V_A \times n_B} \\[10pt]
    &= \frac{V_A^2 - V_B^2}{n_A V_A} + \frac{V_B}{n_B} + \frac{V_B^2}{n_A V_A} \\[10pt]
    &= \frac{V_A}{n_A} + \frac{V_B}{n_B}
    \hspace{1cm} (28)
\end{aligned}
\end{equation}
```
Therefore,
``` math
\begin{equation}
var\left(AME\right)=\frac{V_A}{n_A}+\frac{V_B}{n_B}=\frac{P_AQ_A}{n_A}+\frac{P_BQ_B}{n_B}=\frac{\widehat{p_A}\left(1-\widehat{p_A}\right)}{n_A}+\frac{\widehat{p_B}\left(1-\widehat{p_B}\right)}{n_B}\blacksquare 
\hspace{2cm} (29)
\end{equation}
```
This way I have shown the equivalence of the 2-sample Wald’s z-statistic for comparing proportions with unpooled variances and the average marginal effect of the logistic regression with a single binary predictor distinguishing the compared samples.
Also, since the Estimated Marginal Means (EM-Means) on so defined logistic regression and “re-grided” to the probability scale represent the estimated probabilities, the contrast comparing them through the Wald’s approach yields exactly the same result.
 
Although some minor discrepancies exist due to Maximum Likelihood Estimation, even for so small samples (N=10 and 20) the agreement is just perfect.

**Average Marginal Effect**
```r
> wald_z_test_non_pooled(x1 = 6, n1 = 20, x2 = 10, n2 = 20)
  diff         z    chi2        se  p.value  p.value_1        LCI        HCI
1 -0.2 -1.318761 1.73913 0.1516575 0.187249 0.09362452 -0.4972433 0.09724326
> 
> data <- data.frame(response = factor(c(rep("Success", 6), rep("Failure", 20-6),
+                                        rep("Success", 10), rep("Failure", 20-10))),
+                    grp    = factor(rep(c("B", "A"), each=20)))
> 
> m <- glm(response ~ grp, data = data, family = binomial(link = "logit"))
> data.frame(marginaleffects::avg_slopes(m)) %>% mutate(across(where(is.numeric), ~round(., 6)))
  term contrast estimate std.error statistic  p.value  s.value  conf.low conf.high
1  grp    B - A     -0.2  0.151657 -1.318762 0.187249 2.416973 -0.497243  0.097243
```
![obraz](https://github.com/user-attachments/assets/4f8db144-1a41-4baa-bf10-df8e062bf6ff)

**EM-means**
```r
> library(emmeans)
> update(pairs(emmeans(m, specs = ~grp, regrid="response")), infer = c(TRUE, TRUE)) %>% 
+     data.frame() %>% 
+     mutate(across(where(is.numeric), ~round(., 6)))
  contrast estimate       SE  df asymp.LCL asymp.UCL  z.ratio  p.value
1    A - B      0.2 0.151657 Inf -0.097243  0.497243 1.318761 0.187249
```

![obraz](https://github.com/user-attachments/assets/56c91bf1-f22f-49a3-96ba-38a0470710b7)

The implementation of the z statistic:
```r
wald_z_test_non_pooled <- function(x1, n1, x2, n2, conf.level=0.95) {
    p1 <- x1/n1
    p2 <- x2/n2
    
    se_p1 <- sqrt(p1 * (1 - p1) / n1)
    se_p2 <- sqrt(p2 * (1 - p2) / n2)
    
    se_diff <- sqrt(se_p1^2 + se_p2^2)
    
    z <- (p1 - p2) / se_diff
    p <- 2 * (1 - pnorm(abs(z)))
    hCI <- abs(qnorm((1 - conf.level)/2)) * se_diff

    return(data.frame(diff=p1-p2,
                      z = z, chi2 = z^2,
                      se = se_diff, 
                      p.value = p, p.value_1 =p/2,
                      LCI = (p1-p2) - hCI,
                      HCI = (p1-p2) + hCI,
                      row.names = NULL))
}
```

````````
