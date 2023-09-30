This is a temporary file that I will separate file-per-test.
Now, only to not forget it.

-----
Let me show you how the logistic regression (with a few extensions) can be used to test hypotheses about fractions (%) of successes, repacling the classic "test for proportions".
Namely, it can replicate the results of:

1. [the Wald's (normal approximation) **z test for 2 proportions with non-pooled standard errors**](#wald_2prop_z) (common in clinical trials) via LS-means on the prediction scale or AME (average marginal effect)
2. [the Rao's score (normal appr.) **z test for 2 proportions with pooled standard errors**](#rao_2prop_z) (just what the `prop.test()` does in R)
3. the **z test for multiple (2+) proportions**
4. **ANOVA-like** (joint) test for multiple caterogical predictors (n-way ANOVA). Also (n-way) ANCOVA if you employ numerical covariates.
5. [the **Cochran-Mantel-Haenszel (CMH) for stratified/matched data**](#cmh) via _conditional logistic regression_
7. [the **Breslow-Day test for odds ratios**](#breslow-day) through Rao's ANOVA --> the interaction term
8. [the **Cochran-Armitage test for trend in ordered proportions**](#armitage-trend)
9. [the **McNemar and Cochran Q** test of paired proportions](#mcnemar) via GEE estimation (Generalized Estimating Equations with compound symmetry)
10. [the **Friedman test**](#mcnemar) - as above
11. [the **Mann-Whitney-Wilcoxon and Kruskal-Wallis**](#mww) via Ordinal Logistic Regression (and paired Wilcoxon via GEE)

Actually, the model-based approach to testing hypotheses is not anything new, and lots of other tests can be replicated with the general linear model via Ordinal Least Square (OLS) and Generalized Least Square (GLS) estimation, generalized linear models (GLM) via both Maximum-Likelhiood estimation (MLE) and semi-parametric Generalized Estimating Equations (GEE). Let's add to this also the conditional approach via Mixed-Effect models (both general and generalized). And let's not forget about the Quantile Regression (with mixed effects), robust regression models, survival models (Cox, AFT, Andersen-Gill, frailty models) and dozens of others!

All those models, followed by the Likelihood Ratio testing (LRT) or Wald's testing of model coefficients, especially combined with LS-means (EM-means) will give you incredibly flexible testing framework.

This time we will look at the Logistic Regression, part of the Generalized Linear Model - the binomial regression with logit link. We will also employ certain extensions i generalizations to achieve concrete effects.

---

We will use 3 data sets (defined at the bottom of this file):
* unpaired 2-group data

``` r
> head(unpaired_data)
     sex response    trt
1 female        0 active
2 female        0 active
3 female        0 active
4 female        0 active
5 female        0 active
6 female        0 active
> tail(unpaired_data)
     sex response     trt
101 male        1 placebo
102 male        1 placebo
103 male        1 placebo
104 male        1 placebo
105 male        1 placebo
106 male        1 placebo
```

* paired 2-group data

``` r
> head(paired_data)
  ID Time Treatment Response
1  1  Pre   placebo        0
2  1 Post   placebo        1
3  2  Pre   placebo        0
4  2 Post   placebo        0
5  3  Pre   placebo        0
6  3 Post   placebo        0
> 
> tail(paired_data)
   ID Time Treatment Response
35 18  Pre    active        0
36 18 Post    active        1
37 19  Pre    active        0
38 19 Post    active        0
39 20  Pre    active        0
40 20 Post    active        0
```

* ordered data
``` r
> head(ordered_paired_data)
  ID Time Response TimeUnord
1  1   T1        0        T1
2  2   T1        0        T1
3  3   T1        0        T1
4  4   T1        0        T1
5  5   T1        0        T1
6  6   T1        0        T1
> tail(ordered_paired_data)
   ID Time Response TimeUnord
25  5   T3        1        T3
26  6   T3        1        T3
27  7   T3        1        T3
28  8   T3        0        T3
29  9   T3        1        T3
30 10   T3        1        T3
```

* unpaired 2-group ordinal data (Pain score of the ODI (Oswestry Disability Index) questionnaire; 6-items Likert data.
https://www.lni.wa.gov/forms-publications/F252-130-000.pdf
``` r
> head(ordinal_data)
                 ODIPain Arm Age_centered
1      [2] Moderate pain   B     -6.15315
2            [0] No pain   B     12.84685
3     [1] Very mild pain   A     -9.15315
4      [2] Moderate pain   B     14.84685
5 [3] Fairly severe pain   A     12.84685
6      [2] Moderate pain   B      2.84685
> tail(ordinal_data)
                  ODIPain Arm Age_centered
106    [2] Moderate pain   A   -15.153153
107    [2] Moderate pain   B   -11.153153
108    [2] Moderate pain   A    -4.153153
109 [4] Very severe pain   B    -0.153153
110   [1] Very mild pain   B    -4.153153
111   [1] Very mild pain   B    -7.153153
```

---
Loading necessary packages
```{r}
library(emmeans)
library(broom)
library(survival)
library(marginaleffects)
library(geepack)
```

Defining auxiliary function (to validate the results)
``` r
wald_z_test <- function(table) {
  p1 <- prop.table(table, 1)[1, 1]
  p2 <- prop.table(table, 1)[2, 1]
  n1 <- rowSums(table)[1]
  n2 <- rowSums(table)[2]
  se_p1 <- sqrt(p1 * (1 - p1) / n1)
  se_p2 <- sqrt(p2 * (1 - p2) / n2)
  se_diff <- sqrt(se_p1^2 + se_p2^2)
  z <- (p1 - p2) / se_diff
  p <- 2 * (1 - pnorm(abs(z)))
  return(data.frame(estimate = p1 - p2, z = z, se = se_diff, p.value = p, row.names = NULL))
}
```

---
<a name="wald_2prop_z"></a>
# Wald's z test for 2 proportions (non-pooled SE)

We want to reproduce this result:
``` r
> wald_z_test(xtabs(~ trt + response,data = unpaired_data))
   estimate        z         se     p.value
1 0.2737968 3.047457 0.08984435 0.002307865
```

We will use this logistic regression (LR) model:
``` r
> summary(lr_model <- glm(response ~ trt , data = unpaired_data, family = binomial(link = "logit")))

Call:
glm(formula = response ~ trt, family = binomial(link = "logit"), 
    data = unpaired_data)

Deviance Residuals: 
    Min       1Q   Median       3Q      Max  
-1.7011  -1.1620   0.7325   1.0778   1.1929  

Coefficients:
            Estimate Std. Error z value Pr(>|z|)   
(Intercept) -0.03637    0.26972  -0.135  0.89274   
trtplacebo   1.21502    0.42629   2.850  0.00437 **
---
Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

(Dispersion parameter for binomial family taken to be 1)

    Null deviance: 140.50  on 105  degrees of freedom
Residual deviance: 131.88  on 104  degrees of freedom
AIC: 135.88

Number of Fisher Scoring iterations: 4
```

## Wald's z test via LS-means on re-grided scale (probability scale)
``` r
> pairs(emmeans(lr_model, regrid="response", specs = ~ trt))
 contrast         estimate     SE  df z.ratio p.value
 active - placebo   -0.274 0.0898 Inf  -3.047  0.0023
 ```
Let's look closer at the results:
| Outcome   | LS-means | raw z test | comment |
|-----------|----------|------------|---------|
| estimate  | -0.2737968 | 0.2737968| 👍; swap factor levels to change the sign or ignore |
| SE        | 0.08984432 | 0.08984435 | agreement by 7 dec. digits 👍 |
| statistic | -3.047458 | 3.047457 | sign - as above; agreement by 5 dec. digits 👍 |
| p-value   | 0.002307857 | 0.002307865 | aggrement by 7 dec. digits 👍 |

Excellent agreement!

## Wald's z test via AME (average marginal effect)
``` r
> marginaleffects::avg_slopes(lr_model)

 Term         Contrast Estimate Std. Error    z Pr(>|z|)   S  2.5 % 97.5 %
  trt placebo - active    0.274     0.0898 3.05  0.00231 8.8 0.0977   0.45

Columns: term, contrast, estimate, std.error, statistic, p.value, s.value, conf.low, conf.high 
 ```
Let's look closer at the results:
| Outcome   | AME | raw z test | comment |
|-----------|----------|------------|------|
| estimate  | 0.2737968 | 0.2737968 | | 👍 |
| SE        | 0.08984433 | 0.08984435 | 👍 |
| statistic | 3.047458 | 3.047457 | agreement by 5 dec. digits 👍 |
| p-value   | 0.002307859 | 0.002307865 | agreement by 6 dec. digits 👍 |

Perfect agreement!

---
<a name="rao_2prop_z"></a>
# Rao score z test for 2 proportions (pooled SE)

We want to reproduce this result:
``` r
> prop.test(xtabs(~ trt + response,data=unpaired_data), correct = FALSE)

	2-sample test for equality of proportions without continuity correction

data:  xtabs(~trt + response, data = unpaired_data)
X-squared = 8.4429, df = 1, p-value = 0.003665
alternative hypothesis: two.sided
95 percent confidence interval:
 0.09770511 0.44988848
sample estimates:
   prop 1    prop 2 
0.5090909 0.2352941 
```

We will use the same logistic regression (LR) model as previously

## Rao score z test via ANOVA with Rao test
``` r
> anova(glm(response ~ trt , data = unpaired_data, family = binomial(link = "logit")), test = "Rao")
Analysis of Deviance Table

Model: binomial, link: logit

Response: response

Terms added sequentially (first to last)

     Df Deviance Resid. Df Resid. Dev    Rao Pr(>Chi)   
NULL                   105     140.50                   
trt   1   8.6257       104     131.88 8.4429 0.003665 **
---
Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
```
Let's look closer at the results:
| Outcome   | ANOVA + Rao test | prop.test() | comment |
|-----------|----------|------------|---------|
| statistic | 8.442898 | 8.442897 | agreement by 5 dec. digits 👍 |
| p-value   | 0.003664718 | 0.003664719 | agreement by 8 dec. digits 👍 |

Perfect agreement!

---
<a name="breslow-day"></a>
# Breslow-Day test for odds ratios via ANOVA with Rao test

We want to reproduce this result for treatment and sex:
``` r
> BreslowDayTest(xtabs(~ trt +response + sex, data=unpaired_data), correct = TRUE)

	Breslow-Day Test on Homogeneity of Odds Ratios (with Tarone correction)

data:  xtabs(~trt + response + sex, data = unpaired_data)
X-squared = 1.4905, df = 1, p-value = 0.2221
```
This time add sex to the model and will look at the interaction term

``` r
> as.data.frame(anova(glm(response ~ trt * sex , data = unpaired_data, family = binomial(link = "logit")), test="Rao")[4, ])
        Df Deviance Resid. Df Resid. Dev      Rao  Pr(>Chi)
trt:sex  1 1.498573       102   130.0512 1.496552 0.2212027
```
Let's look closer at the results:
| Outcome   | ANOVA + Rao test | Breslow-Day | comment |
|-----------|----------|------------|---------|
| statistic | 1.496552 |  1.490537  | agreement by 2 dec. digits 👍 |
| p-value   | 0.2212027 | 0.2221331 | agreement bt 2 dec. digits 👍 |

Good agreement!

---
<a name="cmh"></a>
# (Cochrane-) Mantel-Haenszel via conditional logistic regression

We want to reproduce this result for sex strata:
``` r
> mantelhaen.test(unpaired_data$response, unpaired_data$trt, unpaired_data$sex, exact = F, correct = F)

	Mantel-Haenszel chi-squared test without continuity correction

data:  unpaired_data$response and unpaired_data$trt and unpaired_data$sex
Mantel-Haenszel X-squared = 8.3052, df = 1, p-value = 0.003953
alternative hypothesis: true common odds ratio is not equal to 1
95 percent confidence interval:
 1.445613 7.593375
sample estimates:
common odds ratio 
         3.313168 
```
And through the model:
``` r
> summary(clogit(response~trt + strata(sex),data=unpaired_data))$sctest
      test         df     pvalue 
8.30516934 1.00000000 0.00395324 
```
Let's look closer at the results:
| Outcome   | Cond. LR | CMH | comment |
|-----------|----------|------------|---------|
| statistic | 8.30516934 | 8.305169 | 👍 |
| p-value   | 0.00395324 | 0.00395324 | 👍 |

Ideal agreement!

---
<a name="mcnemar"></a>
# McNemar's, Cochran Q, Friedman tests via GEE estimated LR
We want to reproduce this result for sex strata:
``` r
> mcnemar.test(x=paired_data[paired_data$Time == "Pre", "Response"], y=paired_data[paired_data$Time == "Post", "Response"], correct = F)

	McNemar's Chi-squared test

data:  paired_data[paired_data$Time == "Pre", "Response"] and paired_data[paired_data$Time == "Post", "Response"]
McNemar's chi-squared = 10.286, df = 1, p-value = 0.001341

# or this one

> paired_data %>% rstatix::friedman_test(Response ~ Time |ID)
# A tibble: 1 × 6
  .y.          n statistic    df       p method       
* <chr>    <int>     <dbl> <dbl>   <dbl> <chr>        
1 Response    20      10.3     1 0.00134 Friedman test

# or this one

> RVAideMemoire::cochran.qtest(Response ~ Time | ID,data=paired_data)

	Cochran's Q test

data:  Response by Time, block = ID 
Q = 10.2857, df = 1, p-value = 0.001341
alternative hypothesis: true difference in probabilities is not equal to 0 
sample estimates:
proba in group Post  proba in group Pre 
                0.7                 0.1 
```

Through the GEE-estimated model:
``` r
> summary(geepack::geeglm(Response ~ Time, id = ID,data=paired_data, family = binomial(), corstr = "exchangeable"))

Call:
geepack::geeglm(formula = Response ~ Time, family = binomial(), 
    data = paired_data, id = ID, corstr = "exchangeable")

 Coefficients:
            Estimate Std.err   Wald Pr(>|W|)   
(Intercept)   0.8473  0.4880  3.015  0.08249 . 
TimePre      -3.0445  0.9484 10.305  0.00133 **
---
Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

Correlation structure = exchangeable 
Estimated Scale Parameters:

            Estimate Std.err
(Intercept)        1  0.7215
  Link = identity 

Estimated Correlation Parameters:
      Estimate Std.err
alpha  -0.1455  0.2819
Number of clusters:   20  Maximum cluster size: 2 

# or in a more compact form:
> coef(summary(geepack::geeglm(Response ~ Time, id = ID,data=paired_data, family = binomial(), corstr = "exchangeable")))[2,]
        Estimate Std.err  Wald Pr(>|W|)
TimePre   -3.045  0.9484 10.31 0.001327
```

Let's look closer at the results:
| Outcome   | GEE LR | Tests | comment |
|-----------|----------|------------|---------|
| statistic | 10.31 | 10.2857 | agreement by 1 deci. digits 👍 |
| p-value   | 0.001327 | 0.001341 | agreement by 4 dec. digits 👍 |
 
Acceptable agreement!

---
<a name="armitage-trend"></a>
# Cochrane-Armitage test for trend via GLM + ANOVA LRT (Likelihood Ratio Test)
We want to reproduce this result for sex strata:
``` r
> DescTools::CochranArmitageTest(xtabs(~Response + Time,data=ordered_paired_data))

	Cochran-Armitage test for trend

data:  xtabs(~Response + Time, data = ordered_paired_data)
Z = -3.6, dim = 3, p-value = 0.0003
alternative hypothesis: two.sided

# or this one

> rstatix::prop_trend_test(xtabs(~Response + Time,data=ordered_paired_data))
# A tibble: 1 × 6
      n statistic        p p.signif    df method               
* <int>     <dbl>    <dbl> <chr>    <dbl> <chr>                
1    30      12.9 0.000336 ***          1 Chi-square trend test
```

Through the GLM model:
``` r
> as.data.frame(anova(glm(Response ~ Time,data=ordered_paired_data, family = binomial()), test="LRT"))[2,]
     Df Deviance Resid. Df Resid. Dev  Pr(>Chi)
Time  2    14.99        27      26.46 0.0005553
```

Let's look closer at the results:
| Outcome   | GLM + LRT ANOVA | Test | comment |
|-----------|----------|------------|---------|
| statistic | 12.86 | 14.99 | same order of magnitude |
| p-value   | 0.0005553 | 0.000336 | agreement by 3 dec. digits 👍 |
 
Reasonable agreement. (Maybe I'll find a better one).

---
<a name="mww"></a>
# Mann-Whitney (-Wilcoxon) test of stochastic equivalence (vs. stochastic superiority / dominance)
**Note:** This test DOES NOT TEST MEDIANS in general, unless strong distributional assumptions hold:
1) IID samples (same dispersion, variance & same shape - if skewed, then in the same direction)
2) Symmetric around their medians.
For detailed explanations, read my gist and find a rich list of literature (mostly freely accessible) and examples: https://gist.github.com/adrianolszewski/2cec75678e1183e4703589bfd22fa8b2

We want to reproduce this result:
``` r
> (wtest <- wilcox.test(as.numeric(ODIPain) ~ Arm, data = ordinal_data, exact = FALSE, correct = FALSE))

	Wilcoxon rank sum test

data:  as.numeric(ODIPain) by Arm
W = 1472, p-value = 0.68
alternative hypothesis: true location shift is not equal to 0

> wtest$p.value
[1] 0.679575
```
By using the proportional-odds model (ordinal logistic regression) we obtain:
``` r
> coef(summary(m <- MASS::polr(ODIPain ~ Arm , data = ordinal_data, Hess=T)))
                                                   Value Std. Error   t value
ArmB                                            0.141709   0.341471  0.414995
[0] No pain|[1] Very mild pain                 -1.444439   0.299213 -4.827458
[1] Very mild pain|[2] Moderate pain           -0.273260   0.259784 -1.051875
[2] Moderate pain|[3] Fairly severe pain        1.361363   0.291704  4.666935
[3] Fairly severe pain|[4] Very severe pain     2.093502   0.345203  6.064551
[4] Very severe pain|[5] Worst imaginable pain  4.072209   0.736078  5.532306

> pairs(emmeans(m, specs = ~Arm))
 contrast estimate    SE  df z.ratio p.value
 A - B      -0.142 0.341 Inf  -0.415  0.6781

# or
> (mtest <- joint_tests(m))
 model term df1 df2 F.ratio p.value
 Arm          1 Inf   0.172  0.6781

mtest$p.value
[1] 0.678146
```

This time, the two outputs (model vs. test) look very different, but give a very close p-value!
It's not a coincidence. 
You can find detailed explanations and necessary formulas here: [Equivalence of Wilcoxon Statistic and Proportional Odds Model](https://www.fharrell.com/post/powilcoxon/) | [Resources for Ordinal Regression Models](https://www.fharrell.com/post/rpo/) | [If You Like the Wilcoxon Test You Must Like the Proportional Odds Model](https://www.fharrell.com/post/wpo/)

So, like Prof. Harrell, we will check also the concordance index:
``` r
# From the Wilcoxon statistic
> (bind_cols(tidy(wilcox.test(as.numeric(ODIPain) ~ Arm, data = ordinal_data, exact = FALSE, correct = FALSE)),
          ordinal_data %>% 
            group_by(Arm) %>% 
            summarize(n=n()) %>% 
            summarize("n1*n2" = prod(n))) %>% 
  mutate(c = statistic / `n1*n2`) -> concor)

# A tibble: 1 × 6
  statistic p.value method                 alternative `n1*n2`     c
      <dbl>   <dbl> <chr>                  <chr>         <dbl> <dbl>
1     1472.   0.680 Wilcoxon rank sum test two.sided      3080 0.478

> concor$c
0.478084

# From the odds ratio taken from the model:
> (OR <- 1/exp((coef(summary(m)))[1,1]))
[1] 0.867874

and finally
> (c_mod <- OR^0.66 / (1 + OR ^ 0.66))
0.476635

# So we are off by:
> sprintf("%.2f%%",100*(concor$c - c_mod) / concor$c)
[1] "0.30%"

# Isn't this IMPRESSIVE?
```

Let's collect the results closer at the results:
| Outcome   | OLR | Wilcox | comment |
|-----------|----------|------------|---------|
| concordance | 0.478084 | 0.476635 | agreement by 2 dec. digits 👍 |
| p-value   | 0.679575 | 0.678146 | agreement by 2 dec. digits 👍 |

Very good agreement!

Later, in a separate gist, I will shoud you, through simulation, that this equivalence holds very well!

**Think about the the consequences. This way obtain the Mann-Whitney (-Wilcoxon) test adjusted for covariates.**
By the way, this is another interesting example, where the result of a completely non-parametric test can be obtained via parametric method.

---
* _EM-means_ (estimated marginal means) is another name of the well-known in experimental research _LS-means_ (least-square means)
It's a model-based predicted (estimated) mean. If you remember the definition of regression (NO, not the Machine Learning one...)
then you know that regresion gives you a some function of the data conditional to the predictor.
For the linear regression it's E(Y|X=x), for the GLM it is link(E(Y|X=x)), for quantile regression it's median(Y|X=x).
And since the predictor exclusively consists of categorical variables, they form sub-groups in which the (conditional) 
means are calculated. If we include also numerical covariates into the model, the predictions will account for it, giving us so-called "covariate-adjusted means".

----
The datasets for your own experiments:

``` r
unpaired_data <- structure(list(sex = structure(c(1L, 1L, 1L, 1L, 1L, 1L, 1L, 
1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 
1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 
1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 2L, 2L, 2L, 
2L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 
2L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 
2L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 
2L, 2L, 2L), levels = c("female", "male"), class = "factor"), 
response = c(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 
1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 
0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 
1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 
1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1), trt = structure(c(1L, 
1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 
2L, 2L, 2L, 2L, 2L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 
1L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 
2L, 2L, 2L, 2L, 2L, 2L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 
1L, 1L, 1L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 1L, 1L, 1L, 1L, 1L, 
1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 2L, 2L, 2L, 2L, 
2L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 2L
), levels = c("active", "placebo"), class = "factor")), row.names = c(NA, 
-106L), class = "data.frame")

paired_data <- structure(list(ID = c(1L, 1L, 2L, 2L, 3L, 3L, 4L, 4L, 5L, 5L, 
6L, 6L, 7L, 7L, 8L, 8L, 9L, 9L, 10L, 10L, 11L, 11L, 12L, 12L, 
13L, 13L, 14L, 14L, 15L, 15L, 16L, 16L, 17L, 17L, 18L, 18L, 19L, 
19L, 20L, 20L), Time = structure(c(2L, 1L, 2L, 1L, 2L, 1L, 2L, 
1L, 2L, 1L, 2L, 1L, 2L, 1L, 2L, 1L, 2L, 1L, 2L, 1L, 2L, 1L, 2L, 
1L, 2L, 1L, 2L, 1L, 2L, 1L, 2L, 1L, 2L, 1L, 2L, 1L, 2L, 1L, 2L, 
1L), levels = c("Post", "Pre"), class = "factor"), Treatment = structure(c(2L, 
2L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 
2L, 2L, 2L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 
1L, 1L, 1L, 1L, 1L, 1L, 1L), levels = c("active", "placebo"), class = "factor"), 
Response = c(0L, 1L, 0L, 0L, 0L, 0L, 0L, 1L, 0L, 1L, 1L, 
0L, 0L, 1L, 0L, 0L, 1L, 1L, 0L, 1L, 0L, 1L, 0L, 1L, 0L, 1L, 
0L, 1L, 0L, 1L, 0L, 1L, 0L, 1L, 0L, 1L, 0L, 0L, 0L, 0L)), row.names = c(NA, 
-40L), class = "data.frame")

ordered_paired_data <- structure(list(ID = structure(c(1L, 2L, 3L, 4L, 5L, 6L, 7L, 8L, 
9L, 10L, 1L, 2L, 3L, 4L, 5L, 6L, 7L, 8L, 9L, 10L, 1L, 2L, 3L, 
4L, 5L, 6L, 7L, 8L, 9L, 10L), levels = c("1", "2", "3", "4", 
"5", "6", "7", "8", "9", "10"), class = "factor"), Time = structure(c(1L, 
1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 
2L, 2L, 2L, 3L, 3L, 3L, 3L, 3L, 3L, 3L, 3L, 3L, 3L), levels = c("T1", 
"T2", "T3"), class = c("ordered", "factor")), Response = c(0L, 
0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 1L, 1L, 1L, 1L, 1L, 1L, 0L, 0L, 
0L, 0L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 0L, 1L, 1L), TimeUnord = structure(c(1L, 
1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 
2L, 2L, 2L, 3L, 3L, 3L, 3L, 3L, 3L, 3L, 3L, 3L, 3L), levels = c("T1", 
"T2", "T3"), class = "factor")), row.names = c(NA, -30L), class = "data.frame")

ordinal_data <- structure(list(ODIPain = structure(c(3L, 1L, 2L, 3L, 4L, 3L, 
4L, 2L, 3L, 5L, 2L, 5L, 5L, 6L, 2L, 3L, 1L, 2L, 3L, 3L, 1L, 3L, 
3L, 2L, 2L, 5L, 5L, 2L, 5L, 3L, 5L, 1L, 3L, 3L, 3L, 1L, 5L, 3L, 
5L, 1L, 1L, 2L, 1L, 2L, 3L, 2L, 3L, 1L, 2L, 1L, 2L, 4L, 6L, 4L, 
3L, 3L, 3L, 3L, 1L, 4L, 5L, 4L, 3L, 3L, 1L, 3L, 1L, 4L, 3L, 3L, 
2L, 3L, 3L, 3L, 3L, 3L, 2L, 2L, 1L, 2L, 2L, 1L, 3L, 4L, 4L, 3L, 
2L, 2L, 2L, 2L, 2L, 1L, 1L, 3L, 1L, 3L, 1L, 3L, 4L, 4L, 3L, 3L, 
1L, 2L, 3L, 3L, 3L, 3L, 5L, 2L, 2L), levels = c("[0] No pain", 
"[1] Very mild pain", "[2] Moderate pain", "[3] Fairly severe pain", 
"[4] Very severe pain", "[5] Worst imaginable pain"), class = c("ordered", 
"factor")), Arm = structure(c(2L, 2L, 1L, 2L, 1L, 2L, 1L, 2L, 
2L, 2L, 1L, 1L, 1L, 2L, 2L, 1L, 1L, 2L, 1L, 2L, 2L, 1L, 2L, 1L, 
1L, 2L, 1L, 1L, 2L, 2L, 2L, 1L, 1L, 1L, 2L, 1L, 1L, 2L, 1L, 2L, 
2L, 2L, 1L, 1L, 2L, 1L, 1L, 2L, 1L, 2L, 1L, 2L, 1L, 2L, 1L, 2L, 
1L, 2L, 2L, 2L, 2L, 1L, 2L, 2L, 1L, 1L, 1L, 2L, 2L, 2L, 1L, 1L, 
1L, 1L, 1L, 2L, 2L, 2L, 1L, 2L, 2L, 1L, 2L, 1L, 1L, 2L, 1L, 2L, 
1L, 1L, 2L, 1L, 2L, 1L, 1L, 2L, 1L, 2L, 1L, 2L, 1L, 1L, 2L, 2L, 
1L, 1L, 2L, 1L, 2L, 2L, 2L), levels = c("A", "B"), class = "factor"), 
Age_centered = c(-6.15315315315316, 12.8468468468468, -9.15315315315316, 
14.8468468468468, 12.8468468468468, 2.84684684684684, -10.1531531531532, 
-18.1531531531532, -1.15315315315316, 8.84684684684684, -17.1531531531532, 
13.8468468468468, 9.84684684684684, 17.8468468468468, -19.1531531531532, 
-7.15315315315316, -10.1531531531532, -19.1531531531532, 
-7.15315315315316, 0.846846846846844, -17.1531531531532, 
5.84684684684684, -25.1531531531532, -1.15315315315316, -15.1531531531532, 
4.84684684684684, 1.84684684684684, 12.8468468468468, -11.1531531531532, 
5.84684684684684, -6.15315315315316, -0.153153153153156, 
20.8468468468468, 5.84684684684684, -0.153153153153156, 12.8468468468468, 
-19.1531531531532, -11.1531531531532, 1.84684684684684, 0.846846846846844, 
-21.1531531531532, 9.84684684684684, 15.8468468468468, 14.8468468468468, 
-12.1531531531532, -11.1531531531532, -9.15315315315316, 
5.84684684684684, -4.15315315315316, 12.8468468468468, 1.84684684684684, 
-7.15315315315316, -3.15315315315316, 7.84684684684684, 0.846846846846844, 
-4.15315315315316, 5.84684684684684, -0.153153153153156, 
1.84684684684684, -7.15315315315316, 1.84684684684684, -9.15315315315316, 
6.84684684684684, 9.84684684684684, 17.8468468468468, 5.84684684684684, 
9.84684684684684, -10.1531531531532, -5.15315315315316, 18.8468468468468, 
21.8468468468468, -0.153153153153156, 2.84684684684684, -8.15315315315316, 
-5.15315315315316, 5.84684684684684, 2.84684684684684, -15.1531531531532, 
2.84684684684684, 25.8468468468468, -11.1531531531532, 27.8468468468468, 
2.84684684684684, 20.8468468468468, -0.153153153153156, -2.15315315315316, 
12.8468468468468, -0.153153153153156, 0.846846846846844, 
11.8468468468468, -8.15315315315316, 3.84684684684684, 22.8468468468468, 
5.84684684684684, 12.8468468468468, 4.84684684684684, 11.8468468468468, 
-5.15315315315316, -17.1531531531532, -7.15315315315316, 
-16.1531531531532, 0.846846846846844, -13.1531531531532, 
-13.1531531531532, -19.1531531531532, -15.1531531531532, 
-11.1531531531532, -4.15315315315316, -0.153153153153156, 
-4.15315315315316, -7.15315315315316)), row.names = c(NA, 
-111L), class = "data.frame")
```
