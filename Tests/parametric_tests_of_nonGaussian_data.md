# Intro

When learning statistical inference, you've probably heard many (false) claims that all parametric tests require a normal distribution.

This is incomprehensible to me. The very word "parametric" conceals a trivially simple explanation: parametric - based on the parameters of a certain theoretical distribution!

Do you see the words "exclusively the normal distribution" in the sentence above? Of course not! Because it applies to any distribution! If you believe that your data is well described by a certain theoretical distribution, you can try to "fit" it to the data, i.e., find its parameters. You can use, for example, the Maximum Likelihood Method. With these parameters, you can perform inference based on that distribution. For the normal distribution, these are the mean and variance. For the Poisson distribution, the lambda parameter, which acts as the mean and variance. For the beta distribution - two parameters that determine its shape...

And NO - this isn't about "for a sufficiently large N, the central limit theorem (CLT) causes the distribution of means to converge to the normal distribution...". That's not the point at all! It doesn't matter whether N = 5 or 50k, and we're not just talking about the Gosset (Student) t-test and "approximating the distribution of means with the normal distribution"! We're talking about testing *using* a specific distribution --> directly!

Have you ever heard of the Generalized Linear Model (GLM)? If so, you know that it comprises a family of regression models, such as: linear regression (general linear model), Poisson regression, negative binomial regression, quasi-Poisson regression, inverse Gaussian regression, logistic and probit regression, logistic multinomial regression, logistic ordinal regression, gamma regression, beta regression (here - an extension for a distribution with >1 parameter - VGLM [vector GLM])... And you probably know that the beta distribution family includes left-skewed, right-skewed, U-shaped, and other asymmetric distributions? Well, they certainly don't resemble a Gaussian bell curve :)

If you work in epidemiology, you've probably used the Poisson parametric test for standardized mortality ratios (SMRs) more than once? Or perhaps the binomial test for comparing two SMRs? Or for inferring a proportion? But the Poisson or binomial distribution is definitely not a normal distribution! Not only does the Poisson distribution apply to discrete data, and even to counts, but the mean is also dependent on (equal to) the variance, unlike the normal distribution, where both parameters are independent!

Have you heard of parametric survival analysis? If so, you probably know that it uses distributions like the Gompertz, Weibull, lognormal, and exponential? And again - these are far from bell curves!

Have you heard of the log-linear model for analyzing... contingency tables? Yes, yes, it's the chi-square test!

Unfortunately, sometimes even such excellent monographs as Sheskin's Handbook of Parametric and Non-Parametric Methods mislead learners, writing, for example, that logistic regression... is not parametric because... it does not follow the normal distribution.

Professors Nelder and Wedderburn, the creators of the Generalized Linear Model, Professors McFadden and Berkson, Sir David Cox - they all are turning in their graves hearing these nonsenses...

------------------

# Examples

Okay, so let's see some parametric tests on data with a very... NOT NORMAL distributions.

(red curve/line - density/QQ of the Gaussian distribution)


``` r
library(ggplot2)
library(dplyr)
library(ggh4x)
library(patchwork)
library(qqplotr)
library(marginaleffects)
library(rlang)

N <- 100

plot_data <- function(data, title, bins = 10) {
  (data %>% 
     ggplot(aes(x=sample)) + 
     geom_histogram(bins = 10, aes(y=after_stat(density))) + 
     theme_bw() + 
     facet_wrap(~ group) +
     ggtitle(title) + 
     stat_theodensity(distri="norm", col="red")) /
    (data %>% 
       ggplot(aes(sample=sample)) + 
       theme_bw() + 
       stat_qq_band(bandType = "ts") + 
       stat_qq_point() + 
       stat_qq_line(col="red") + 
       geom_label(data = data %>%
                    group_by(group) %>%
                    summarize(SW = sprintf("Shapiro-Wilk: %s",
                                           scales::pvalue(shapiro.test(sample)$p.value, add_p = TRUE)),
                              .groups = "drop"),
                  aes(x=I(0.1), y=I(0.9), label = SW, group=group), inherit.aes = FALSE, hjust="left") +
       facet_wrap(~ group) +
       labs(caption = "Aldor-Noiman tail-sensitive simult. conf. band",
            x = "Theoretical Quantiles", y = "Sample Quantiles"))
}
```
------------------

## Poisson

``` r
set.seed(1000)
poisson_data <- data.frame(sample = c(rpois(N/2, 0.5), rpois(N/2, 1.2)),
                           group  = rep(c("A", "B"), each=N/2))

plot_data(poisson_data, title = "Poisson")
```

<img width="672" height="480" alt="obraz" src="https://github.com/user-attachments/assets/c3531c71-a35b-4ddb-b8dc-82c11099a23c" />

``` r
model_full <- glm(sample ~ group, family = poisson(link = "log"), data = poisson_data)
model_null <- glm(sample ~ 1, family = poisson(link = "log"), data = poisson_data)

# Wald:
coef(summary(model_full))[2,]
```

```
##     Estimate   Std. Error      z value     Pr(>|z|) 
## 0.9808292482 0.2558739238 3.8332520701 0.0001264603
```

``` r
# LRT
anova(model_null, model_full, test = "LRT")
```

```
## Analysis of Deviance Table
## 
## Model 1: sample ~ 1
## Model 2: sample ~ group
##   Resid. Df Resid. Dev Df Deviance  Pr(>Chi)    
## 1        99     116.89                          
## 2        98     100.38  1   16.508 4.845e-05 ***
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
```

``` r
# Rao
anova(model_null, model_full, test = "Rao")
```

```
## Analysis of Deviance Table
## 
## Model 1: sample ~ 1
## Model 2: sample ~ group
##   Resid. Df Resid. Dev Df Deviance    Rao  Pr(>Chi)    
## 1        99     116.89                                 
## 2        98     100.38  1   16.508 15.909 6.644e-05 ***
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
```
``` r
#---------------------------------------------------
# Exemplary diagnostics
# do NOT use deviance residuals for discrete data:
car::qqPlot(residuals(model_full, type = "deviance"))
```
<img width="672" height="480" alt="obraz" src="https://github.com/user-attachments/assets/c0830198-3872-45b6-a224-650a8f9c5914" />

```
## [1] 41 34
```

``` r
# only the DHARMa (permuted quantile residuals)
DHARMa::testResiduals(DHARMa::simulateResiduals(model_full))
```
<img width="672" height="480" alt="obraz" src="https://github.com/user-attachments/assets/ef913c31-3428-422f-b1de-00ed7befddf4" />

```
## $uniformity
## 
##  Asymptotic one-sample Kolmogorov-Smirnov test
## 
## data:  simulationOutput$scaledResiduals
## D = 0.040587, p-value = 0.9965
## alternative hypothesis: two-sided
## 
## 
## $dispersion
## 
##  DHARMa nonparametric dispersion test via sd of residuals fitted vs.
##  simulated
## 
## data:  simulationOutput
## dispersion = 0.92236, p-value = 0.776
## alternative hypothesis: two.sided
## 
## 
## $outliers
## 
##  DHARMa bootstrapped outlier test
## 
## data:  simulationOutput
## outliers at both margin(s) = 0, observations = 100, p-value = 1
## alternative hypothesis: two.sided
##  percent confidence interval:
##  0.00 0.01
## sample estimates:
## outlier frequency (expected: 0.0028 ) 
##                                     0
```
------------------

### Gamma 1

``` r
set.seed(1000)
gamma1_data <- data.frame(sample = c(rgamma(N/2, shape = 1, rate = 2), rgamma(N/2, shape = 0.5, rate = 1)),
                           group  = rep(c("A", "B"), each=N/2))


plot_data(gamma1_data, title = "Gamma")
```

<img width="672" height="480" alt="obraz" src="https://github.com/user-attachments/assets/f0cfc237-b615-491f-905a-da8686ff4c60" />

``` r
model_full <- glm(sample ~ group, family = Gamma(link = "inverse"), data = gamma1_data)
model_null <- glm(sample ~ 1,     family = Gamma(link = "inverse"), data = gamma1_data)

# Wald:
coef(summary(model_full))[2,]
```

```
##   Estimate Std. Error    t value   Pr(>|t|) 
##  0.1357780  0.4575729  0.2967353  0.7672966
```

``` r
# LRT
anova(model_null, model_full, test = "LRT")
```

```
## Analysis of Deviance Table
## 
## Model 1: sample ~ 1
## Model 2: sample ~ group
##   Resid. Df Resid. Dev Df Deviance Pr(>Chi)
## 1        99     196.73                     
## 2        98     196.59  1  0.13367   0.7664
```

``` r
# Rao
anova(model_null, model_full, test = "Rao")
```

```
## Analysis of Deviance Table
## 
## Model 1: sample ~ 1
## Model 2: sample ~ group
##   Resid. Df Resid. Dev Df Deviance     Rao Pr(>Chi)
## 1        99     196.73                             
## 2        98     196.59  1  0.13367 0.13358   0.7665
```
------------------

### Beta 1

``` r
set.seed(1000)
beta1_data <- data.frame(sample = c(rbeta(N/2, shape1 = 0.1, shape2 = 0.1), rbeta(N/2, shape1 = 1, shape2 = 0.1)),
                           group  = rep(c("A", "B"), each=N/2))


plot_data(beta1_data, title = "Beta 1")
```

<img width="672" height="480" alt="obraz" src="https://github.com/user-attachments/assets/b15af46e-0c20-4038-9787-b9f2157b1808" />

``` r
model_full <- betareg::betareg(sample ~ group, data = beta1_data)
model_null <- betareg::betareg(sample ~ 1    , data = beta1_data)
```

```
## Warning in betareg.fit(X, Y, Z, weights, offset, link, link.phi, type,
## control): no valid starting value for precision parameter found, using 1
## instead
```

``` r
# Wald:
lmtest::waldtest(model_null, model_full)
```

```
## Wald test
## 
## Model 1: sample ~ 1
## Model 2: sample ~ group
##   Res.Df Df  Chisq Pr(>Chisq)    
## 1     98                         
## 2     97  1 22.045  2.663e-06 ***
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
```

``` r
# LRT
lmtest::lrtest(model_null, model_full)
```

```
## Likelihood ratio test
## 
## Model 1: sample ~ 1
## Model 2: sample ~ group
##   #Df LogLik Df  Chisq Pr(>Chisq)    
## 1   2 550.35                         
## 2   3 561.14  1 21.574  3.405e-06 ***
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
```
------------------

### Beta 2

``` r
set.seed(1000)
beta2_data <- data.frame(sample = c(rbeta(N/2, shape1 = 1, shape2 = 0.1), rbeta(N/2, shape1 = 0.1, shape2 = 1)),
                           group  = rep(c("A", "B"), each=N/2))

beta2_data[beta2_data$sample == 1, "sample"] <- beta2_data[beta2_data$sample == 1, "sample"] - 0.0000000001

plot_data(beta2_data, title = "Beta 2")
```

<img width="672" height="480" alt="obraz" src="https://github.com/user-attachments/assets/bbbaafaf-d015-4f44-a85f-7fc5180296e6" />

``` r
model_full <- betareg::betareg(sample ~ group, data = beta2_data, )
model_null <- betareg::betareg(sample ~ 1    , data = beta2_data)
```

```
## Warning in betareg.fit(X, Y, Z, weights, offset, link, link.phi, type,
## control): no valid starting value for precision parameter found, using 1
## instead
```

``` r
# Wald:
lmtest::waldtest(model_null, model_full)
```

```
## Wald test
## 
## Model 1: sample ~ 1
## Model 2: sample ~ group
##   Res.Df Df  Chisq Pr(>Chisq)    
## 1     98                         
## 2     97  1 110.32  < 2.2e-16 ***
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
```

``` r
# LRT
lmtest::lrtest(model_null, model_full)
```

```
## Likelihood ratio test
## 
## Model 1: sample ~ 1
## Model 2: sample ~ group
##   #Df LogLik Df  Chisq Pr(>Chisq)    
## 1   2 515.26                         
## 2   3 568.10  1 105.69  < 2.2e-16 ***
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
```
------------------

### Bernoulli

``` r
set.seed(1000)
bern_data <- data.frame(sample = c(rbinom(N/2, 1, 0.2), rbinom(N/2, 1, 0.5)),
                           group  = rep(c("A", "B"), each=N/2))


plot_data(bern_data, title = "Bernoulli")
```

<img width="672" height="480" alt="obraz" src="https://github.com/user-attachments/assets/0f4d0d4d-293e-46c9-b4dd-f6660f37f8b8" />

``` r
model_full <- glm(sample ~ group, family = binomial(link = "logit"), data = bern_data)
model_null <- glm(sample ~ 1,     family = binomial(link = "logit"), data = bern_data)

# Wald:
coef(summary(model_full))[2,]
```

```
##     Estimate   Std. Error      z value     Pr(>|z|) 
## 2.1171818527 0.5498307037 3.8506068114 0.0001178255
```

``` r
#... maybe also on the probability scale?
update(pairs(emmeans::emmeans(model_full, specs = ~ group, regrid="response"), reverse = TRUE), infer=c(TRUE, TRUE))
```

```
##  contrast estimate     SE  df asymp.LCL asymp.UCL z.ratio p.value
##  B - A        0.38 0.0824 Inf     0.218     0.542   4.611  <.0001
## 
## Confidence level used: 0.95
```

``` r
#... or this way:
bern_data$grp = bern_data$group
avg_slopes(glm(sample ~ grp, family = binomial(link = "logit"), data = bern_data[, -2]))
```

```
## 
##  Estimate Std. Error    z Pr(>|z|)    S 2.5 % 97.5 %
##      0.38     0.0824 4.61   <0.001 17.9 0.218  0.542
## 
## Term: grp
## Type: response
## Comparison: B - A
```

``` r
# LRT
anova(model_null, model_full, test = "LRT")
```

```
## Analysis of Deviance Table
## 
## Model 1: sample ~ 1
## Model 2: sample ~ group
##   Resid. Df Resid. Dev Df Deviance Pr(>Chi)    
## 1        99     120.43                         
## 2        98     101.74  1   18.687 1.54e-05 ***
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
```

``` r
# Rao
anova(model_null, model_full, test = "Rao")
```

```
## Analysis of Deviance Table
## 
## Model 1: sample ~ 1
## Model 2: sample ~ group
##   Resid. Df Resid. Dev Df Deviance    Rao  Pr(>Chi)    
## 1        99     120.43                                 
## 2        98     101.74  1   18.687 17.533 2.824e-05 ***
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
```
------------------

### Inverse Gaussian

``` r
set.seed(1000)
invg_data <- data.frame(sample = c(rmutil::rinvgauss(N/2, .5, 2), rmutil::rinvgauss(N/2, 1, 2)),
                           group  = rep(c("A", "B"), each=N/2))
```

```
## Registered S3 method overwritten by 'rmutil':
##   method         from
##   print.response httr
```

``` r
plot_data(invg_data, title = "Inverse Gaussian")
```

<img width="672" height="480" alt="obraz" src="https://github.com/user-attachments/assets/cc061a1d-b412-4242-bddc-a1ff511f0716" />

``` r
model_full <- glm(sample ~ group, family = inverse.gaussian(link = "1/mu^2"), data = invg_data)
model_null <- glm(sample ~ 1,     family = inverse.gaussian(link = "1/mu^2"), data = invg_data)

# Wald:
coef(summary(model_full))[2,]
```

```
##    Estimate  Std. Error     t value    Pr(>|t|) 
## -3.38491902  1.42179641 -2.38073398  0.01921173
```

``` r
# LRT
anova(model_null, model_full, test = "LRT")
```

```
## Analysis of Deviance Table
## 
## Model 1: sample ~ 1
## Model 2: sample ~ group
##   Resid. Df Resid. Dev Df Deviance Pr(>Chi)   
## 1        99     209.28                        
## 2        98     195.28  1   13.993 0.007933 **
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
```

``` r
# Rao
anova(model_null, model_full, test = "Rao")
```

```
## Analysis of Deviance Table
## 
## Model 1: sample ~ 1
## Model 2: sample ~ group
##   Resid. Df Resid. Dev Df Deviance    Rao Pr(>Chi)  
## 1        99     209.28                              
## 2        98     195.28  1   13.993 12.857  0.01093 *
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
```
------------------

### data suitable for quasi binomial (fractional logit)

``` r
set.seed(1000)
quasibin_data <- data.frame(sample = c(runif(N/2, 0, 0.6), runif(N/2, 0.1, 0.6)),
                           group  = rep(c("A", "B"), each=N/2))


plot_data(quasibin_data, title = "quasi binomial")
```

<img width="672" height="480" alt="obraz" src="https://github.com/user-attachments/assets/64d57f10-c0cc-4928-82c8-8ff7c87273be" />

``` r
model_full <- glm(sample ~ group, family = quasibinomial(link = "logit"), data = quasibin_data)
model_null <- glm(sample ~ 1,     family = quasibinomial(link = "logit"), data = quasibin_data)

# Wald:
coef(summary(model_full))[2,]
```

```
##   Estimate Std. Error    t value   Pr(>|t|) 
##  0.2236396  0.1459759  1.5320313  0.1287366
```

``` r
# LRT
anova(model_null, model_full, test = "LRT")
```

```
## Analysis of Deviance Table
## 
## Model 1: sample ~ 1
## Model 2: sample ~ group
##   Resid. Df Resid. Dev Df Deviance Pr(>Chi)
## 1        99     12.576                     
## 2        98     12.307  1  0.26843   0.1251
```

``` r
# Rao
anova(model_null, model_full, test = "Rao")
```

```
## Analysis of Deviance Table
## 
## Model 1: sample ~ 1
## Model 2: sample ~ group
##   Resid. Df Resid. Dev Df Deviance     Rao Pr(>Chi)
## 1        99     12.576                             
## 2        98     12.307  1  0.26843 0.26823   0.1252
```
------------------

### Negative binomial

``` r
set.seed(1000)
negbin_data <- data.frame(sample = c(MASS::rnegbin(N/2, 5, 1), MASS::rnegbin(N/2, 10, 1)),
                           group  = rep(c("A", "B"), each=N/2))


plot_data(negbin_data, title = "Negative binomial")
```

<img width="672" height="480" alt="obraz" src="https://github.com/user-attachments/assets/e6f4b324-0d71-4aed-8c74-c3ae60622e7f" />

``` r
model_full <- MASS::glm.nb(sample ~ group, data = negbin_data)
model_null <- MASS::glm.nb(sample ~ 1, data = negbin_data)

# Wald:
coef(summary(model_full))[2,]
```

```
##   Estimate Std. Error    z value   Pr(>|z|) 
## 0.53179625 0.20760009 2.56163791 0.01041799
```

``` r
# LRT
anova(model_null, model_full)
```

```
## Likelihood ratio tests of Negative Binomial Models
## 
## Response: sample
##   Model     theta Resid. df    2 x log-lik.   Test    df LR stat.    Pr(Chi)
## 1     1 0.9885627        99       -631.6261                                 
## 2 group 1.0570542        98       -625.3037 1 vs 2     1 6.322369 0.01192242
```

------------------

### A bonus: contingency table analysed via the Chi2 test vs Poisson GLM + Rao score test

``` r
set.seed(1000)
group   <- factor(rep(c("A","B","C"), each = 100))
outcome <- factor(sample(c("Yes","No"), 300, replace = TRUE, prob = c(0.2, 0.6)))

(tab <- table(group, outcome))
```

```
##      outcome
## group No Yes
##     A 79  21
##     B 71  29
##     C 75  25
```

``` r
chisq.test(tab)
```

```
## 
## 	Pearson's Chi-squared test
## 
## data:  tab
## X-squared = 1.7067, df = 2, p-value = 0.426
```

``` r
df <- as.data.frame(tab)
model_null <- glm(Freq ~ group + outcome, family = poisson, data = df)
model_full <- glm(Freq ~ group * outcome, family = poisson, data = df)
anova(model_null, model_full, test = "Rao")
```

```
## Analysis of Deviance Table
## 
## Model 1: Freq ~ group + outcome
## Model 2: Freq ~ group * outcome
##   Resid. Df Resid. Dev Df Deviance    Rao Pr(>Chi)
## 1         2     1.7124                            
## 2         0     0.0000  2   1.7124 1.7067    0.426
```
------------------

### A bonus: multinomial logistic regression (yes, it's still parametric!):

``` r
set.seed(1000)

mult_data <- data.frame(sample <- factor(sample(c("red","blue","green"), 2*N, replace=TRUE)),
                        group   <- factor(rep(c("A","B"), each=N)))

library(nnet)
```

```
## Warning: package 'nnet' was built under R version 4.2.3
```

``` r
model_full <- multinom(sample ~ group, data = mult_data)
```

```
## # weights:  9 (4 variable)
## initial  value 219.722458 
## final  value 217.156883 
## converged
```

``` r
model_null <- multinom(sample ~ 1, data = mult_data)
```

```
## # weights:  6 (2 variable)
## initial  value 219.722458 
## final  value 218.759943 
## converged
```

``` r
broom::tidy(model_full, conf.int = TRUE, exponentiate = TRUE)
```

```
## # A tibble: 4 × 8
##   y.level term        estimate std.error statistic p.value conf.low conf.high
##   <chr>   <chr>          <dbl>     <dbl>     <dbl>   <dbl>    <dbl>     <dbl>
## 1 green   (Intercept)    0.943     0.243    -0.243  0.808     0.586      1.52
## 2 green   groupB         1.86      0.353     1.75   0.0795    0.930      3.70
## 3 red     (Intercept)    0.914     0.245    -0.366  0.714     0.566      1.48
## 4 red     groupB         1.55      0.362     1.21   0.226     0.762      3.15
```

``` r
anova(model_null, model_full, test = "Chisq")
```

```
## Likelihood ratio tests of Multinomial Models
## 
## Response: sample
##   Model Resid. df Resid. Dev   Test    Df LR stat.   Pr(Chi)
## 1     1       398   437.5199                                
## 2 group       396   434.3138 1 vs 2     2 3.206119 0.2012797
```
