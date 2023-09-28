# Introduction
Testing hypotheses through statistical models opens a universe of new possibilities. Learn how to improve your daily work with this approach.

[In this document](https://github.com/adrianolszewski/Logistic-regression-is-regression/blob/main/Testing%20hypotheses%20about%20proportions%20using%20logistic%20regression.md) I showed how the model-based approach, namely
the logistic regression (with extensions) can replicate numerous classic tests of proportions: 
Wald's and Rao z test for 2+ proportions (+ ANOVA-like), Cochran-Mantel-Haenszel (CMH), Breslow-Day, Cochran-Armitage, McNemar, Cochran Q, Friedman, Mann-Whitney (-Wilcoxon), Kruskal-Wallis.

But this extends much more! General and Generalized Linear Models fit via Generalized Least Square, Generalized Estimating Equations and through Generalized Linear Mixed Models, (Mixed-effect) Quantile regression, Cox proportional-hazard regression, Accelerated Failure Time, Tobit regression, Simplex regression and dozens of others models will be at your disposal to test hypotheses in a way you probably were never told.

## Just one question: WHY?
Well, classic tests are "simple" and fast. But simple method is for simple scenarios.
A more advanced inferential analysis often goes FAR beyond that these tests can do.

For our considerations it's important to say, that **by applying Likelihood Ratio, Rao's, or Wald's testing procedure to an appropriate model you will be (mostly) able to reproduce what the classic "plain" tests do (and much more!)**.

/ PS: Oh my, so many related topics!
I'm not going to write a textbook of statistical methods, so if you're curious about the Wilk's Likelihood Ratio (LR), Wald's and Rao testing, google these terms.
You can start from:
- [Likelihood-ratio test or z-test?](https://stats.stackexchange.com/questions/48206/likelihood-ratio-test-or-z-test)
- [Mathematical Statistics — Rigorous Derivations and Analysis of the Wald Test, Score Test, and Likelihood Ratio Test. Derivations of the Classic Trinity of Inferential Tests with Full Computational Simulation](https://towardsdatascience.com/mathematical-statistics-a-rigorous-derivation-and-analysis-of-the-wald-test-score-test-and-6262bed53c55)
- [FAQ: How are the likelihood ratio, Wald, and Lagrange multiplier (score) tests different and/or similar?](https://stats.oarc.ucla.edu/other/mult-pkg/faq/general/faqhow-are-the-likelihood-ratio-wald-and-lagrange-multiplier-score-tests-different-andor-similar/)
- (warning, math!) [Handouts for Lecture Stat 461-561 Wald, Rao and Likelihood Ratio Tests](https://www.cs.ubc.ca/~arnaud/stat461/lecture_stat461_WaldRaoLRtests_handouts_2008.pdf)
- (warning, math!) [STAT 713 MATHEMATICAL STATISTICS II - Lecture Notes](https://people.stat.sc.edu/Tebbs/stat713/s18notes.pdf)

**Take home messages:**
1. In our considerations we'll primarily employ two statistical methods: Wald's test and Wilk's Likelihood Ratio Test (LRT). We'll closely examine how these methods behave in various scenarios.

2. For you now it's important to know, that asymptotically these methods yield equivalent results. Asymptotically means "at infinite sample size". In real-world scenarios they will always differ, but unless you hit the "edge case", they will be rather consistent.
  
3. But the MAY also diverge **noticeably** in "edge cases", where the log-likelihood curve of a model parameter in the parameter space deviates from a parabolic shape. If you read either of the 2 first links, you will know what I mean. We will observe it...

4. Wald's test is extremely useful for testing hypotheses in both
   a. covariate-adjusted planned comparisons via contrasts - **that's the ESSENCE of work in experimental research!**
   b. AN[C]OVA-type joint testing of model coefficients assessing the main and interaction effects (like the classic ANOVA does)

5. Wald's is faster than LR as it needs only a single model fitting.
  
6. LRT is often considered more conservative than Wald's test. It's less likely to reject a null hypothesis when it's true (i.e., it has a lower Type I error rate). This conservativeness can be advantageous when you want to avoid making false-positive errors, such as in clinical trials or other critical applications. In contrast, Wald's test can sometimes be more liberal, leading to a higher likelihood of false positives. It may not perform well when sample sizes are small or when the distribution of the parameter estimates deviates from normality. But often Wald's is the best we can do.
   
7. Wald's is also the ONLY available method in non-likelihood models, like GEE estimation or quantile regression. Depiste it weaknesses, it gains importance!
  
8. But sometimes Wald's testing fails to calculate (e.g. estimation of covariance fails), while the likelihood is still obtainable and then the LRT is the only method that works. Who said the world is simple? :)

9. Likelihood Ratio is employed through model comparisons: one model with the term we want to assess and one model without it (you can also compare other settings, like covariance structures) and is often found as **more accurrate** than Wald's.

10. LRT allows one for performing AN[C]OVA-type analyse (which requires careful specification of the model terms!) but doesn't help in covarite-adjusted planned comparisons.

11. LRT involves fitting two models, doubling the time needed compared to Wald's test. This additional computational cost can be a consideration. Sometimes it really MAKES an inssue, especially if you need to test several predictor variables, so you need to compare more and more nested models!

12. Wald's approach takes full advantage of the estimated parameter and covariance matrix, which means that "sky is the limit" when testing.
    
13. To conclude: for simple testing 2+ groups (single n-level categorical predictor) the LRT is usually preferred. For planned and exploratory assessment of various contrasts - Wald's is often the method of choice. But there is no univerally good scenario.

14. "_When the sample size is small to moderate, the Wald test is the least reliable of the three tests. We should not trust it for such a small n as in this example (n = 10). Likelihood-ratio inference and score-test based inference are better in terms of actual error probabilities coming close to matching nominal levels. A marked divergence in the values of the three statistics indicates that the distribution of the ML estimator may be far from normality. In that case, small-sample methods are more appropriate than large-sample methods._ (Agresti, A. (2007). An introduction to categorical data analysis (2nd edition). John Wiley & Sons.)"

15. Now, a mandatory reading! [Wald vs likelihood ratio test](https://thestatsgeek.com/2014/02/08/wald-vs-likelihood-ratio-test/), which perfectly summarized our above considerations: "_In conclusion, although the likelihood ratio approach has clear statistical advantages, computationally the Wald interval/test is far easier. In practice, provided the sample size is not too small, and the Wald intervals are constructed on an appropriate scale, they will usually be reasonable (hence their use in statistical software packages). In small samples however the likelihood ratio approach may be preferred._"
/

This way, by employing appropriate model followed by the metioned testing procedures, you can obtain perform complex planned and exploratory hypothesis testing, perform  Type-II and Type-III ANOVA-like analysis over Poisson or logistic regression... Actually, you can plug-in any single model.

So what?

Think about it. **By using a more general model, that can handle several "problems" in data, you can bypass various assumptions of classic parametric tests RATHER than switching to non-parametric methods** (which have often complex interpretation and entirely change your null hypotheses!).

Namely, you can address:
- **mean-variance dependency and general heteroscedasticity**? How about the _Generalized Linear Model (GLM)_, _Generalized Least Square (GLS)_ or _Generalized Estimating Equation (GEE)_ estimation?
- **dependency in residuals** (repeated observations, clusetered data)? How about the _GLS_, _GEE_ and _mixed models_ (LMM, GLMM)?
- **non-normal residuals**? Maybe the _GLM_ or _GEE_ (doesn't make distributional assumptions about the residuals)? Oh, OK, GEE is semi-parametric, but still retains your null hypotheses about conditional means!
- **conditional distributions are skewed**? (related to the above one) Again - GLM: gamma or beta regression, maybe Poisson or negative binomial? Or fractional logit?
- **categorical data**? Again, GLM (with extensions) will offer you the ordinal logistic regression (proportional-odds model), partial ordinal regression, multinomial logistic regression
- **relationships (hierarchy) in the response**? How about the nested logistic regression?
- **survival data**? I guess you know that Log-rank test is equivalent to using Cox regression with a 2-level categorical predictor.
And that's not all!

How about the _quantile regression_ (especially combined with random effects) can handle even more! And yes - it's a non-parametric method, BUT at least well interpretable (what cannot be said about Mann-Whitney (-Wilcoxon), for instance.

**But that's not all - try doing THIS with classic tests!**

1. ordinary tests won't allow you to control for covariates. Goodbye more advanced analyses.
2. most classic tests cannot handle more complex m x n-way designs. Want to test some outcome across multiple groups rendered by `sex * visit * treatment_arm`? Forget! You will likely need to run multiple simple tests and they won't be able to detect inter-variable relationships. 
3. most of the classic tests won't be able to test interactions between multiple factors (a few modern ones, like ATS (ANOVA-Type Statistic) or WTS (Wald-Type Statistic) can do this, but only in a limited scope (1-level interaction between just 2 factors).
4. classic tests won't allow you to test simple effects via contrasts, e.g.: "Visit: 3, Arm: Treatment | Male vs Female effect" vs. "Visit 3, Arm: Control | Male vs Female effect". **For the model-based testing it's a piece of cake**.
5. you may simply NOT KNOW which test to use! Believe me or not, there are 850+ (EIGHT HUNDRED FIFTY!) statistical tests and counting. A colleague of mine has been counting them for years (with my little support). With a model - you don't care about all these "version", just provide the formula, set some parameters and test the hypotheses you need. Need a trest for trend? Just user ordinal factor for your time variable.
6. and you will obtain standard errors and confidence intervals for these comparisons!
7. want to test some hypotheses jointly? Forget with classic tests!
8. want to effectively adjust for multiple testing using parametric exact method employing the estimated effects and covariaces through the multivariate t distribution ("MVT")? This is far better than Bonferroni :-) But forget this flexibility when using plain tests!

Fair enough?

## But the model-based testing is SLOW!
Let's be honest - model-based testing can be (and IS) **SLOWER** than running a plain good-old test, especially if you perform then under multiple imputation approach to missing data (where you repeat the same analysis on each of dozens of imputed datasets and then pool the results.), especially if you also employ the Likelihood Ratio testing (LRT).

But if you have just a few ten-to-hundred of observations, then I wouldn't be much concerned.

**But - well - did you expect a free lunch?
As I said in the beginning - simple tests do simple jobs and cannot do anything more. If you want more, you have to pay for it. That's as simple.**

----
## Not only it's slow. It may also fail to converge!
Yes, they may and sometimes they do.

And while tests may fail computationally too, it happens incomparably rarer than with models.
Just recall how many times you saw messages like "failed to converge", "did not converge", "negative variance estimate" or similar, and how often the plain test failed? Well, that's it.

Models are non-trivial procedures and may computationally fail under "edge" conditions. For exmaple, if you fit logistic regression with all responses equal to TRUE or FALSE, then - well... it's rather obvious that it will fail, right?

Let's take an example of the ordinal logistic regression (aka proportional-odds model) to replicate the Mann-Whitney (-Wilcoxon) test
Let me show you a trivial case: we will compare two samples from a normal distribution: N(mean=0, SD=1) vs. N(mean=5, SD=1).
What can be simpler than that?

/ PS: yes, YOU CAN use the ordinal logistic regression for continuous numerical data. You will have just dozens of intercepts :)
Read: [If you like the Wilcoxon test you must like the proportional odds model](https://www.fharrell.com/post/wpo/), [Equivalence of Wilcoxon Statistic and Proportional Odds Model](https://www.fharrell.com/post/powilcoxon/), and (best repository on this topic I've ever seen!) [Resources for Ordinal Regression Models](https://www.fharrell.com/post/rpo/) /

``` r
> set.seed(1000)
> stack(
+     data.frame(arm1 = rnorm(50, mean=0),
+                arm2 = rnorm(50, mean=5))) %>% 
+     mutate(values_ord = ordered(values), ind = factor(ind)) %>% 
+     rms::orm(values ~ ind, data=.)
Did not converge in 12 iterations
Unable to fit model using  “orm.fit” 
Error in switch(x$family, logistic = "Logistic (Proportional Odds)", probit = "Probit",  : 
  EXPR must be a length 1 vector
```
As you can see, the model failed to converge. The Mann-Whitney (-Wilcoxon) test did well, however:
``` r
> set.seed(1000)
> stack(
+      data.frame(arm1 = rnorm(50, mean=0.0),
+                 arm2 = rnorm(50, mean=5))) %>% 
+      mutate(values_ord = ordered(values), ind = factor(ind)) %>% 
+      wilcox.test(values ~ ind, data=., conf.int = TRUE)

	Wilcoxon rank sum test with continuity correction

data:  values by ind
W = 0, p-value < 2.2e-16
alternative hypothesis: true location shift is not equal to 0
95 percent confidence interval:
 -5.772126 -4.977111
sample estimates:
difference in location 
              -5.34864 
```
Since we used normal distributions, let's also check the classic t-test:
``` r
> set.seed(1000)
> stack(
+          data.frame(arm1 = rnorm(50, mean=0.01),
+                     arm2 = rnorm(50, mean=5))) %>% 
+          mutate(values_ord = ordered(values), ind = factor(ind)) %>% 
+          t.test(values ~ ind, data=., var.equal = TRUE) 

	Two Sample t-test

data:  values by ind
t = -26.803, df = 98, p-value < 2.2e-16
alternative hypothesis: true difference in means between group arm1 and group arm2 is not equal to 0
95 percent confidence interval:
 -5.735385 -4.944653
sample estimates:
mean in group arm1 mean in group arm2 
        -0.1486302          5.1913886
```
As expected. In parametric case Wilcoxon and t-test are naturally consistent.
Sometimes a model that failed can converge if we a little "perturb" the data. In our case let's add just 0.01 to the mean.
``` r
> set.seed(1000)
> stack(
+     data.frame(arm1 = rnorm(50, mean=0.01),
+                arm2 = rnorm(50, mean=5))) %>% 
+     mutate(values_ord = ordered(values), ind = factor(ind)) %>% 
+     rms::orm(values ~ ind, data=.)
Logistic (Proportional Odds) Ordinal Regression Model
 
 rms::orm(formula = values ~ ind, data = .)
 
                           Model Likelihood               Discrimination    Rank Discrim.    
                                 Ratio Test                      Indexes          Indexes    
 Obs              100    LR chi2     133.10    R2                  0.736    rho     0.865    
 Distinct Y       100    d.f.             1    R2(1,100)           0.733                     
 Median Y    2.677509    Pr(> chi2) <0.0001    R2(1,100)           0.733                     
 max |deriv|   0.0009    Score chi2   99.69    |Pr(Y>=median)-0.5| 0.487                     
                         Pr(> chi2) <0.0001                                                  
 
          Coef   S.E.   Wald Z Pr(>|Z|)
 ind=arm2 9.1781 1.7356 5.29   <0.0001 
```
----

## You said that that different method of testing (Wald's, Rao's, Wilk's LRT) may yield a bit different results?
Yes. But the differences will be mostly mostly mostly under the ALTERNATIVE hypothesis and mostly at lower p-values.

This is my observation for numerous examples. I guess that's because (if you imagine the concave (ideally quadratic) curve of log-likelihood for the model parameter representing the group indicator)
under the null hypothesis:
- **Wald's:** the difference on the horizontal axis (parameter space) = the difference between compared means (for example) -> 0
- **Wilk's LRT** the difference on the vertical axis (log-likelihoods) between both models -> 0. So the likelihood ratio -> 1.
- **Rao's score** the slope of the tangential line approaches zero (Rao' score). And this reflects the fact that under H0 the derivative (score) of the log-likelihood function with respect to the parameter should be close to zero.

While, at the alternative hypothesis, they will diverge from each other as the likelihood curve deviates from the quadratic curve.
See? It all makes sense!

The good news is that typically the differences will "manifestte" itself below most common significance levels (<0.001)

For instance, for 0-5 Likert items sampled with some predefined probability:
```r
# Under H0
# simulate_wilcox_olr_LR(samples = 100, n_group = 50, set = 0:5, 
#                        arm_1_prob =  c(20, 10, 5, 2, 2, 2),
#                        arm_2_prob =  c(20, 10, 5, 2, 2, 2),
#                        which_p = "Wald") #... and Rao, LRT 
# Under H1
# simulate_wilcox_olr_LR(samples = 100, n_group = 50, set = 0:5, 
#                        arm_1_prob =  rev(c(20, 10, 5, 2, 2, 2)),
#                        arm_2_prob =  c(20, 10, 5, 2, 2, 2),
#                        which_p = "Wald") #... and Rao, LRT
```
![obraz](https://github.com/adrianolszewski/model-based-testing-hypotheses/assets/95669100/c6b14283-caa3-4855-9024-cfcdfbebe705)

## Even worse! You just said that different method of testing can yield OPPOSITE results!

Yes, Wald's and LRT may yield opposite results, e.g. Wald's p-value = 0.8 vs. LRT p-value < 0.001.

But that doesn't happen very often, in "**edge cases**". That's why doing statistical analysis is not just a brainless "click-and-go" task!

This is summarized in the article I already cited ([Wald vs likelihood ratio test](https://thestatsgeek.com/2014/02/08/wald-vs-likelihood-ratio-test/)): "Further, a situation in which the Wald approach completely fails while the likelihood ratio approach is still (often) reasonable is when testing whether a parameter lies on the boundary of its parameter space." 

If you anticipate that something bad may happen, VALIDATE your calculations using a different method. It doesn't have to bring SAME results (it's a different method) but at least you will be able to assess if they are approximately consistent.

If, for instance, you run the ordinal logistic regression over Likert data and obtain Wald's p-value = 0.3 and LRT p-value = 0.0003, then you can also simplify it for a while and check with the Mann-Whitney (-Wilcoxon) test or even run the OLS with ranked response (it's a very good method!).

**And always - always start with the EDA (exploratory data analysis), plot your data, know it. Otherwise you won't even suspect what could go wrong.**

**Take-home messages #1:** 
1. A model can fail to converge if the data are "not nice" :-)
2. Wald's approach may fail, while the LRT will do well.
3. A test may still work in this case (like the LRT).
4. Wald's may differ from LRT under the alternative hypothesis

By the way! You may ask: _but cannot we just use a better implementations able to complete the estimation?_

Good question! Let's find a different implementation, which will converge in this problematic case. There's one MASS::polr()
Again, we start with N(0, 1) vs. N(5, 1):
```r
> set.seed(1000)
> stack(
+     data.frame(arm1 = rnorm(50, mean=0.0),
+                arm2 = rnorm(50, mean=5))) %>% 
+     mutate(values_ord = ordered(values), ind = factor(ind)) %>% 
+     MASS::polr(values_ord ~ ind, data=., Hess = TRUE) %>% summary() %>% coef() %>% as.data.frame() %>% slice(1)
           Value Std. Error  t value
indarm2 21.30274   19.75443 1.078378

# p-value
> 2*pnorm(q=1.078378, lower.tail=FALSE)
[1] 0.2808651
```
Good! This one converged!  But... wait a minute! What?! Is this a joke?! Why so big p-value?! 

**Because the p-value is simply wrong. We are exactly in the "edge case" and Wald's approach failed evidently.**

And what about the LRT? Is this any better?
```r
> set.seed(1000)
> stack(
+  data.frame(arm1 = rnorm(50, mean=0.0),
+             arm2 = rnorm(50, mean=5))) %>% 
+      mutate(values_ord = ordered(values), ind = factor(ind)) -> ord_data
 
> model_full <-  MASS::polr(values_ord ~ ind, data=ord_data, Hess = TRUE)
> model_intercept <- MASS::polr(values_ord ~ 1, data=ord_data, Hess = TRUE)
 
> anova(model_full, model_intercept)
Likelihood ratio tests of ordinal regression models

Response: values_ord
  Model Resid. df Resid. Dev   Test    Df LR stat. Pr(Chi)
1     1         1   921.0340                              
2   ind         0   782.4094 1 vs 2     1 138.6246       0
```
Yes, the LRT did better. This result is consistent with the Mann-Whitney (-Wilcoxon test).

Now let's again increase the mean in one group by 0.01. Not only the model converges properly, but also the Wald's and LRT are in a good agreement!
``` r
> set.seed(1000)
> stack(
+     data.frame(arm1 = rnorm(50, mean=0.01),
+                arm2 = rnorm(50, mean=5))) %>% 
+     mutate(values_ord = ordered(values), ind = factor(ind)) %>% 
+     MASS::polr(values_ord ~ ind, data=., Hess = TRUE) %>% summary() %>% coef() %>% as.data.frame() %>% slice(1)
           Value Std. Error  t value
indarm2 9.179891   1.735655 5.289008

# p-value
> 2*pnorm(q=5.289008, lower.tail=FALSE)
[1] 1.229815e-07

# Let's confirm with the previous method rms::orm()
# -------------------------------------------------
>  set.seed(1000)
>  stack(
+   data.frame(arm1 = rnorm(50, mean=0.01),
+              arm2 = rnorm(50, mean=5))) %>% 
+   mutate(values_ord = ordered(values), ind = factor(ind)) %>% 
+   rms::orm(values ~ ind, data=.)
Logistic (Proportional Odds) Ordinal Regression Model
 
 rms::orm(formula = values ~ ind, data = .)
 
                           Model Likelihood               Discrimination    Rank Discrim.    
                                 Ratio Test                      Indexes          Indexes    
 Obs              100    LR chi2     133.10    R2                  0.736    rho     0.865    
 Distinct Y       100    d.f.             1    R2(1,100)           0.733                     
 Median Y    2.677509    Pr(> chi2) <0.0001    R2(1,100)           0.733                     
 max |deriv|   0.0009    Score chi2   99.69    |Pr(Y>=median)-0.5| 0.487                     
                         Pr(> chi2) <0.0001                                                  
 
          Coef   S.E.   Wald Z Pr(>|Z|)
 ind=arm2 9.1781 1.7356 5.29   <0.0001 
```
**And now that's fine.**

So, as I said - even if a model *actually converges*, the estimate may be unreliable. Always remember to check your outcomes.

**Take-home messages #2:** 

3. Having multiple options, always choose the method that alerts you that something went wrong rather than method that silently pretends nothing wrong happened and happily continues. **It's always better to have NO result rather than having WRONG result.**

You. Have. Been. Warned.
----
## Interpretation
For some models the equivalence is strict, I mean - it could be shown by formulas (I was shown it in the past yet  I don't take the challenge to replicate :-) ).
For others it's a close approximation, coming from different approaches (thus, different set of calculations), testing different (yet related!) hypotheses and sometimes even rounding.
When you take a simple model like `response ~ categorical_predictor`, then, using the _treatment coding_ (default in R) the interpretation comes out of the box:
- value of the intercept is the (conditional) mean, or median, or any other function of data in the group rendered by  the reference level of your predictor (e.g. non-smokers). Note, that the mean can be either "raw"
or transformed via link function in case of the GLM (Generalized Linear Model).
- value of the 2nd model coefficient -assigned to the 2nd level of the predictor- is the difference between the 2nd group and the intercept (1st) group. 
- ... this extends to multiple predictors with >2 levels and their interactions. The calculations complicate quickly (especially in presence of interactions), but the idea remains the same.

Whether it is directly and easily interpretable depends on the formulation of your model. 

For the general linear model (say, just the linear regression) the interpretation will be easy: intercept stands for the mean in the reference group,
and the value of the "beta" (assigned to the 2nd level of predictor) represents the difference in means between both groups.
Similarly for the quantile regression yet now you get differences in medians, first, third or any other quantile you specified.

The situation complicates if you applied transformation of the dependent variable (response) or in case of the Generalized Linear Models (GLM) transforming -in contrast- the E(Y|X=x), 
like the logistic, ordinal logistic, beta and gamma regressions (+ a few more). First, the interpretation of the "conditional mean" depends
on the conditional distribution of the data (e.g. for the Bernoulli's distribution in the logistic regression the E(Y|X=x) means is probability of success). Second, it depends on the used link function, where - depending on the used transformation - differences may turn into ratios, and (say) probabilities turn into log-odds. No surprise that the interpretation of model coefficients will complicate (a bit). 

And this poses another issue. Let's take the logistic regression. While all these terms: probabilities, odds, log-odds, odds-ratios, testing hypotheses are mutually related through the well-known formulas, performing inference
on different scales may yield a little bit different results. Well, at the end of the day, you test *different* hypotheses (though related, so the results cannot be wildly discrepant!) and you will operate on different scales.
And this will complicate EVEN MORE if you will ignore some terms (average over them), because it DOES matter whether you average on the original (logit) or the transformed (probability) scale.
Change on the scale of the log-odds is not the same as change on the probability scale.
[Why is the p value by average marginal effects different than the p value of the coefficients?](https://stats.stackexchange.com/questions/464115/why-is-the-p-value-by-average-marginal-effects-different-than-the-p-value-of-the)

With the logistic regression you can compare outcomes on both logit and probability scales. The former comes out of the box, the latter requires either calculating the AME (average marginal effect) or employing the LS-means on appropriate scale.
[Transformations and link functions in emmeans](https://cran.r-project.org/web/packages/emmeans/vignettes/transformations.html)

/ If you're curious about the AME, especially in terms of the logistic regression, visit: [Usage Note 22604: Marginal effect estimation for predictors in logistic and probit models](http://support.sas.com/kb/22/604.html) | [A Beginner’s Guide to Marginal Effects](https://library.virginia.edu/data/articles/a-beginners-guide-to-marginal-effects) |
[How can I use the margins command to understand multiple interactions in logistic regression? Stata FAQ](https://stats.oarc.ucla.edu/stata/faq/how-can-i-use-the-margins-command-to-understand-multiple-interactions-in-logistic-regression-stata/) |
[Week 13: Interpreting Model Results: Marginal Effects and the margins Command](https://clas.ucdenver.edu/marcelo-perraillon/sites/default/files/attached-files/week_13_margins.pdf) | [Stata 18 - Help for margins](https://www.stata.com/help.cgi?margins) |
[An introduction to margins in R](https://cran.r-project.org/web/packages/margins/vignettes/Introduction.html)

If you are curious about the LS-means (EM-means), that's another story and just use Google.
Ideally add "SAS" to your queries, as commercial packages (SAS, Stata, SPSS, NCSS) have the best documentation I've ever seen (sorry R).
[ Getting started with emmeans ](https://aosmith.rbind.io/2019/03/25/getting-started-with-emmeans/)
[Shared Concepts and Topics - LSMEANS Statement](https://documentation.sas.com/doc/en/pgmsascdc/9.4_3.4/statug/statug_introcom_sect040.htm)
[Least-Squares Means: The R Package lsmeans](https://www.jstatsoft.org/article/view/v069i01) /

----

## Conclusions
So you can see with your own eyes, that model-based testing has LOTS of advantages. But sometimes you will need just a simple, classic test that runs FAST. Especially, if you have lots of tests to do under multiple imputation conditions, with lengthy data, and you are approaching a deadline :-)

**So the choice depends really on your needs.  I only want to show you that this is doable and how well it performs.
**
