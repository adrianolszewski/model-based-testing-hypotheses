# Introduction
Testing hypotheses through statistical models opens a universe of new possibilities. Learn how to improve your daily work with this approach.

[In this document](https://github.com/adrianolszewski/Logistic-regression-is-regression/blob/main/Testing%20hypotheses%20about%20proportions%20using%20logistic%20regression.md) I showed how the model-based approach, namely
the logistic regression (with extensions) can replicate numerous classic tests of proportions: 
Wald's and Rao z test for 2+ proportions (+ ANOVA-like), Cochran-Mantel-Haenszel (CMH), Breslow-Day, Cochran-Armitage, McNemar, Cochran Q, Friedman, Mann-Whitney (-Wilcoxon), Kruskal-Wallis

## What for?
Well, classic tests are "simple" and fast. But simple method is for simple scenarios.
A more advanced inferential analysis often goes FAR beyond that these tests can do.

For our considerations it's important to say, that **by applying Likelihood Ratio, Rao's, or Wald's testing procedure to an appropriate model you will be (mostly) able to reproduce what the classic "plain" tests do (and much more!)**.

/ PS: Oh my, so many related topics!
I'm not going to write a textbook of statistical methods, so if you're curious about the Wilk's Likelihood Ratio (LR), Wald's and Rao testing, google these terms.
You can start from:
- [Mathematical Statistics — Rigorous Derivations and Analysis of the Wald Test, Score Test, and Likelihood Ratio Test. Derivations of the Classic Trinity of Inferential Tests with Full Computational Simulation](https://towardsdatascience.com/mathematical-statistics-a-rigorous-derivation-and-analysis-of-the-wald-test-score-test-and-6262bed53c55)
- [FAQ: How are the likelihood ratio, Wald, and Lagrange multiplier (score) tests different and/or similar?](https://stats.oarc.ucla.edu/other/mult-pkg/faq/general/faqhow-are-the-likelihood-ratio-wald-and-lagrange-multiplier-score-tests-different-andor-similar/)
- [Handouts for Lecture Stat 461-561 Wald, Rao and Likelihood Ratio Tests](https://www.cs.ubc.ca/~arnaud/stat461/lecture_stat461_WaldRaoLRtests_handouts_2008.pdf)
- [STAT 713 MATHEMATICAL STATISTICS II - Lecture Notes](https://people.stat.sc.edu/Tebbs/stat713/s18notes.pdf)

Take home messages:
1. We will observe at work mostly two methods: Wald's and LRT, and see how they behave.

2. For you now it's important to know, that asymptotically they yield equivalent results.
  
3. But the MAY also diverge **noticeably** in "edge cases", where the curve representing log-likelihood of a model parameter in the parameter space doesn't form a parabola. If you read either of the 2 first links, you will know what I mean. We will observe it...

4. Wald's test is extremely useful for testing hypotheses in both
   a. covariate-adjusted planned comparisons via contrasts - **that's the ESSENCE of work in experimental research!**
   b. AN[C]OVA-type joint testing of model coefficients assessing the main and interaction effects (like the classic ANOVA does)

5. Wald's is fast (single model fitting).
  
6. Wald's may be too liberal or too conservative. There are many simulations over the Internet showing, that sometimes it's the "poorest" method and sometimes - the superior one. No easy judgement.
   
7. Wald's it's the ONLY available method in non-likelihood models, like GEE estimation of quantile regression. So better learn about it.

8. Likelihood Ratio is employed through model comparisons: one model with the term we want to assess and one model without it (you can also compare other settings, like covariance structures) and is often found as better than Wald's. It allows one for performing AN[C]OVA-type analyses, but doesn't help in covarite-adjusted planned comparisons. Requires careful specification of the model terms!

9. LRT fits 2 models, so it doubles the time needed to complete compared to Wald's. Sometimes it really MAKES an inssue.

10. For simple testing 2+ groups (single n-level categorical predictor) the LRT is usually preferred. For planned and exploratory assessment of various contrasts - Wald's is the method of choice.
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

## But it's SLOW!
Let's be honest - model-based testing can be (and IS) **SLOWER** than running a plain good-old test, especially if you perform then under multiple imputation approach to missing data (where you repeat the same analysis on each of dozens of imputed datasets and then pool the results.), especially if you use the Likelihood Ratio testing (LRT).

**But - well - did you expect a free lunch?
As I said in the beginning - simple tests do simple jobs and cannot do anything more. If you want more, you have to pay for it. That's as simple.**

----

## Not only it's slow. It may also fail to converge!
Yes, they may and sometimes they do.

And while tests may fail computationally too, it happens incomparably rarer than with models.
Just recall how many times you saw messages like "failed to converge", "did not converge", "negative variance estimate" or similar, and how often the plain test failed? Well, that's it.

Models are non-trivial procedures and may computationally fail under "edge" conditions. For exmaple, if you fit logistic regression with all responses equal to TRUE or FALSE, then - well... it's rather obvious that it will fail, right?

----

## Even worse, different method of testing can yield OPPOSITE results!

Yes, Wald's and LRT may yield opposite results, e.g. Wald's p-value = 0.8 while LRT p-value < 0.001.

But that doesn't happen very often, rather in "edge cases". That's why doing statistical analysis is not just brainless "click-and-go" task!
If you anticipate that something bad may happen, VALIDATE your calculations using a different method. It doesn't have to bring SAME results (well, it's different :-) ) but at least you will be able to assess if they are approximately consistent.

If you run the ordinal logistic regression over Likert data and obtain Wald's p-value = 0.3 and LRT p-value = 0.0003, then you can also simplify it for a while and check with the Mann-Whitney (-Wilcoxon) test or even run the OLS with ranked response (it's a very good method!).

**And always - always start with the EDA (exploratory data analysis), plot your data, know it. Otherwise you won't even suspect what could go wrong. **

Let's take an example of the ordinal logistic regression (aka proportional-odds model) to replicate the Mann-Whitney (-Wilcoxon) test
Let me show you a trivial case: we will compare two samples from a normal distribution: N(mean=0, SD=1) vs. N(mean=5, SD=1).
What can be simpler than that?
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
The model failed to converge. While Mann-Whitneh (-Wilcoxon) test did well:
``` r
> set.seed(1000)
> stack(
+     data.frame(arm1 = rnorm(50, mean=0.0),
+                arm2 = rnorm(50, mean=5))) %>% 
+     mutate(values_ord = ordered(values), ind = factor(ind)) %>% 
+     wilcox.test(values ~ ind, data=.)

	Wilcoxon rank sum test with continuity correction

data:  values by ind
W = 0, p-value < 2.2e-16
alternative hypothesis: true location shift is not equal to 0
```
Since we work with both normal distributions, let's also check the classic t-test:
``` r
> set.seed(1000)
> stack(
+     data.frame(arm1 = rnorm(50, mean=0.01),
+                arm2 = rnorm(50, mean=5))) %>% 
+     mutate(values_ord = ordered(values), ind = factor(ind)) %>% 
+     t.test(values ~ ind, data=.) %>% tidy() %>% select(p.value)
# A tibble: 1 × 1
   p.value
     <dbl>
1 1.73e-46
```
Sometimes a model CAN converge if we only a little change "perturb" the data.
In our case let's add just 0.01 to the mean.
It won't change much, still the probability of superiority (dominance) is 1 but this time it will work.

PS: I'm still investigating how it works exactly, but it does not surprise me, that altering data a little can make a model work.

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

**Take-home messages #1:** 
1. A model can fail to converge if the data are "not nice" :-)
2. A test may still work in this case

By the way! You may ask: _but cannot we just use a better implementations able to complete the estimation?_
OK, let's find a different implementation, which will converge in this problematic case
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

Cool! This one converged! 
But... wait, what?! Is this a joke?! So big p-value?! 

**Yes. Because the p-value is simply wrong. The model calculated something but the estimation is "worse than poor".**

Now let's again increase the mean in one group by 0.01
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
```
And now that's fine.
So, as I said - even if a model *actually converges*, the estimate may be unreliable. Always remember to check your outcomes.

**Take-home messages #2:** 

3. Having multiple options, always choose the method that alerts you that something went wrong rather than method that silently pretends nothing wrong happened and happily continues. **It's always better to have NO result rather than having WRONG result.**

You. Have. Been. Warned.

PS: **Wait!** You tested totally different hypotheses with these methods, so maybe that's the issue?
**NO.** The ordinal logistic regression with just 2-level categorical predictor and no covariates is equivalent to Mann-Whitney (-Wilcoxon), which tests for stochastic dominance (_NO, it's NOT about medians_). Under IID it reduces to pseudo-median difference. Under symmetry of the distributions of [Walsh averages](https://stats.stackexchange.com/questions/215889/prove-the-relationship-between-walsh-averages-and-wilcoxon-signed-rank-test) it gives median difference.
And, under normality of both distributions, the median difference = mean difference. Which is equal to difference in means (in the normal distribution mean=median)
Following me?

``` r
> set.seed(1000); x1 <- rnorm(1000000); x2 <- rnorm(1000000, mean=10)
> sprintf("mean(diff)=%f, diff(means)=%f, median(diff)=%f, diff(medians)=%f, ps-median(diff)=%f",
+ mean(x1 - x2), mean(x1) - mean(x2), median(x1 - x2), median(x1) - median(x2), wilcox.test(x1, x2, conf.int = TRUE, exact = FALSE, adjust = FALSE)$estimate)
[1] "mean(diff)=-9.999591, diff(means)=-9.999591, median(diff)=-10.000615, diff(medians)=-9.998330, ps-median(diff)=-9.999202"
```

Now let's look from another perspective. The ordinal logistic regression naturally fits the Mann-Whitney (-Wilcoxon) null hypothesis. Citing the `rms` package from Prof. Harrell: "_[orm] fits ordinal cumulative probability models for continuous or ordinal response variables [...]. The ordinal cumulative probability models are stated in terms of exceedance probabilities (P rob[Y ≥ y|X]) so that as with OLS larger predicted values are associated with larger Y._"

Check this link (may change in future, when I reorganize things) for an example of their close agreement not only in terms o p-values, but also the measure of concordance:
[Mann-Whitney (-Wilcoxon) test of stochastic equivalence (vs. stochastic superiority / dominance)](https://github.com/adrianolszewski/Logistic-regression-is-regression/blob/main/Testing%20hypotheses%20about%20proportions%20using%20logistic%20regression.md#mann-whitney--wilcoxon-test-of-stochastic-equivalence-vs-stochastic-superiority--dominance)

So actually - with just different methods - I tested the same hypotheses. Just indirectly.

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

So the choice depends really on your needs.  I only want to show you that this is doable and how well.

