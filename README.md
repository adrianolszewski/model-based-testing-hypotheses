# Introduction
Testing hypotheses through statistical models opens a universe of new possibilities. Learn how to improve your daily work with this approach.

[In this document](https://github.com/adrianolszewski/Logistic-regression-is-regression/blob/main/Testing%20hypotheses%20about%20proportions%20using%20logistic%20regression.md) I showed how the model-based approach, namely
the logistic regression (with extensions) can replicate numerous classic tests of proportions: 
Wald's and Rao z test for 2+ proportions (+ ANOVA-like), Cochran-Mantel-Haenszel (CMH), Breslow-Day, Cochran-Armitage, McNemar, Cochran Q, Friedman, Mann-Whitney (-Wilcoxon), Kruskal-Wallis

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

Now let's back on the track.
**By applying Wald's, LR (Likelihood Ratio), or Rao testing procedure to an appropriate model you will be (mostly) able to reproduce the same (or pretty equivalent) hypotheses as the ordinary statistical tests do.**

## But why? It's SLOW!
Let's be honest - model-based testing is **MUCH SLOWER** than running a plain good-old test. Just try. When I performed these simulations, at these sample sizes
tests completed calculations immediately, while for the model it took 2-3 seconds. Now add to this calculation of the AME, emmeans, some adjustments for degrees of freedom and you'll get about 5-10 seconds at just few dozens of observations. Now, imagine you run these procedures over a multiply imputed (MICE) dataset.
For just 20 imputed datasets x 5-10 seconds it's **100-200 seconds**. At work, on larger data, with more complex models, not rarely it took **5-10 minutes** to complete. In contrast, ordinary testing in this setting may take **no more than just few seconds**.

## Not only it's slow. It may also fail to converge!
And while tests may fail computationally too, it happens incomparably rarer than with models.
Just recall how many times you saw messages like "failed to converge", "did not converge", "negative variance estimate" or similar, and how often the plain test failed? Well, that's it.

Let me show you a trivial case: we will compare two samples from a normal distribution. One will have mean = 0, the other = 5.
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
As expected!

OK, so the model failed. Let's slligthly adjust the mean in the first arm by adding 0.01
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
And now it worked well.

Sure, you may say, but cannot we just use a better implementations able to complete the estimation? Aren't there any such methods?
OK, let's try with a different implementation!
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

Now let's increase the mean in one group by 0.01 and again use the more "robust" method:
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

Look at the above example again. Our model-based testing (via proportional-odds model) found N(0, 1) vs. N(5, 1) not statistically significantly different,
but when we changed the comparison to N(0.01, 1) vs. N(5, 1) - it "magically" worked well!

So, as I said - even if a model *actually converges*, the estimate may be unreliable. Always remember to check your outcomes.

**AN IMPORTANT LESSON:** Having multiple options, always choose the method that alerts you that something went wrong rather than method that silently pretends nothing wrong happened and happily continues. **It's always better to have NO result rather than having WRONG result.**

You. Have. Been. Warned.

PS: **Wait!** You tested totally different hypotheses with these methods, so maybe that's the issue?
**NO.** The ordinal logistic regression with just 2-level categorical predictor and no covariates is equivalent to Mann-Whitney (-Wilcoxon), which tests for stochastic dominance (_NO, it's NOT about medians_). Under IID it reduces to pseudo-median difference. Under symmetry of the distributions of [Walsh averages](https://stats.stackexchange.com/questions/215889/prove-the-relationship-between-walsh-averages-and-wilcoxon-signed-rank-test) it gives median difference.
And, under normality of both distributions, the median difference = mean difference. Which is equal to difference in means (in the normal distribution mean=median)
Following me?

``` r
set.seed(1000); x1 <- rnorm(100000000); x2 <- rnorm(100000000, mean=10)
sprintf("mean(diff)=%f, diff(means)=%f, median(diff)=%f, diff(medians)=%f, ps-median(diff)=%f",
mean(x1 - x2), mean(x1) - mean(x2), median(x1 - x2), median(x1) - median(x2), wilcox.test(x1, x2, conf.int = TRUE, exact = FALSE, adjust = FALSE)$estimate)
```

Now let's look from another perspective. The ordinal logistic regression naturally fits the Mann-Whitney (-Wilcoxon) null hypothesis. Citing the `rms` package from Prof. Harrell: "_[orm] fits ordinal cumulative probability models for continuous or ordinal response variables [...]. The ordinal cumulative probability models are stated in terms of exceedance probabilities (P rob[Y ≥ y|X]) so that as with OLS larger predicted values are associated with larger Y._"

So actually - with just different methods - I tested same hypotheses. Just indirectly.

## OK, if it's so problematic, then why even bother?
Well, classic tests are "simple" and fast. But simple method is for simple scenarios.
A more advanced inferential analysis often goes FAR beyond that these tests can do.

**Try doing THIS with classic tests!**

1. ordinary tests won't allow you to control for covariates. Goodbye more advanced analyses.
2. most classic tests cannot handle more complex m x n-way designs. Want to test some outcome across multiple groups rendered by `sex * visit * treatment_arm`? Forget! You will likely need to run multiple simple tests and they won't be able to detect inter-variable relationships. 
3. most of the classic tests won't be able to test interactions between multiple factors (a few modern ones, like ATS (ANOVA-Type Statistic) or WTS (Wald-Type Statistic) can do this, but only in a limited scope (1-level interaction between just 2 factors).
4. classic tests won't allow you to test simple effects via contrasts, e.g.: "Visit: 3, Arm: Treatment | Male vs Female effect" vs. "Visit 3, Arm: Control | Male vs Female effect". **For the model-based testing it's a piece of cake**.
5. you may simply NOT KNOW which test to use! Believe me or not, there are 850+ (EIGHT HUNDRED FIFTY!) statistical tests and counting. A colleague of mine has been counting them for years (with my little support). With a model - you don't care about all these "version", just provide the formula, set some parameters and test the hypotheses you need. Need a trest for trend? Just user ordinal factor for your time variable.
6. and you will obtain standard errors and confidence intervals for these comparisons!
7. Want to test some hypotheses jointly? Forget with classic tests!
8. Want to effectively adjust for multiple testing using parametric exact method employing the estimated effects and covariaces through the multivariate t distribution ("MVT")? This is far better than Bonferroni :-) But forget this flexibility when using plain tests!

Fair enough?
  
## Conclusions
So you can see with your own eyes, that model-based testing has LOTS of advantages. But sometimes you will need just a simple, classic test that runs FAST. Especially, if you have lots of tests to do under multiple imputation conditions, with lengthy data, and you are approaching a deadline :-)

So the choice depends really on your needs.  I only want to show you that this is doable and how well.

