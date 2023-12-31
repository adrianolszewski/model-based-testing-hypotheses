Testing hypotheses through statistical models opens a universe of new possibilities. Learn how to improve your daily work with this approach.

It's not a new idea. In my old Polish books (~2007-2012) in statistics ANOVA and t-test were mentioned as special cases of the general linear model. That was the first time I realized that every parametric test (and numerous non-parametric ones) are inferential procedures **applied on a top** of various models. Later I found this approach in other books too.

![obraz](https://github.com/adrianolszewski/model-based-testing-hypotheses/assets/95669100/bd74461d-f695-434d-b755-f481fadb89ff)

**Then it turned out, that also many non-parametric tests can be accurately replicated with "classic" parametric or semi-parametric models!**

Sometimes the relationship is easy to find (like ANOVA vs. assessment of the main and interaction effects of the general linear model = reduction of the residual variance when comparing nested models), and sometimes it's not that direct (like Wald's z-test for proportions or Breslow-Day test and the logistic regression). Sometimes the properties connecting non-paramertic tests and parametric models are really surprising and breathtaking (like the Mann-Whitney and Ordinal Logistic Regression)!

You can find excellent resources on the Internet, like [Common statistical tests are linear models (or: how to teach stats)](https://lindeloev.github.io/tests-as-linear/), but I wanted to summarize my own experience with model-based testing, used in my work (biostatistics in clinical trials) on daily basis.

So far I examined Wald's and Rao z test for 2+ proportions (+ ANOVA-like), Cochran-Mantel-Haenszel (CMH), Breslow-Day, Cochran-Armitage, McNemar, Cochran Q, Friedman, Mann-Whitney (-Wilcoxon), Kruskal-Wallis replicated with the logistic regression.

But this extends much more and beyond just classic testing! General and Generalized Linear Models fit via Generalized Least Square, Generalized Estimating Equations and through Generalized Linear Mixed Models, (Mixed-effect) Quantile regression, Cox proportional-hazard regression, Accelerated Failure Time, Tobit regression, Simplex regression and dozens of others models will be at your disposal to test hypotheses in a way you probably were never told.

This document is incomplete and dynamic. I will update it over time, so stay tuned.

You will find descriptions of subsequent tests in the [Test subdirectory](https://github.com/adrianolszewski/model-based-testing-hypotheses/tree/main/Tests). So far I started describing the [Mann-Whitney-Wilcoxon vs. Ordinal Logistic Regression](https://github.com/adrianolszewski/model-based-testing-hypotheses/blob/main/Tests/Mann-Whitney%20(-Wilcoxon).md). Still a lot is to be done: https://github.com/adrianolszewski/model-based-testing-hypotheses/blob/main/Tests/Various%20tests%20to%20describe%20one%20by%20one.md

## Just one question: WHY?
Well, classic tests are "simple" and fast. But simple method is for simple scenarios.
A more advanced inferential analysis often goes FAR beyond that these tests can do.

For our considerations it's important to say, that **by applying Likelihood Ratio, Rao's, or Wald's testing procedure to an appropriate model you will be (mostly) able to reproduce what the classic "plain" tests do but also go far beyond that** (find below a section dedicated to these 3 methods).

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

Let's summarize the key facts:

1. Ordinary tests won't allow you to control for covariates. Goodbye more advanced analyses. Models will do it naturally (born for that).
2. Most classic non-parametric tests cannot handle more complex m x n-way designs. Want to test some outcome across multiple groups rendered by `sex * visit * treatment_arm`? Forget! You will likely need to run multiple simple tests and they won't be able to detect inter-variable relationships. 
3. Most of the classic tests won't be able to test interactions between multiple factors (a few modern ones, like ATS (ANOVA-Type Statistic) or WTS (Wald-Type Statistic) can do this, but only in a limited scope (1-level interaction between just 2 factors).
4. Classic tests won't allow you to test and compare simple effects via contrasts, e.g.: "Visit: 3, Arm: Treatment | Male vs Female effect" vs. "Visit 3, Arm: Control | Male vs Female effect". **For the model-based testing it's a piece of cake**.
5. You may simply NOT KNOW which test to use! Believe me or not, there are 850+ (EIGHT HUNDRED FIFTY!) statistical tests and counting. A colleague of mine has been counting them for years (with my little support). With a model - you don't care about all these "version", just provide the formula, set some parameters and test the hypotheses you need. Need a trest for trend? Just user ordinal factor for your time variable.
6. You will obtain standard errors and confidence intervals for these comparisons.
7. Want to test some hypotheses jointly? Model-based approach followed by Wald's testing (sometimes also LRT) will do the job.
8. Want to effectively adjust for multiple testing using parametric exact method employing the estimated effects and covariace\s through the multivariate t distribution ("MVT")? This is far better than Bonferroni :-) But forget this flexibility when using plain tests!
9. Need to run type-2 or type-3 AN(C)OVA-type analysis for binary data? Or counts? Or anything else? Use appropriate model capable of handling such data and run the necessary testing on the top this model! Depending on what kind of model you use and what link transformation was used, the interpretation will be more complex, but - at the end of the day - you will have the way to assess the impact of your predictors on the response.
10. The assumptions for AN(C)OVA on your data fail? As above - fit appropriate model addressing these issues and follow it with approprite testing procedure!

Fair enough?

### A longer by-the-way about the ANOVA/ANCOVA-like analysis

The AN(C)OVA you learned from textbooks is a special kind of a much more general procedure: _the assessment of main and interaction effects of a model containing categorical predictor variables_. This can be done in two ways:
1. by joint testing of the model coefficients (Wald's approach)
2. by comparing nested models (Wilks's Likelihood Ratio testing), one of which has the predictor of interest and the other does not. This actually assesses the reduction in residual variance between the models. In the Generalized Linear Model the deviance is assessed instead.

For the general linear model (under the hood of the "textbook AN(C)OVA" and t-test) the two approaches are equivalent. The whole procedure reduces to just comparing means across the groups rendered by the levels of your categorical predictors (adjusted for numerical covariates, if present) through appropriate contrasts (functions of means). If you test contrasts related to a certain predictor jointly, you obtain the "global" test for the main effect of this predictor. It enhances also to interactions. And that's all. Exactly the same way you can obtain main and interaction effects of any model you want - only circumstances may differ (e.g. LRT is not available for non-likelihood models, so you need to stick with Wald's).

Or differently. Let's consider a simple model with just one 3-level categorical predictor variable. You can test the impact of this variable in 2 ways:
- the estimated model coefficients are shifts (differences) between the corresponding group means and the intercept-related mean, so you can test theset coefficiets jointly for the effect of the whole variable. That's the Wald's approach.
  
- if the variable of interest has statistically significant impact on response (means), then removing it from the model will worsen it (now "throwing" everything to theoveral mean - the intercept), so the residual variance will increase. You can test this change in variance by comparing the two nested model. This is the "omnibus" F-test you likely know from the ANOVA textbooks.

- In case of the GLM we replace the residual variance with deviance, directly related to the likelihood ratio. It's still based on model comparison - this time through the model likelihood ratio. If LR <> 1, then it means that removing the variable of interest worsened the fit. That's the Wilk's LRT aproach, employing the $\chi^2$ distribution.

**That's exactly how the R's anova(), car::Anova(), emmeans::joint_tests() and a few more methods work.**
Other statistical packages do it too.
See? **The model-based approach to ANOVA is not just "some observation made in some books", that's the common approach to it!**

PS: let's recall, by the way, that Chi2 is a limiting distribution for F under infinite denominator degrees of freedom:
If $X \sim F(n_1, n_2)$, the limiting distribution of $n_1X$ as $n_2 \rightarrow \infty$ is the $\chi^2$ distribution with $n_1$ degrees of freedom.
https://www.math.wm.edu/~leemis/chart/UDR/PDFs/FChisquare.pdf

You will find the $\chi^2$ distribution in at least two key places:
- when comparing models via LRT: under H0, as the sample size $n \rightarrow \infty$, the test statistic $\lambda _{\text{LR}}$ will be asymptotically $\chi^2$  distributed with degrees of freedom = to the difference in dimensionality of $\Theta$ - $\Theta_0$.

- when you perform Wald's testing on multiple parameters, then the distribution of the test statistic converges asymptotically with distribution to $\chi^2$.

You'll also meet the name _Wald's_ in 3 contexts, when the Wald's approach is applied to
- a single-parameter inference (simple H0), e.g. when doing simple comparisons (contrasts), and the degrees of freedom are known, you will see statistical packages reporting "Wald t". If you look at the formula, it looks much like a single-sample t-test (nominator is the difference, denominator - standard deviation), but that's a different creature and don't confuse the two. It's rather "t ratio". It's t-distributed only under normal sampling distribution of the parameter.
  
- a single-parameter inference, but the degrees of freedom aren't known (or are hard to guess) and assumed to be infinite. Then you will see statistical packages reporting "Wald z". You may find also "Wald chi2" - but that's the same, asymptotically $z^2 = \chi^2$.
  
- a multi-pararameter (like in ANOVA) asymptotic (under infinite degrees of freedom) inference employing composite H0. Then you will see statistical packages reporting "Wald chi2" or similar just "chi2". You may also find something like "Score chi2", which refers to Rao's approach to testing.

## What will you present in this repository?

We will do 3 things.

1. We will asesss how well the classic tests can be replicated with models. OK, I hear you: but why, if tests are faster nad rarerly fail? Because sometimes there may be no test that fits your needs in your statistical package.

2. We will extend the classic tests, for example to more than 1 factor with interactions. That's exactly the place where models enter the scene! Want to compare % of succeseses between males and females? Smokers and non-smokers? Check their interactions? In other words - do you want "ANOVA" for the binary response case? OK, that's exactly where we will go.

3. We will run some basic analyses of contrasts by employing severeal models and estimation methods and the awesome emmeans package. So, we will try GLM, GEE-GLM, GLS (MMRM) to compare % of success, counts, ratios, maybe concentrations (via log-linked gamma regression) via contrasts.

That's my plans. But first I need to cover the basic tests.

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
Read these excellent resources by Prof. Frank Harrell: [If you like the Wilcoxon test you must like the proportional odds model](https://www.fharrell.com/post/wpo/), [Equivalence of Wilcoxon Statistic and Proportional Odds Model](https://www.fharrell.com/post/powilcoxon/), and (best repository on this topic I've ever seen!) [Resources for Ordinal Regression Models](https://www.fharrell.com/post/rpo/) /

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

And sometimes it does not work at all and then you cannot proceed with this model. Try OLS with signed-ranked response instead.
Should work and be consistent with Mann-Whitney (-Wilcoxon).

----
## Wald's, Wilks' Likelihood Ratio, and Rao's testing. Why should I care?

Oh my, so many related topics! But nobody's said it will be a piece of cake.

I'm not going to write a textbook of statistical methods, so if you're curious about the Wilk's Likelihood Ratio (LR), Wald's and Rao testing, google these terms.
You can start from:
- [Likelihood-ratio test or z-test?](https://stats.stackexchange.com/questions/48206/likelihood-ratio-test-or-z-test)
- [FAQ: How are the likelihood ratio, Wald, and Lagrange multiplier (score) tests different and/or similar?](https://stats.oarc.ucla.edu/other/mult-pkg/faq/general/faqhow-are-the-likelihood-ratio-wald-and-lagrange-multiplier-score-tests-different-andor-similar/)
- [Wald vs likelihood ratio test](https://thestatsgeek.com/2014/02/08/wald-vs-likelihood-ratio-test/)
- [Score, Wald, and Likelihood Ratio](https://myweb.uiowa.edu/pbreheny/7110/f21/notes/10-25.pdf)

Cannot live wihout a solid dose of math? Look there:
- [Mathematical Statistics — Rigorous Derivations and Analysis of the Wald Test, Score Test, and Likelihood Ratio Test. Derivations of the Classic Trinity of Inferential Tests with Full Computational Simulation](https://towardsdatascience.com/mathematical-statistics-a-rigorous-derivation-and-analysis-of-the-wald-test-score-test-and-6262bed53c55)
- [Handouts for Lecture Stat 461-561 Wald, Rao and Likelihood Ratio Tests](https://www.cs.ubc.ca/~arnaud/stat461/lecture_stat461_WaldRaoLRtests_handouts_2008.pdf)
- [STAT 713 MATHEMATICAL STATISTICS II - Lecture Notes](https://people.stat.sc.edu/Tebbs/stat713/s18notes.pdf)

Briefly, if you imagine the IDEALIZED log-likelihood curve for some model parameter:
![obraz](https://github.com/adrianolszewski/model-based-testing-hypotheses/assets/95669100/d85738a7-f5a8-40d4-9463-5433c3f386d7)

then 
- **Wald's** is about the difference between the true value of some parameter p0 (e.g. the difference between groups) and the parameter p_hat estimated from data, placed on the horizontal axis (called the parameter space).
- **Wilk's LRT** is about the difference between two log-likelihoods (or ratio of two likelihoods, thus the name: "likelihood ratio testing") of a "full" model with the parameter of interest and a reduced model without it. This difference is placed on the vertical axis
- **Rao's score** is about the slope of the tangential line (Rao's score) to the log-likelihood curve.

OK, now let's make some notes:

1. Wald's, Wilks' and Rao's approach are just different ways to test hypotheses about some model parameter of interest. They offer 3 different ways to answer the question _whether constraining parameters of interest to zero (leaving out corresponding predictor variables from the model) reduces the fit of the model?_

2. In our considerations we'll primarily employ two statistical methods: Wald's test and Wilk's Likelihood Ratio Test (LRT). We'll closely examine how these methods behave in various scenarios.

3. Under the null hypothesis they provide asymptotically equivalent results. Asymptotically means "at infinite sample size". In real-world scenarios they will always differ, but unless you hit the "edge case", they will be rather consistent (at the end of the day, they assess a single thing). But still, Wald ≥ LRT ≥ Rao or, equivalently, $p-value_{WALD} \le p-value_{LRT} \le p_value_{RAO}$. Even if you don't notice that.
  
4. They may **noticeably** diverge in "edge cases", where the log-likelihood curve of a model parameter in the parameter space deviates from a parabolic shape. If you read either of the 3 first links, you will know what I mean.

5. Wald's test is extremely useful for testing hypotheses in the following cases:
   a. any kind of comparisons between group means and their combinations - via contrasts - **that's essential in experimental research employing planned comparisons**
   b. AN[C]OVA-type joint testing (of model coefficients) to assess the main and interaction effects in non-likelihood models, where it's the only option (__"why is it best? because it's the only available"), like GEE estimated GLM models (accounting for dependency and heteroscedasticity) or in quantile regression. So, despite it's critique (in small samples!) **it's important enough to be the widespread in all statistical packages and that's why you should know it**.

6. Wald's is also faster than LRT as it needs only a single model fitting. On the contrary, LRT needs at least 2 fits to compare 2 nested models: with and without the term you want to assess (you can also compare other settings, like covariance structures). With the increasing number of variables to assess, the number of compared models increases as well.

I hear you: _Yes, but computers nowadays are faster and faster, minimizing the importance of this very issue_. OK, but nature abhors a vacuum - whenever something drops, something else appears: 
- analyses employing **multiple imputation of missing data** (like MICE). In such setting the inference is repeated on each imputed dataset several times (20-50 and more) and then the results are pooled. The amount of necessary time may be not acceptable, especially if you have more analyses to complete.
- in small-sample inference, especially in repeated-observation (longitudinal) studies, **computationally-demanding adjustments to degrees of freedom** (like **Satterthwaite** and **Kenward-Roger** (DoF + covariance)) are in common use. They are standard in fields like mine (clinical trials). They can slow down the calculations several times!
- **employing complex covariance structures in dependent-observation model** (repeated, clustered) requires many parameters (co-variances) to be estimated, which increases the time of fitting noticeably!
- MLE, IWLS and similar methods are iterative, so let's hope the implementation you use is optimal, or it will slow down the LRT inference as well.

And now: ALL the four scenarios often **occur together**. In my work it's a daily reality to work with models accounting for repeated-observation, employinig small-sample adjustments to degrees of freedom, and working under the multiple imputation setting.

**LRT inference in this case can be really painful. If only the sample size > 15-20 observations, Wald's approach with all ciriticism can win with the complexity of LRT**.

9. LRT is found more conservative than Wald's approach in small samples, because of the relationship between the value of statistic obtained with these three approaches: Wald >= LR >= Rao ([Ranking of Wald, LR and score statistic in the normal linear regression model](https://stats.stackexchange.com/questions/449494/ranking-of-wald-lr-and-score-statistic-in-the-normal-linear-regression-model). It means, that it's less likely to reject a null hypothesis when it's true (i.e., LR has a lower Type-I error rate). This conservativeness can be advantageous when you want to avoid making false-positive errors, such as in clinical trials or other critical applications. In contrast, Wald's test can be more liberal, leading to a higher likelihood of false positives. That's why typically it's advised to select LRT over Wald's - as long as you are free to choose.

Let's summarize it with: "_When the sample size is small to moderate, the Wald test is the least reliable of the three tests. We should not trust it for such a small n as in this example (n = 10). Likelihood-ratio inference and score-test based inference are better in terms of actual error probabilities coming close to matching nominal levels. A marked divergence in the values of the three statistics indicates that the distribution of the ML estimator may be far from normality. In that case, small-sample methods are more appropriate than large-sample methods._ (Agresti, A. (2007). An introduction to categorical data analysis (2nd edition). John Wiley & Sons.)"

and

"_In conclusion, although the likelihood ratio approach has clear statistical advantages, computationally the Wald interval/test is far easier. In practice, provided the sample size is not too small, and the Wald intervals are constructed on an appropriate scale, **they will usually be reasonable** (hence their use in statistical software packages). In small samples however the likelihood ratio approach may be preferred._" [Wald vs likelihood ratio test](https://thestatsgeek.com/2014/02/08/wald-vs-likelihood-ratio-test/)

9. Wald's may not perform well when sample sizes are small or when the distribution of the parameter estimates deviates from normality. In other words, Wald tests assume that the log-likelihood curve (or surface) is quadratic. LRT does not, thus is more robust to the non-normality of the sampling distribution of the parameter of interest.

I quickly wrote a simulation using the general linear model with a single 3-level categorical predictor (just for fun), response sampled from the normal distribution, with parameters chosen empirically so that the p-values span several orders of magnitude under the inreasing sample size. The statements made in point 8 and 9 reflected by the results. Under normality of the parameter sampling (natural in this scenario), f**or N>10-20 the differences betweeh the methods are negligible so there's no need to blame Wald's and pray LRT**.

```r
set.seed(1000)
do.call(rbind,
        lapply(c(5, 10, 50, 100, 200), function(N) {
          print(N)
          do.call(rbind,
                  lapply(1:100, function(i) {
                      tmp_data <- data.frame(val = c(rnorm(N, mean=0.9), rnorm(N, mean=1), rnorm(N, mean=1.3)),
                                             gr = c(rep("A", N), rep("B", N), rep("C", N)))
                      m1 <- lm(val ~ gr, data=tmp_data)
                      m2 <- lm(val ~ 1, data=tmp_data)
                      data.frame(Wald = lmtest::waldtest(m1, m2)$`Pr(>F)`[2],
                                 LRT =  lmtest::lrtest  (m1, m2)$`Pr(>Chisq)`[2])
                      })) %>% cbind(N=N)
          })) %>% 
  mutate(N = factor(N), 
         "Wald to LRT" = Wald / LRT) -> results

library(gghalves)
library(patchwork)
  (results %>% 
    ggplot(aes(x=Wald, y=LRT, col=N)) +
    geom_point() +
    geom_line() +
    theme_bw() +
    scale_y_log10() +
    scale_x_log10() +
    geom_abline(slope=1, intercept = 0) +
    labs(title="Wilks' LRT vs. Wald's p-values; General Linear Model")) +
  (results %>% 
    ggplot(aes(x=N, y=`Wald to LRT`)) +
    geom_half_boxplot(nudge = 0.1, outlier.color = NA) +
    geom_half_violin(side = "r", trim = TRUE, nudge = .1, scale = "width") +
    geom_point() +
    scale_y_continuous(breaks=c(0,0.5,1,1.5,2,3,4,5, max(round(results$`Wald to LRT`))), trans="log1p") +
    geom_hline(yintercept=1, col="green") +
    theme_bw() +
    labs(title="Wald : LRT ratio of p-values; General Linear Model"))
```
![obraz](https://github.com/adrianolszewski/model-based-testing-hypotheses/assets/95669100/833e505d-2aa9-4a88-bc12-d193a5add1b9)

PS: Something weird here - Wald is MORE conservative in this scenario! P-values are larger than for LRT, which means lower test statistic. Need to investigate it!

If we move away from the quadrative log-likelihood (typically for binomial-distribution problem), we observe discrepancies, as anticipated by the theory. But when the differences matter, it's far below any common significance level. So again, Wald's not that bad here!

```r
set.seed(1000)
do.call(rbind,
        lapply(c(5, 10, 50, 100, 200), function(N) {
          print(N)
          do.call(rbind,
                  lapply(1:100, function(i) {
                    
                    tmp_data <- data.frame(val = c(sample(x = 0:1, size=N, prob=c(5,5), replace=TRUE),
                                                   sample(x = 0:1, size=N, prob=c(6,4), replace=TRUE), 
                                                   sample(x = 0:1, size=N, prob=c(5,5), replace=TRUE)),
                                           gr = c(rep("A", N), rep("B", N), rep("C", N)))
                    
                    m1 <- glm(val ~ gr, family = binomial(link="logit"), data=tmp_data)
                    m2 <- glm(val ~ 1, family = binomial(link="logit"), data=tmp_data)
                    
                    data.frame(Wald = lmtest::waldtest(m1, m2)$`Pr(>F)`[2],
                               LRT =  lmtest::lrtest  (m1, m2)$`Pr(>Chisq)`[2])
                    
                  })) %>% cbind(N=N)
        })) %>% 
  mutate(N = factor(N), 
         "Wald to LRT" = Wald / LRT,
         Magnitude = abs(Wald - LRT)) -> results

results %>% 
  group_by(N) %>% 
  summarize(median_magn = median(Magnitude)) -> diff_magnit

(results %>% 
    ggplot(aes(x=Wald, y=LRT, col=N)) +
    geom_point() +
    geom_line() +
    theme_bw() +
    scale_y_log10() +
    scale_x_log10() +
    geom_abline(slope=1, intercept = 0) +
    labs(title="Wilks' LRT vs. Wald's p-values; Generalized Linear Model")) +
  (results %>% 
     ggplot(aes(x=N, y=`Wald to LRT`)) +
     geom_half_boxplot(nudge = 0.1, outlier.color = NA) +
     geom_half_violin(side = "r", trim = TRUE, nudge = .1, scale = "width") +
     geom_point() +
     scale_y_continuous(breaks=c(0,0.5,1,1.5,2,3,4,5, max(round(results$`Wald to LRT`))), trans="log1p") +
     geom_hline(yintercept=1, col="green") +
     theme_bw() +
     labs(title="Wald : LRT ratio of p-values; Generalized Linear Model") +
     geom_text(data=diff_magnit, aes(y=-0.1, x=N, label=sprintf("Me.Δ=%.4f", median_magn)), size=3))
```

![obraz](https://github.com/adrianolszewski/model-based-testing-hypotheses/assets/95669100/c4a39713-9b7e-4952-9a32-e98bdfe08217)
  
10. Sometimes Wald's testing fails to calculate (e.g. estimation of covariance fails), while the likelihood is still obtainable and then the LRT is the only method that works. Who said the world is simple? And sometimes the LRT is not available, as mentioned above. Happy those, who have both at their disposal.

11. LRT allows one for performing AN[C]OVA-type analyse (which requires careful specification of the model terms!) but mat not be available in your software for a more complex analysis of contrasts (e.g. in planned covariate-adjusted comparisons), so you'll have to specify the nested models by hand. At the same time, Wald's approach takes full advantage of the estimated parameter and covariance matrix, **which means that "sky is the only limit" when testing contrasts of any complexity**. With LRT it will be harder to obtain appropriate parametrization, but it's still doable if only you carefully set contrasts! (TODO: show example with `glmglrt`)
    
12. **Wald's (and Rao!) inference is not transformation (or re-parametrization) invariant**, using derivatives, sensitive to it. If you use Wald's to calculate p-value or confidence interval on two different scales, e.g. probability and logit transformed back to the probability scale, you will get different results. Often they are consistent, but discrepancies may occur at the boundary of significance and then you're in trouble. By the way, Wald's on the logit scale will return sensible results, while Wald's applied to probability scale may yield negative probabilities (e.g. -0.12) or exceeding 100% (e.g. 1.14). This is very important when employing LS-means (EM-means) on probability-regrid scale(!).
On the contrary, the **LRT IS transformation invariant** and will have exactly the same value regardless of any monotononus transformation of the parameter.

/ Side note: Wald's assumes normally distributed parameter sampling, which briefly means **symmetry**. 0-1 truncated data will never be so, that's why you may easily obtain negative or >1 bounds of the confidence interval. See? I believe that's the manifestation of _Hauck–Donner effect in binomial models_ when the estimated parameter is close to the boundary of the parameter space (here close to 0 or 1 and at a very small sample size). 
```r
> binom::binom.asymp(x=1, n=5)
      method x n mean     lower    upper
1 asymptotic 1 5  0.2 -0.150609 0.550609   # -0.15 ?
> binom::binom.asymp(x=4, n=5)
      method x n mean    lower    upper
1 asymptotic 4 5  0.8 0.449391 1.150609    # 1.15?
```
/

13. There are another problematic areas for Wald's testing, like mixed models (ironically, numerous statistical packages offer only (or mostly) such inference due to big complexity of the solutions!). For example read: [Paper 5088 -2020; A Warning about Wald Tests](https://support.sas.com/resources/papers/proceedings20/5088-2020.pdf)

## You said that that different method of testing (Wald's, Rao's, Wilk's LRT) may yield a bit different results?

**Yes. But serious differences mostly happens under the ALTERNATIVE hypothesis and at low p-values.**
Let's consider a 2-sample test comparing means (or any other quantity of interest).
A corresponding (general linear with treatment contrast coding) model will have just the response (dependent) variable and a single, 2-level categorical predictor variable which levels form the groups we want to compare: `response ~ group`.

The fitted model will have two estimated parameters: the intercept coefficient (representing the 1st group mean, at the reference level of the group indicator) and the coefficient representing the shift (in means or anything) between the 2nd group and the intercept. It's value added to the intercept term will give us the 2nd group mean.

Now, imagine the log-likelihood curve for the "shift" parameter. Under the null hypothesis the difference between the means is 0, so is the parameter. This generalized directly to comparing multiple groups through a predictor variable with multiple levels (forming the groups). 

Now, if removing this predictor from the model (effectively - if the parameter is ~zero) doesn't change the fit, than - well - this predictor is "useless", as  it didn't affect the response. Which means that the groups distinguished by this predictor levels are mutually similar: respective "shifts" are zero, all means are similar enough to say they come from the same population. Which means that H0 cannot be rejected.

![obraz](https://github.com/adrianolszewski/model-based-testing-hypotheses/assets/95669100/8c3bd8ed-df10-4da1-bfde-d96574c24da5)

If you operate under H0, then regardless of its shape (ideally quadratic, but here let's agree on a concave one, with single maximum):
- **Wald's:** the difference between groups = 0, so the distance in the parameter space approaches 0.
- **Wilk's LRT** the difference in both log-likelihoods between both models approaches 0 (then the likelihood ratio approaches 1) - no gain.
- **Rao's score** the slope of the tangential line (Rao's score) approaches zero as well. Which reflects the fact that under H0 the derivative (score) of the log-likelihood function with respect to the parameter should approach zero.

The word "approaches 0" means, that we operate in its neighbourhood, where (in case of concave shape) changes in all 3 measures are relatively small and close to each other (they follow each other). This explains why mostly you will obtain similar results under the nul hypothesis, and why they start diverging as you go away from this "safe area" towards the alternative hypothesis.

This is a simple graphical explanation, but formal proofs can be found in many places, e.g.: [Score Statistics for Current Status Data:
Comparisons with Likelihood Ratio and Wald Statistics](https://sites.stat.washington.edu/jaw/JAW-papers/NR/jaw-banerjee-05IJB.pdf) or [Mathematical Statistics — Rigorous Derivations and Analysis of the Wald Test, Score Test, and Likelihood Ratio Test](https://towardsdatascience.com/mathematical-statistics-a-rigorous-derivation-and-analysis-of-the-wald-test-score-test-and-6262bed53c55)

In simpler words, under "_no difference_" we are approaching the maximum of this curve which must be reflected by either perspective: "horizontal" (Wald), "vertical" (Wilks LRT) and "slope" (Rao).

It's ideally observed if the distribution of parameter is normal, which is necessary for the Wald's approach to work properly.
But remember that normality never holds exactly, as it's only an idealized, theoretical construct! And because the normality is always approximate (at finite sample size), the Wald's and LRT p-values may be **very close** to each other, but not exactly (yet you may not be able to distinguish them in the "good case").

/ BTW: That's why you saw, many times, the "t" (if degrees of freedom are easy/possible to derive) or "z" (otherwise; at infinite degrees of freedom) statistic reported for coefficients of most regression models and when testing contrasts. That's because the sampling distribution for these parameters is theoretically assumed to be normal (via Central Limit Theorem). But it's not always the case. For less "obvious" models, like quantile regression, such approximation may be incorrect (or limited only to "good cases"), and then resampling techniques are typically employed (e.g. bootstrap, jacknife). /

**But the more the log-likehood curve deviates from the quadratic shape, the more they will diverge from each other even under the null hypothsesis!**

And the shape may diverge from the desired one also by using common transformations, e.g. from logit scale to probabiltiy scale, when employing the logistic regression family.

The mentioned, excellent article [Wald vs likelihood ratio test](https://thestatsgeek.com/2014/02/08/wald-vs-likelihood-ratio-test/) by Prof. Jonathan Bartlett shows additional examples that speak louder than words!

![obraz](https://github.com/adrianolszewski/model-based-testing-hypotheses/assets/95669100/4b89d9d8-c571-4ddf-a137-8a3f1f8889e0)
(cited with Prof. Bartlett's permission).

**But what will happen under alternative hypothesis?**
Briefly - the 3 measures will increase. Whether this will happen in a more or less consistent manner (you remember the relationship: Wald ≥ LRT ≥ Rao) and depends on the log-likelihood curve. You may observe a good agreement in the beginning (near to H0) and then quick divergence.

**The good news is that the discrepancies will rather occur at small p-values, typically below common levels of statistical significance level, so for us, users of the test, the differences may be barely noticeable**, if we are lucky. But read the next section...

Let me show you an example. We will compare data sampled from two normal distributions in four settings:
1. Both groups are sampled from the same N(50, 20). We operate under the null hypothesis. Practically full overlap of empirical distributions.
```r
set.seed(1000)
x1 <- rnorm(100, 50, 20)
x2 <- rnorm(100, 50, 20)
effectsize::p_overlap(x1, x2)
Overlap |       95% CI
----------------------
0.96    | [0.86, 1.00]
```

2. One group sampled to have a smaller mean, with variances unchanged: N(50, 20) vs. N(35, 20). We operate under alternative hypothesis, but the data overlap mostly (this is important for the ordinal logistic regression, as the log-likelihood data should be more or less symmetric)
```r
# Exemplary sampling
set.seed(1000)
x1 <- rnorm(100, 50, 20)
x2 <- rnorm(100, 35, 20)

effectsize::p_overlap(x1, x2)
Overlap |       95% CI
----------------------
0.73    | [0.62, 0.84]

> tidy(t.test(x1, x2, var.equal = TRUE))
# A tibble: 1 × 10
  estimate estimate1 estimate2 statistic    p.value parameter conf.low conf.high method            alternative
     <dbl>     <dbl>     <dbl>     <dbl>      <dbl>     <dbl>    <dbl>     <dbl> <chr>             <chr>      
1     13.3      50.3      37.0      4.92 0.00000179       198     7.98      18.7 Two Sample t-test two.sided  
```

3. Like in case 2, but now variance is reduced: N(50, 10) vs. N(35, 10). Now the data are less overlapping (bigger separation). This will distort the log-likelihood curve.
```r
# Exemplary sampling
set.seed(1000)
x1 <- rnorm(100, 50, 10)
x2 <- rnorm(100, 35, 10)
effectsize::p_overlap(x1, x2)
Overlap |       95% CI
----------------------
0.46    | [0.37, 0.56]
```

4. Like in case 3, with further reduced variance: N(50, 5) vs. N(35, 5).
```r
set.seed(1000)
x1 <- rnorm(100, 50, 5)
x2 <- rnorm(100, 35, 5)
effectsize::p_overlap(x1, x2)
Overlap |       95% CI
----------------------
0.13    | [0.08, 0.19]
```

And now let's have a look at the p-values obtained from the ordinal logistic regression (as if we applied the Mann-Whitney test to normally distributed data) by
employing the three metehods:

```r
compare_3methods_distr(n_group = 50, samples=50,
                       distr_pair_spec = list("N(50,20) vs. N(50,20)" =
                                                c(arm_1_distr = "rnorm(n_group, mean=50, sd=20)",
                                                  arm_2_distr = "rnorm(n_group, mean=50, sd=20)"),
                                              "N(50,20) vs. N(35,20)" =
                                                c(arm_1_distr = "rnorm(n_group, mean=50, sd=20)",
                                                  arm_2_distr = "rnorm(n_group, mean=35, sd=20)",
                                                  log = "TRUE"),
                                              "N(50,10) vs. N(35,10)" =
                                                c(arm_1_distr = "rnorm(n_group, mean=50, sd=10)",
                                                  arm_2_distr = "rnorm(n_group, mean=35, sd=10)",
                                                  log = "TRUE"),
                                              "N(50,5) vs. N(35,5)" =
                                                c(arm_1_distr = "rnorm(n_group, mean=50, sd=5)",
                                                  arm_2_distr = "rnorm(n_group, mean=35, sd=5)",
                                                  log = "TRUE")))
```
![obraz](https://github.com/adrianolszewski/model-based-testing-hypotheses/assets/95669100/8c684e6e-cb79-4c79-80db-0579044a36e6)

As expected, the methods gave very similar p-values (look especially at the 1st, top-left plot for H0). Even, if differences occurred, they were too small to be observed. That's because the p-values spanned several orders of magnitude! And we could end here: the situation is perfect!
1. ideal agreement under H0 (where we are worrying the most for wrong rejections; type-1 error)
2. potencial discrepancies under H1 at so low magnitudes, so we can totally ignore it (in this case).

Let's log-transform the axes to have a better view:
OK, now we can observe the differences. No surprise they were barely noticeable at these magnitudes!
Generally we can observe, that the lower p-value, the bigger disrepancy. But why?

That's simple - because we more farther and farther from H0. With the convex curve, at some point, far from the maximum, small changes in argument (X) will result in big changes of value (Y). Here I have no idea what's the shape of the log-likelihood curve, which can even amplify the problem. So let's ignore it.

![obraz](https://github.com/adrianolszewski/model-based-testing-hypotheses/assets/95669100/62977972-f69e-49ba-9431-d331edfda0f4)

By this occasion, let's have a lookg how the p-values differ with respect to the Mann-Whitney (-Wilcoxon) test.
It's interesting to observe the opposite directions in which the p-values diverge.

/ PS: Remember, that the way the p-value is calculated in Wilcoxon depends on the implementation! There can be a normal approximation, exact method, ties may be unhandled or handled via different methods (Wilcoxon, Pratt). /

```r
compare_3methods_distr_model_vs_test(n_group = 50, samples=50,
                        distr_pair_spec = list("N(50,20) vs. N(50,20)" =
                                                 c(arm_1_distr = "rnorm(n_group, mean=50, sd=20)",
                                                   arm_2_distr = "rnorm(n_group, mean=50, sd=20)"),
                                               "N(50,20) vs. N(35,20)" =
                                                 c(arm_1_distr = "rnorm(n_group, mean=50, sd=20)",
                                                   arm_2_distr = "rnorm(n_group, mean=35, sd=20)",
                                                   log = "TRUE"),
                                               "N(50,10) vs. N(35,10)" =
                                                 c(arm_1_distr = "rnorm(n_group, mean=50, sd=10)",
                                                   arm_2_distr = "rnorm(n_group, mean=35, sd=10)",
                                                   log = "TRUE"),
                                               "N(50,5) vs. N(35,5)" =
                                                 c(arm_1_distr = "rnorm(n_group, mean=50, sd=5)",
                                                   arm_2_distr = "rnorm(n_group, mean=35, sd=5)",
                                                   log = "TRUE")))
```
![obraz](https://github.com/adrianolszewski/model-based-testing-hypotheses/assets/95669100/51718316-4a2e-41d6-a951-249617404a23)

## Even worse! You just said that different method of testing can yield OPPOSITE results!

Well, yes, Wald's and LRT **may yield opposite results**, e.g. Wald's p-value = 0.8 vs. LRT p-value < 0.001.

But that doesn't happen very often, in "**edge cases**". The problem is that you may NOT know when it's an "edge case".
That's why doing statistical analysis is not just a brainless "click-and-go" task!

If only you suspect something is wrong, I mean - that the results of some analysis does not reflect what you observe in data, validate it with another analysis.

This is summarized in the article I already cited ([Wald vs likelihood ratio test](https://thestatsgeek.com/2014/02/08/wald-vs-likelihood-ratio-test/)): "Further, a situation in which the Wald approach completely fails while the likelihood ratio approach is still (often) reasonable is when testing whether a parameter lies on the boundary of its parameter space." 

If you anticipate that something bad may happen, VALIDATE your calculations using a different method. It doesn't have to bring SAME results (it's a different method) but at least you will be able to assess if they are approximately consistent.

If, for instance, you run the ordinal logistic regression over Likert data and obtain Wald's p-value = 0.3 and LRT p-value = 0.0003, then you can also simplify it for a while and check with the Mann-Whitney (-Wilcoxon) test or even run the OLS with ranked response (it's a very good approximation).

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

**You. Have. Been. Warned.**

## Any words of comfort...?

Yes :)

First, for several tests, to reproduce it closely with a model, you just need one approach and don't have to overthink the problem.
A few examples of perfect agreement you can find in this gist (temporarily; I will move it to this repo soon): https://github.com/adrianolszewski/Logistic-regression-is-regression/blob/main/Testing%20hypotheses%20about%20proportions%20using%20logistic%20regression.md

If you use R:
- the prop.test() (z test with pooled standard errors) will exactly match Rao ANOVA over logistic regression
- the classic z test for 2 proportions with non-pooled standard errors will exactly match Wald's test over AME over logistic regression
- Cochrane-Mantel-Haenszel will exactly match Wald's test applied to the conditional logistic regression
- McNemar and Cochrane Q will be quite close with GEE-estimated logistic regression (mostly with exchangeable covariance) tested with Wald's
- Breslow-Day will agree with Rao ANOVA over the logistic regression.

And then you don't have to to split hairs. Remember also, that as long, as your model contains multiple categorical predictors and numerical covariates, you cannot compare it with most classic tests (limited to just 1 categorical predictor).

## I give up! Too many information! Briefly - is Wald just bad and LRT just awesome?

No. At small samples (N~10-30) LRT will be less biased than Wald's. Above they will do approximately equally well (unless edge cases and problematic transformations). Wald's will be eaiser to work with and much faster. That's it.

More precisely:

**LRT - advantages:**
1) At small samples you rather won't test complex hypotheses and contrasts, but rather the effect off (a) single variable(s), e.g. "treatment arm". LRT for such comparison will be avaiable in your statistical package and will be likely better than Wald's
2) Whatever scale you choose (e.g. in binomial models: log-ods or probability) - LRT will do equally well (unlike Wald's)
3) In simple scenarios you may not notice the increased testing time due to multiple fitting.
4) May be available even when Wald's fails to calculate.
5) It's possible to analyze complex contrasts (including difference-in-difference) with LRT as well as with Wald's. It won't be as simple as with the `emmeans` package, and you will need to construct the contrasts carefully selecting appropriate model coefficients.

LRT - disadvantages or no clear benefit:
0) When you work with small N, like N<20, then ANY analysis you do may be unreliable. If you believe that "tools and tricks" will manage, you cheat yourself. Think twice, then pick the LRT.
1) If you need to test contrasts accounting for appropriate covariance(s) and degrees of freedom adjustments, like in repeated-data setting, the LRT rather won't be available in your stat. software.
2) For testing more complex contrasts it won'y be as easy as with, say, the `emmeans` package, where you have all LS-means (EM-means) listed nicely, so just pick the right ones. Model coefficients will depend on the parametrization and for the default, `treatment contrast` (dummy) you will need to be very careful about the reference levels represented by the model coefficients! Also, only a very few models are supported in R package (mostly GLM, Cox), but GLS (MMRM) require extra work...
4) If you need to test something over a non-likelihood model (e.g. GEE-estimated) - LRT will no be available.
5) If you work under a computational-demanding setting (multiple imputation + Kenward-Roger + unstructured covariance + multiple variables, you may prefer not to use it due to unacceptably long time of calculations.
6) Support in statistical packages may be seriously limited. Lots' of necessary work may be up to you.
7) No test for LS-means (EM-means). Recall, testing model coefficients is NOT the same as testing LS-means in general. Moreover, you may want to test at CERTAIN value of the numerical covariate.
   
Wald's - advantages:
1) As fast as fitting the model sigle time. Matters when run under multiple imputation.
2) OK if you pick apporpriate scale (e.g. log-odds rather than probability; Note, the probability scale CAN BE USED TOO, just with care and confirm it with the analaysis on log-odds)
3) In models assuming conditional normality and at N>20-50 (depending on model) - not much worse than LRT. I guess you have enough data, if you start playing with complex contrasts?
4) The farther from the null hypothesis, the less you care about the discrepancies (occurring typically at low p-values). Who cares if Wald p=0.0000000034 and LRT p=0.000000084? Don't make it absurd!
5) You can compare any contrasts you want, accounting for covariances, degrees of freedom small-sample adjustments, robust estimators of covariance, any kind of covariance structure, etc.
6) rellted with 5 - so you can use the multivariate-t distribution ("mvt") based exacty adjustment for multiple comparisons
7) Available for non-likelihood models (GEE-estimated GLM; and such models usually NEED more data!)
8) Available in all statistical packages.
9) Allows for testing LS-means (EM-means) and at CERTAIN values of the numerical covariate(s).

Wald's - disadvantages:
1) At N<20-50 may be anti-conservative, worse than LRT
2) May fail to calculate
3) May give different estimates when used on different scale (e.g. log-odds vs. probability). You can calculate it in both and - if consistent - report the one you prefer. Otherwise report the one you find more relevant.
4) In edge case can yield very different results - typically that's a sign of small sample or very "edge case". In such case don't trust either, but prefer LRT.

There is no one-size-fits-all solution in statistics.

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

**So the choice depends really on your needs.  I only want to show you that this is doable and how well it performs.**
