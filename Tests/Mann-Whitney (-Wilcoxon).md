# Mann-Whitney (-Wilcoxon)

## Table of contents
1. Description of the test
2. Model-based implementation via Ordinal Logistic Regression (aka Proportional-Odds model) - categorical ordinal data
3. Model-based implementation via Ordinal Logistic Regression (aka Proportional-Odds model) - numerical data
4. References
5. [Simulations - how well the model-based approach replicates the classic test](#simulations)

----
<a name="description"></a>
## Description of the test
This is a non-parametric test for testing H0: stochastic equality vs. H1: stochastic superiority (aka dominance).

Surprised, that I did not wrote "for comparing medians"? Well, I did NOT, because it does NOT - in general, unless additional, strong distributional assumptions are met:
1. **data in both compared samples are IID**. It means: have same dispersion (scale parameter, more/less variance) and similar shape (e.g. if skewed, then in same direction). Think about it. Stochastic superiority can be imposed by 1) shift in locations, by 2) unequal variances and by 3) differences in shapes of the distributin. It's like saying that Mann-Whitney (-Wilcoxon) test is sensitive to ALL THREE factors. If you want to test for the shift in location parameter, you have to address the two other differences. In other words, if both dipersions are similar, if both shapres are similar as well, then what remains is the possible difference in locations, right? 

And if your data meet this assumption, then Mann-Whitney gives you the _Hodges–Lehmann estimator of pseudo-median of differences_ between the compared groups.
WHAT!? _Pseudo-median difference_? Not "difference in medians"?

No. But **if the distribution of differences betwene the groups is symmetric, then it reduces to median difference.**

Quick and dirty check (you can use bigger numbers for better agreement, but it may take whole lot of time!)
```r
# set.seed(1000); 
# x1 <- rnorm(1000000); 
# x2 <- rnorm(1000000, mean=10)
# sprintf("mean(diff)=%f, diff(means)=%f, median(diff)=%f, diff(medians)=%f, ps-median(diff)=%f",
# mean(x1 - x2), mean(x1) - mean(x2), median(x1 - x2), median(x1) - median(x2), 
 wilcox.test(x1, x2, conf.int = TRUE, exact = FALSE, adjust = FALSE)$estimate)

[1] "mean(diff)=-9.999591, diff(means)=-9.999591, median(diff)=-10.000615, diff(medians)=-9.998330, ps-median(diff)=-9.999202"
```

2. **data in both compared samples are symmetric around their medians**. Only now it reduces to difference in medians.

In my other respository [Mann-Whitney fails as a test of medians in general](https://github.com/adrianolszewski/Mann-Whitney-is-not-about-medians-in-general) you will find numerous resources (mostly free and some available via Google Books) to learn about this problem and simulations proving that treating this test as a test of medians can really be dangerous.

If your teacher, statistical reviewer or boss will "force" you to use Mann-Whitney (-Wilcoxon) as a test of medians withut checking the 2 assumptions, show them the literature I recommend (full of formulas and examples) and the results of my simulations. If this does not convinvce them, quit or resign (if you can) or at least put "On teacher's/client's demand" somewhere in your work, to make yourself less guilty of harming research... Yes, seriously - it's easy to reject the null hypothesis under exactly equal medians with this test (which invalidates it as a test of medians) - and type-1 error rate can go far beyond 5%, 10% and even 20%.

**Take-home message:**
If you want to compare quantiles (including medians), use the quantile regression (will show you how elsewhere), the Brown-Mood test of medians or make sure the above assumptions hold in your case. And always visually assess your data.

BTW: you will find articles, that "Mood test of medians should be ditched". And the authors say? That this test lacks power, so they propose Mann-Whitney, because it has greater power. But they don't say A WORD about the assumptions. So indeed, Mann-Whitney (-Wilcoxon) is MORE powerful, but at the cost of being sensitive to not only locations, but also dispersions and shapes.

If you are fine with it, if it's your actual goal - then OK, Mann-Whitney (-Wilcoxon) will do its job. But only then. If you want to compare medians, as I wrote, choose quantile regression (ideally). It is designed to deal with quantiles, so no better method will ever exist.

----
<a name="simulations"></a>
## Simulations
Let's check how good is the approximation by performing some simulations. Nothing will give us better feeling of situaion than that.

### Methodoloy
I set up a loop, in which data were sampled (according to some patterns) and the p-values from both ordinary test and its model counterpart were collected.

Then I display these p-values in a series of graphs:
1. Test vs. Model to observe how well they follow each other and observe discrepancies
2. Ratios of p-values - the closer to 1 the better (this works well regardless of the magnitude of p-values)
3. Raw differences in p-values. These differences do depend on the magnitude of p-values, and that's exactly what I wanted to see. At very small p-values (p<0.001) I would not care about differences like p=0.00032 vs. p=0.00017 because both results are far from the typical signfiicance levels (0.1, 0.01, 0.05, 0.01 and also 0.001).
For practical purposes I would care if they were discrepant near 0.05.
4. Comparison of % of cases where Test > Model vs. Test < Model, with additional descriptie statistics of the p-values. Because raw differences depend on the magnitude, it's difficult to observe them on a common scale. Log transformation doesn't help there (0, negative values).
So this approach gives me some feelings about the situation
5. How many times both methods disagreed on the rejection of the null hypothesis and what were the actual p-values.

I repeated the simulations sample sizes:
- n_group = 20 - typical small data size; smaller samples don't make much sense in a reliable research,
- n_group = 30 - the "magical number" so many statisticians find "magical". Of course it's not, but this constitutes a reasonable "lower limit of data"
- n_group = 50 - still small data but in a safer area
- n_group = 100 - safe area (optionally)

Actually, the type-1 error and power is totally off topic in this simulation, we only check how well the methods follow each other, but having the opportunity - we will look at it too. We only have to remember that this is a long-run property and at 100 replications it will be "just close" rather than "nominal", especiallly at smaller samples. It's normal to have, say 7/100 (=0.07), 12/200 (=0.06) to finally reach 15/300 (=0.05).

Whenever possible I will provide both Wald's and LRT results.

---
### Simulations!
#### The data
For the ordinal logistic regression I made a dataset consisting of patients' responses to question about the intensity of physical pain they feel, part of the ODI (Oswestry Disability Index). Patients are split into two groups: those taking a new investigated medicine vs. those taking active control (existing standard of care).
So the model is `ODIPain ~ Arm` (arm stands for the treatment arm).

The response was recorded on a 6-level Likert ordinal item (don't confuse with Likert scale, where responses to items are summed together giving a numerical outcome).
I programmed it so in one group dominate low responses and in the other group - high responses.
Example:
```r
set.seed(1000)
lapply(1:100, function(i) {
  stack(
    data.frame(arm1 = sample(0:5, size=100, replace=TRUE, prob = c(20, 10, 5, 2, 2, 2)),
               arm2 = sample(0:5, size=100, replace=TRUE, prob = rev(c(20, 10, 5, 2, 2, 2))))) %>% 
    mutate(values_ord = ordered(values), ind = factor(ind))
}) -> data
```
Although the data size will vary during the simulations, I will visualize only this one case to save space. The stacked bar chart clearly shows that in one group dominate "dark" elements, while in the other group - "bright" ones. The density plots confirm the opposite nature of responses across both arms.

Our comparisons will get more power at increasing sample size, which is an ideal situation as we want to observe both high and low p-values.

``` r
lapply(1:100, function(i) {
  print(i)
  data[[i]] %>% 
    ggplot() +
    geom_bar(aes(x=ind, group=values, fill=values), position = "fill", show.legend = FALSE) +
    theme_void()
}) -> plots1

lapply(1:100, function(i) {
  print(i)
  data[[i]] %>% 
    ggplot() +
    geom_density(aes(x=values, group=ind, fill=ind, col=ind), show.legend = FALSE, alpha=0.6, adjust = 1.5)+
    xlab(NULL) +ylab(NULL) +
    theme_void()
}) -> plots2

(patchwork::wrap_plots(plots1, ncol = 10, nrow = 10) |
    patchwork::wrap_plots(plots2, ncol = 10, nrow = 10)) +
  patchwork::plot_annotation(title = "Patient-reported ODI pain scores across both treatment arms")
```
![obraz](https://github.com/adrianolszewski/Logistic-regression-is-regression/assets/95669100/a75158cb-c47e-44d8-8e33-acb98c6534a2)

#### The comparison engine
``` r
simulate_wilcox_olr <- function(samples, n_group, set, arm_1_prob, arm_2_prob, which_p) {
  set.seed(1000)
  
  mixed_hyp <- is.null(arm_1_prob) | is.null(arm_2_prob)
  
  data.frame( 
    do.call( 
      rbind, 
      lapply(1:samples, function(i) {
        print(i)
        
        if(mixed_hyp) {
          arm_1_prob <- sample(x = 1:10, size = 6, replace = FALSE)
          arm_2_prob <- sample(x = 1:10, size = 6, replace = FALSE)
        }
        
        stack(
          data.frame(arm1 = sample(set, size=n_group, replace=TRUE, prob = arm_1_prob),
                     arm2 = sample(set, size=n_group, replace=TRUE, prob = arm_2_prob))) %>% 
          mutate(values_ord = ordered(values), ind = factor(ind)) -> data
        
        #m <- joint_tests(MASS::polr(values_ord ~ ind , data = data, Hess=TRUE))$p.value)
        m <- rms::orm(values_ord ~ ind , data = data)
  
        if(!m$fail)
          case_when(which_p == "Wald" ~ as.numeric(2*pnorm(abs(m$coefficients["ind=arm2"] / sqrt(m$var["ind=arm2","ind=arm2"])), lower.tail = FALSE)),
                    which_p == "LRT" ~ as.numeric(m$stats["P"]),
                    which_p == "Rao" ~ as.numeric(m$stats["Score P"])) -> Model_p
        else
            Model_p <- NA
        
        c(iter = i,
          n_group = n_group,
          Test  = tidy(wilcox.test(values~ind, data=data, exact = FALSE, adjust = FALSE))$p.value,
          Model = Model_p)
      }))) -> result
  
  attr(result, "properties") <- list(arm_1_prob = arm_1_prob, 
                                     arm_2_prob = arm_2_prob, 
                                     under_null = ifelse(mixed_hyp, NA, identical(arm_1_prob, arm_2_prob)),
                                     samples = samples,
                                     n_group = n_group,
                                     title = "",
                                     method = which_p)
  return(result)
}
```
#### Results
##### Under H0 - same probabilities of obtating each score in both groups = same ordering of observations
###### 20 observations per group
``` r
simulate_wilcox_olr(samples = 100, n_group = 20, set = 0:5, 
                    arm_1_prob =  c(20, 10, 5, 2, 2, 2),
                    arm_2_prob =  c(20, 10, 5, 2, 2, 2),
                    which_p = "Wald") %>%  
  plot_differences_between_methods()
```
![obraz](https://github.com/adrianolszewski/model-based-testing-hypotheses/assets/95669100/e072f1a3-e3d3-48bb-b1ce-b754ab48a77c)

OK, let's make some notes:
1. look at the most left plot. The p-values from both methods are well aligned to the line of slope = 1, though some bias is visible towards Test p-values.
2. It means, that test was more conservative than model (higher p-values).
3. The fact, that practically all p-values coming a test were larger from the p-values obtained from the model is confirmed also by the bar plot (central-bottom).
4. When we look at the area below the 0.05 significance level, we notice 5-6 observations. Remembering we "work under the null" it means, that both methods exceeded the 5% significance level, reaching 9% in this simulation. At this sample size it's fine for exploratory analyses.
5. There was just single "contradiction" at very small difference of p-values 0.049 vs 0.051. 

With increased group size the situation will get only better. The type-1 error is well maintained in all cases (remembering what I said about sm

###### 30 observations per group
![obraz](https://github.com/adrianolszewski/model-based-testing-hypotheses/assets/95669100/0643dff4-a2f8-4ac0-9988-809f25b698e0)

###### 50 observations per group
![obraz](https://github.com/adrianolszewski/model-based-testing-hypotheses/assets/95669100/dd8f83d8-f804-4b8f-89a8-33260ad20c1d)

##### 100 observations per group
![obraz](https://github.com/adrianolszewski/model-based-testing-hypotheses/assets/95669100/4514867b-627c-43bf-865b-3cbf3d6fe34a)

##### Under H1 - in one group the probabilities of obtating each score are reversed (as in the figure explaining the data above)
###### 20 observations per group
``` r
simulate_wilcox_olr(samples = 100, n_group = 20, set = 0:5, 
                    arm_1_prob =  rev(c(20, 10, 5, 2, 2, 2)),
                    arm_2_prob =  c(20, 10, 5, 2, 2, 2),
                    which_p = "Wald") %>%  
  plot_differences_between_methods(log_axes = TRUE)
```
![obraz](https://github.com/adrianolszewski/model-based-testing-hypotheses/assets/95669100/b8664220-96d1-4348-a73d-c4172f75a29c)

For Wald's we can notice the already described issue, where some noticeable discrepancies occur. This happens at this sample size, nothing unusual.
In general - at this data size - it's not that bad!

But can read in many textbooks, that at lower sample size Wald's may be unreliable (it does), so the LRT or Rao score is preferred.
Well, sometimes you just NEED Wald's, but if you can ask for the LRT or Rao, do it. Let's have a look:
![obraz](https://github.com/adrianolszewski/model-based-testing-hypotheses/assets/95669100/e6b4e3e2-d4f8-4dc9-b0ca-99c941b1a824)
Indeed, it did bettter! Now it deviates the other way, but definitely in the safe direction - towards lower p-values, farther from the common significant levels.

And Rao's score did even better than LRT, resulting in smaller differences (and ratios):
![obraz](https://github.com/adrianolszewski/model-based-testing-hypotheses/assets/95669100/cbac1d16-54dc-4cb7-8a8f-ee60769fbac1)

###### 30 observations per group
Just 10 more observations makes a big difference and this time Wald's does the job!
I'll skip the other approaches sa they give result consistent with the above. We will compare them later.
![obraz](https://github.com/adrianolszewski/model-based-testing-hypotheses/assets/95669100/bcca869a-330a-4cea-8632-1210575524a4)

###### 50 observations per group
![obraz](https://github.com/adrianolszewski/model-based-testing-hypotheses/assets/95669100/564d2d47-98bf-4b3b-a969-050b710bb526)

##### Under H0 vs H1 - comparison of three estimation methods
Let's make some feeling on how the 3 methods behave.
###### 20 observations per group
![obraz](https://github.com/adrianolszewski/model-based-testing-hypotheses/assets/95669100/c4a96f36-2b38-4cc0-bc03-2fc9ef1e4da0)

###### 30 observations per group
![obraz](https://github.com/adrianolszewski/model-based-testing-hypotheses/assets/95669100/98ed9c48-dca4-40c6-b10a-b49fc4fca5d2)

###### 50 observations per group
![obraz](https://github.com/adrianolszewski/model-based-testing-hypotheses/assets/95669100/a1de3100-d82c-44d4-8f85-75e909c2b42d)

###### 100 observations per group
![obraz](https://github.com/adrianolszewski/model-based-testing-hypotheses/assets/95669100/aac5871e-f7b2-460e-82fa-96a941a4e4b1)

##### Again under H1 - another, less extreme setting
###### 20 observations per group
``` r
simulate_wilcox_olr(samples = 100, n_group = 20, set = 0:5, 
                    arm_1_prob =  rev(c(20, 10, 20, 5, 5, 5)),
                    arm_2_prob =  c(20, 10, 20, 5, 5, 5),
                    which_p = "Wald") %>%  
plot_differences_between_methods(log_axes = TRUE)
```
If the differences between compared groups are not that extreme, the log-likelihood curve becomes not only concave but approaches quadratic shape, so the
differences between the 3 methods and the test (kind of Rao test) become milder.
![obraz](https://github.com/adrianolszewski/model-based-testing-hypotheses/assets/95669100/690fd701-c717-4e97-a062-8f51c356b9dc)

###### 30 observations per group
As the sample size increases, the pattern remains - the discrepancy occurs, but - naturally - at smaller and smaller p-values, where it's completely not harmful. So let's skip checking at bigger samples to save time.
![obraz](https://github.com/adrianolszewski/model-based-testing-hypotheses/assets/95669100/20ce5c77-6b4a-463f-b3d0-e3025207ecc8)

##### Under H0 vs H1 - comparison of three estimation methods
Again, let's make some feeling on how the 3 methods behave. As the differences are less extreme, we expect the methods to give a bit closer results.
But don't expect a spectacular improvement.
###### 20 observations per group
![obraz](https://github.com/adrianolszewski/model-based-testing-hypotheses/assets/95669100/cae95c76-ec7d-4dea-9935-bb21410e376f)

###### 30 observations per group
![obraz](https://github.com/adrianolszewski/model-based-testing-hypotheses/assets/95669100/4e9481e0-dfcd-4076-b0a6-6060d8ae2441)

###### 50 observations per group
![obraz](https://github.com/adrianolszewski/model-based-testing-hypotheses/assets/95669100/36429657-32e5-4fd3-b76d-b5f92f9c2b2e)

###### 100 observations per group
![obraz](https://github.com/adrianolszewski/model-based-testing-hypotheses/assets/95669100/640d7a24-a14f-4bd3-a752-3314801a4ddf)

Well, I must say that despite the inevitable divergence, the improvement actually happened. It's always nice to see the theory agreeing with observation :)

----

#### Comparing numerical data
So far we were comparing Wilcoxon test and ordinal logistic regression model applied to categorical ordinal data.
The agreement was very close, so the selection depended exclusively on the concrete scenario.

The situation complicates when we want to extend comparisons to numeric continuous data.
The ordinal logistic regression can easily handle numerical data, treating them like classes, and esitmating as many intercepts as many unique values exist.

But it won't go so smoothly. For some distributions the procedure may fail to converge, regardless of the implementation, under the alternative hypothesis. Shifting mean in one gruop just by a little (e.g. 0.01) may solve the problem while not affecting the interpretation of the test.

Why is that we explained in file "readme.md". Briefly: the proportional-odds model assesses so-called stochastic superiority, which is probability that for a random pair of observations sampled from 2 groups {a1 ∈ A1, a2 ∈ A2}, a1 > a2. Now, if you have two numeric variables with non-overlapping empirical distributions, this will always hold, so the model won't be able to converge.

Notice, this will almost never happen with Likert items - both variables share same set of values (e.g. 0..5), so they will overlap anyway, except just a few marginal cases, where all observations in one group are lower than in t he other group, e.g. A = {0,0,0,1,1,1,2,2}, B = {3,3,3,4,4,4,5,5}

**Warning: Better choose the implementation that raises errors rather than returns unreliable estimates(!)**

For other distributions it works very well. So let's explore a few examples.

We only need to change a bit the function used for simulations:
```r
simulate_wilcox_olr_distr <- function(samples, n_group, arm_1_distr, arm_2_distr, title) {
  set.seed(1000)
  
  data.frame( 
    do.call( 
      rbind, 
      lapply(1:samples, function(i) {
        print(i)
        stack(
          data.frame(arm1 = eval(parse(text = gsub("n_group", n_group, arm_1_distr))),
                     arm2 = eval(parse(text = gsub("n_group", n_group, arm_2_distr))))) %>% 
          mutate(values_ord = ordered(values), ind = factor(ind)) -> data
        
        c(iter = i,
          n_group = n_group,
          Test  = tidy(wilcox.test(values~ind, data=data, exact = FALSE, adjust = FALSE))$p.value,
          Model = tryCatch(joint_tests(rms::orm(values_ord ~ ind, data))$p.value,
                           error = function(e) {cli::cli_alert_danger("Failed to converge; skip"); NA}))
      }))) -> result
  
  attr(result, "properties") <- list(under_null = identical(arm_1_distr, arm_2_distr),
                                     title = title,
                                     samples = samples,
                                     n_group = n_group)
  return(result)
}
```
##### Under H0:
###### Standard normal distribution. N(0, 1) vs. N(0, 1); n_group = 50 (to save time; for bigger value it's even better)
![obraz](https://github.com/adrianolszewski/model-based-testing-hypotheses/assets/95669100/5481b23d-b715-4e16-b26b-06710dbe49c9)

###### Beta right-skewed β(1, 10) vs. β(1, 10); n_group = 50
![obraz](https://github.com/adrianolszewski/model-based-testing-hypotheses/assets/95669100/01301273-2567-4b20-af4b-25639a4d7bbf)

###### Beta left-skewed β(10, 1) vs. β(10, 1); n_group = 50
![obraz](https://github.com/adrianolszewski/model-based-testing-hypotheses/assets/95669100/94bd9a3e-cc6c-4f36-81c1-ef24605e6415)

###### Beta U-shaped β(0.5, 0.5) vs. β(0.5, 0.5); n_group = 100
![obraz](https://github.com/adrianolszewski/model-based-testing-hypotheses/assets/95669100/b64ebbf8-4e92-440b-9a49-7f7afd238552)

###### Gamma right-skewed Γ(3, 1) vs. Γ(3, 1); n_group = 50
![obraz](https://github.com/adrianolszewski/model-based-testing-hypotheses/assets/95669100/a401d6c7-921b-4c4e-9dcb-f1fe29548438)

**Remarks:**
1. OK, as long as we test similar distributions (under the null), the agreement is perfect regardless of the distribution (Mann-Whitney is a non-parametric test).
2. We can observe, that bias -test giving a little bigger p-values- is common to all cases, but it does not harm at this magnitude. It's just observable.

##### Under H1 (various cases)
###### Standard normal distribution. N(5, 1) vs. N(0, 1); n_group = 20 
```r
simulate_wilcox_olr_distr(samples = 100, n_group = 20, 
                          arm_1_distr = "rnorm(n_group, 5, 1)",
                          arm_2_distr = "rnorm(n_group, 0, 1)", 
                          title = "N(5,1) vs N(0,1)",
                          which_p = "Wald") %>% 
  filter(complete.cases(.)) %>% 
  plot_differences_between_methods(log_axes = TRUE)
```
![obraz](https://github.com/adrianolszewski/model-based-testing-hypotheses/assets/95669100/d8a15792-3e38-4317-812e-2712750fe5c2)

Well, a disaster! It converged(?) for just a few cases and resulted in weird p-values. All the same p-values!

OK, what about Rao?
![obraz](https://github.com/adrianolszewski/model-based-testing-hypotheses/assets/95669100/0eba5c40-246d-4da7-94c6-e006f67fed79)

Well, let's try bigger sample. But I suspect it will just throw an error...

###### Standard normal distribution. N(5, 1) vs. N(0, 1); n_group = 100 
```r
> simulate_wilcox_olr_distr(samples = 100, n_group = 100, 
+                           arm_1_distr = "rnorm(n_group, 5, 1)",
+                           arm_2_distr = "rnorm(n_group, 0, 1)", 
+                           title = "N(5,1) vs N(0,1)",
+                           which_p = "Wald") %>% 
+   filter(complete.cases(.)) %>% 
+   plot_differences_between_methods(log_axes = TRUE)
[1] 1
Did not converge in 12 iterations
Unable to fit model using  “orm.fit” 
[1] 2
Did not converge in 12 iterations
Unable to fit model using  “orm.fit” 
[1] 3
Did not converge in 12 iterations
Unable to fit model using  “orm.fit” 
```
As I expected. The model failed itself.

There may be other implementation that are "silent", but won't make a miracle too:

```r
> set.seed(1000)
> 
> stack(data.frame(arm1 = rnorm(n = 20, mean=5, sd = 1),
+                  arm2 = rnorm(n = 20, mean=0, sd = 1))) %>% 
+           mutate(values_ord = ordered(values), ind = factor(ind)) %>% 
+     MASS::polr(values_ord ~ ind , data = ., Hess = TRUE)
Call:
MASS::polr(formula = values_ord ~ ind, data = ., Hess = TRUE)

Coefficients:
  indarm2 
-19.78993 

Intercepts:
    -1.78384413571217|-1.766198610056     -1.766198610056|-1.34835905258526   -1.34835905258526|-1.22701600631211 
                        -22.732515115                         -21.986947341                         -21.524205629 
[... CUT ...]

Residual Deviance: 239.6635 
AIC: 319.6635 
Warning: did not converge as iteration limit reached
```

It will just not work for this data, sorry.
