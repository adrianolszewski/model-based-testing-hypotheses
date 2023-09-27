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
3. Raw differences in p-values. These differences do depend on the magnitude of p-values, and that's exactly what I wanted to see. At very small p-values (p<0.001) I would not care
about differences like p=0.00032 vs. p=0.00017 because both results are far from the typical signfiicance levels (0.1, 0.01, 0.05, 0.01 and also 0.001).
For practical purposes I would care if they were discrepant near 0.05.
4. Comparison of % of cases where Test > Model vs. Test < Model, with additional descriptie statistics of the p-values. Because raw differences depend on the magnitude, it's difficult to observe them on a common scale. Log transformation doesn't help there (0, negative values).
So this approach gives me some feelings about the situation
5. How many times both methods disagreed on the rejection of the null hypothesis and what were the actual p-values.

I repeated the simulations sample sizes:
- n_group = 20 - typical small data size; smaller samples don't make much sense in a reliable research,
- n_group = 30 - the "magical number" so many statisticians find "magical". Of course it's not, but this constitutes a reasonable "lower limit of data"
- n_group = 50 - still small data but in a safer area
- n_group = 100 - safe area

Actually, the type-1 error and power is totally off topic in this simulation, we only check how well the methods follow each other, but having the opportunity - we will look at it too. We only have to remember that this is a long-run property and at 100 replications it will be "just close" rather than "nominal", especiallly at smaller samples. It's normal to have, say 7/100 (=0.07), 12/200 (=0.06) to finally reach 15/300 (=0.05).

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
               arm2 = sample(0:5, size=100, replace=TRUE, prob = c(2, 2, 2, 5, 10, 20)))) %>% 
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
simulate_wilcox_olr <- function(samples, n_group, set, arm_1_prob, arm_2_prob) {
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
        
        c(iter = i,
          n_group = n_group,
          Test  = tidy(wilcox.test(values~ind, data=data, exact = FALSE, adjust = FALSE))$p.value,
          Model = joint_tests(MASS::polr(values_ord ~ ind , data = data, Hess=TRUE))$p.value)
      }))) -> result
  
  attr(result, "properties") <- list(arm_1_prob = arm_1_prob, 
                                     arm_2_prob = arm_2_prob, 
                                     under_null = ifelse(mixed_hyp, NA, identical(arm_1_prob, arm_2_prob)),
                                     samples = samples,
                                     n_group = n_group,
                                     title = "")
  return(result)
}
```
#### Results
##### Under H0 - same probabilities of obtating each score in both groups = same ordering of observations
###### 20 observations per group
``` r
simulate_wilcox_olr(samples = 100, n_group = 20, set = 0:5, 
                    arm_1_prob =  c(20, 10, 5, 2, 2, 2),
                    arm_2_prob =  c(20, 10, 5, 2, 2, 2)) %>%  
plot_differences_between_methods()
```
![obraz](https://github.com/adrianolszewski/Logistic-regression-is-regression/assets/95669100/22ab75ac-11b4-40d2-ac82-18f440b52706)

OK, let's make some notes:
1. look at the most right plot. The p-values from both methods are well aligned to the line of slope = 1, though some bias is visible towards Test p-values.
2. It means, that test was more conservative than model.
3. The fact, that practically all p-values coming a test were larger from the p-values obtained from the model is confirmed also by the bar plot (central-bottom).
4. When we look at the area below the 0.05 significance level, we notice 9 obervations. Remembering we "work under the null" it means, that both methods exceeded the 5% significance level, reaching 9% in this simulation. At this sample size it's fine for exploratory analyses.
5. What's nice, the model and test did not contradict each other. 

###### 30 observations per group
![obraz](https://github.com/adrianolszewski/model-based-testing-hypotheses/assets/95669100/4c1a0b6b-3217-42e8-904c-ee9a99d6c17c)

1. The p-values from both methods are very well aligned. A small bias is visible towards test
2. It means, that test was more conservative than model. But the discrepancies are small: between 0.005 and 0.01
3. When we look at the area below the 0.05 significance level, we notice 7 false rejections for the model and 6 for test. Remembering we "work under the null" it means, that both methods exceeded the 5% significance level, but not that much. At this sample size it's a pretty good result.
4. The test and model did not contradict each other except for 1 case, where the difference in p-values very small - 0.002. 

###### 50 observations per group
![obraz](https://github.com/adrianolszewski/model-based-testing-hypotheses/assets/95669100/0816b225-26ab-42d2-872c-7f090fc874d9)

1. Perfect alignment, yet a minimal bias towards test (more conservative) is visible and permanent (100%)
2. From the other side, the differences between test and model are mostly <0.006 (under the null!) - about 2%, so this bias is completely off importance.
3. Both classic and the model-based tests kept the nominal 5% type-1 error

##### 100 observations per group
![obraz](https://github.com/adrianolszewski/model-based-testing-hypotheses/assets/95669100/95c92419-3c47-4379-af16-f0a079ddf41b)

We can end here, as at this sample size no big improvement can be found.
Yes, the bias exists and will continue to exist, but the differences between p-values from both methods mostly don't differ by more than 1% and we can totally ignore it.

##### Under H1 - in one group the probabilities of obtating each score are reversed (as in the figure explaining the data above)
###### 20 observations per group
``` r
simulate_wilcox_olr(samples = 100, n_group = 20, set = 0:5, 
                    arm_1_prob =  rev(c(20, 10, 5, 2, 2, 2)),
                    arm_2_prob =  c(20, 10, 5, 2, 2, 2)) %>%  
plot_differences_between_methods()
```
![obraz](https://github.com/adrianolszewski/model-based-testing-hypotheses/assets/95669100/031ff44c-288f-454c-8774-f5dab2062615)

Well! Let's make some notes:
1. Now we went to small and very small p-values, so I log-transformed both axes.
1. There was a relationship between the methods, but the agreement was not perfect, showing a noticeable bias.
2. This time the model was giving larger p-values, which was confirmed also by the bar plot (central-bottom).
3. There were two cases where the discrepancy was **total**: the model gave p-values close to 1, while the test resulted in p-values much lower than 0.000001. Discrepancy of such a magnitude is rather unusual! Maybe the model failed to convetge reliably.
4. Below 0.00-1 the consistency broke completely, discrepancy reached a 1-2 orders of magnitude
5. **Luckily, for practicaly purposes such discrepancy (at this magnitude of p-values) is totally negligible.**

###### 30 observations per group
![obraz](https://github.com/adrianolszewski/model-based-testing-hypotheses/assets/95669100/71c27faa-0f94-40d6-950d-d45336ae7880)

1. Compared to the 20-sample case, the agreement between both methods is a little bit better
2. As previously, 0.000001 the relationship breaks
3. This time there were no cases where the p-values were opposite (~1 vs. ~0)

###### 50 observations per group
![obraz](https://github.com/adrianolszewski/model-based-testing-hypotheses/assets/95669100/acf2fbcb-e33c-4835-ac24-540d5bf5516d)

1. The divergence pattern is kept and big (2 order of magnitude). Again, it happens at extremely small p-values, where it is completely negligible.

##### Again under H1 - another (less extreme) setting
###### 20 observations per group
``` r
simulate_wilcox_olr(samples = 100, n_group = 20, set = 0:5, 
                    arm_1_prob =  rev(c(20, 10, 20, 5, 5, 5)),
                    arm_2_prob =  c(20, 10, 20, 5, 5, 5)) %>%  
  plot_differences_between_methods(log_axes = TRUE)
```
![obraz](https://github.com/adrianolszewski/model-based-testing-hypotheses/assets/95669100/538fb810-d1fb-4bf8-9981-2535b5c284f7)

1. The consistency is MUCH better.
2. Below 0.0001 it brokes again, but of much smaller magnitude (<1 order of)
3. The closer to H0, the more "peacefully" both methods behave.

###### 30 observations per group
![obraz](https://github.com/adrianolszewski/model-based-testing-hypotheses/assets/95669100/16ef31ac-e370-450e-8065-371e64ea1a11)

The consistency is noticeably better compared to the 20-sample case.

###### 50 observations per group
![obraz](https://github.com/adrianolszewski/model-based-testing-hypotheses/assets/95669100/dcfb1373-3264-484a-9b7a-600d453d5616)

##### Under mixed conditions
###### 20 observations per group
The probabilties for the scores were sampled with varying probabilities
``` r
simulate_wilcox_olr(samples = 100, n_group = 20, set = 0:5, 
                    arm_1_prob =  NULL,
                    arm_2_prob =  NULL) %>% 
  plot_differences_between_methods(log_axes = TRUE)
```
![obraz](https://github.com/adrianolszewski/model-based-testing-hypotheses/assets/95669100/faa99b07-c5b4-4161-b0c3-ca887876dec1)

1. Very good consistency of results!
2. In all cases model was more sensitive than test
3. There were 3 discrepancies but they were caused by minimal differences, so the situation is very good.

###### 30 observations per group
![obraz](https://github.com/adrianolszewski/model-based-testing-hypotheses/assets/95669100/127c5b2c-5d75-451c-bae3-2d79f6ad12f4)

###### 50 observations per group
![obraz](https://github.com/adrianolszewski/model-based-testing-hypotheses/assets/95669100/b708dcff-6535-4a32-948f-06a9ce0f290d)

----

#### Comparing numerical data
So far we were comparing Wilcoxon test and ordinal logistic regression model applied to categorical ordinal data.
The agreement was very close, so the selection depended exclusively on the concrete scenario.

The situation complicates when we want to extend comparisons to numeric continuous data.
The ordinal logistic regression can easily handle numerical data, treating them like classes, and esitmating as many intercepts as many unique values exist.

But it won't go so smoothly. For some distributions the procedure may fail to converge, regardless of the implementation, under the alternative hypothesis. Adding just a little value to certain parameters (e.g. 0.01 to the mean parameter) may solve the problem.

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
###### Standard normal distribution. N(5, 1) vs. N(0, 1); n_group = 50 
```r
simulate_wilcox_olr_distr(samples = 100, n_group = 50, 
                          arm_1_distr = "rnorm(n_group, 5, 1)",
                          arm_2_distr = "rnorm(n_group, 0, 1)",
                          title = "N(5,1) vs N(0,1)") %>% 
  filter(complete.cases(.)) %>% 
  plot_differences_between_methods(log_axes = FALSE)

[1] 1
Did not converge in 12 iterations
Unable to fit model using  “orm.fit” 
✖ Failed to converge; skip
[1] 2
Did not converge in 12 iterations
Unable to fit model using  “orm.fit” 
✖ Failed to converge; skip
[1] 3
Did not converge in 12 iterations
Unable to fit model using  “orm.fit” 
✖ Failed to converge; skip
```

Wow! A disaster! It converged for just a few cases and resulted in weird p-values.
BTW: the MASS::polr() function won't report any issue, but the p-values will be totally wrong.
So, again, it's better to see the problem.

Look below, if you are not convinced (I temporarily replaced rms::orm() with MASS::polr())
![obraz](https://github.com/adrianolszewski/model-based-testing-hypotheses/assets/95669100/c6df164c-23c7-4dbb-9dac-372ed1ac89f0)

Convinced now? OK, let's try another families.

###### Beta right-skewed β(1, 10) vs. β(10, 1); n_group = 50  (for simplicity; It's gives oppositely skewed distributions, like in the example with Likert items)
``` r
> simulate_wilcox_olr_distr(samples = 100, n_group = 50, 
+                           arm_1_distr = "rbeta(n_group, 10, 1)",
+                           arm_2_distr = "rbeta(n_group, 1, 10)",
+                           title = "β(10,1) vs β(1,10)") %>% 
+   filter(complete.cases(.)) %>% 
+   plot_differences_between_methods(log_axes = TRUE)
[1] 1
Did not converge in 12 iterations
Unable to fit model using  “orm.fit” 
✖ Failed to converge; skip
[1] 2
Did not converge in 12 iterations
Unable to fit model using  “orm.fit” 
✖ Failed to converge; skip
[1] 3
Did not converge in 12 iterations
Unable to fit model using  “orm.fit” 
✖ Failed to converge; skip
```

###### Gamma right-skewed Γ(3, 1) vs. Γ(1, 3); n_group = 50
Surprise! Just a few errors. The methods agreed that the null hypothesis should be rejected, but the p-values were a "free-ride". Luckily, MUCH below any classic significance level.
![obraz](https://github.com/adrianolszewski/model-based-testing-hypotheses/assets/95669100/4ce571e9-e3f7-4030-afa7-6bc6730ce5be)

But it seems these distributions were too similar.
When I took two more differing from each other, it failed again:

``` r
simulate_wilcox_olr_distr(samples = 100, n_group = 50, 
                          arm_1_distr = "rgamma(n_group, 1, 10)",
                          arm_2_distr = "rgamma(n_group, 10, 10)",
                          title = "Γ(1,10) vs Γ(10,10)") %>% 
  filter(complete.cases(.)) %>% 
  plot_differences_between_methods(log_axes = TRUE)

[1] 1
[1] 2
Did not converge in 12 iterations
Unable to fit model using  “orm.fit” 
✖ Failed to converge; skip
[1] 3
Did not converge in 12 iterations
Unable to fit model using  “orm.fit” 
✖ Failed to converge; skip
```
Although the still - whenever possible - both methods rejected the H0.
![obraz](https://github.com/adrianolszewski/model-based-testing-hypotheses/assets/95669100/75d16b9d-6b3f-462b-a124-5e8bc3dcd163)

It looks like under H1 the model-based approach fails completely for the numerical variables.

###### Further investigation
But why?
- It's not about too many intercepts (under H0 it worked well)
- It's not about negative values (under H0 it worked well)

So why? I don't know yet.

When I increased variance a lot and rounded to integers, it started behaving much better:
``` r
simulate_wilcox_olr_distr(samples = 100, n_group = 50, 
                          arm_1_distr = "round(rnorm(50, 10, 5), 0)",
                          arm_2_distr = "round(rnorm(50, 20, 5), 0)",
                          title = "[N(10,1)] vs [N(20,1)]") %>% 
  filter(complete.cases(.)) %>% 
  plot_differences_between_methods(log_axes = TRUE)
```
No errors. Quite a good agreement between the test and model results.

![obraz](https://github.com/adrianolszewski/model-based-testing-hypotheses/assets/95669100/f5524268-6749-48be-9981-74c86826a5a9)

So maybe it's the problem with variance?

I returned to exact (not rounded) values, only increased standard deviation and...
```r
simulate_wilcox_olr_distr(samples = 100, n_group = 50, 
                          arm_1_distr = "rnorm(50, 0, 5)",
                          arm_2_distr = "rnorm(50, 5, 5)",
                          title = "N(0,5) vs N(5,5)") %>% 
  filter(complete.cases(.)) %>% 
  plot_differences_between_methods(log_axes = TRUE)
```
... no errors!
... and perfect agreement!
![obraz](https://github.com/adrianolszewski/model-based-testing-hypotheses/assets/95669100/d8501939-0361-4342-9c04-4f1d0d3227d7)

So what's going on?!

Now we also know it's not a matter of comparing fractional values too (if anyone asked).

**It's probably a matter of dispersion (variance)?**

Let's return to gamma distribution, but not make it wider (increase variance):
``` r
simulate_wilcox_olr_distr(samples = 100, n_group = 50, 
                          arm_1_distr = "rgamma(n_group, 2, .4)",
                          arm_2_distr = "rgamma(n_group, 4, .4)",
                          title = "Γ(2,.4) vs Γ(4,.4)") %>% 
  filter(complete.cases(.)) %>% 
  plot_differences_between_methods(log_axes = TRUE)
```
![obraz](https://github.com/adrianolszewski/model-based-testing-hypotheses/assets/95669100/4cc8cf5a-d119-476d-b5a7-cb9aeef0bcce)

**And one more: wider N(4,4) vs. Γ(4,.4)**
![obraz](https://github.com/adrianolszewski/model-based-testing-hypotheses/assets/95669100/f6756416-b662-4c62-9dbb-0f1171331597)

Perfect agreement.

So now let's scale both distributions to a unit-variance:
```r
simulate_wilcox_olr_distr(samples = 100, n_group = 50, 
                          arm_1_distr = "scale(rgamma(n_group, 2, .4), center=FALSE)",
                          arm_2_distr = "scale(rnorm(n_group, 4, 4), center=FALSE)",
                          title = "sc(Γ(2,0.4)) vs sc(Γ(4,0.4))") %>% 
  filter(complete.cases(.)) %>% 
  plot_differences_between_methods(log_axes = TRUE)
```
... no errors and perfect agreement. Now because we made the distributions too similar to each other (similar mean, variance, only shapes differ).

Let's return to the example with beta distribution, now with both distributions scaled to a unit variance:
```r
simulate_wilcox_olr_distr(samples = 100, n_group = 50, 
                          arm_1_distr = "scale(rbeta(n_group, 10, 1), center=FALSE)",
                          arm_2_distr = "scale(rbeta(n_group, 1, 10), center=FALSE)",
                          title = "sc(β(10,1)) vs sc(β(1,10))") %>% 
  filter(complete.cases(.)) %>% 
  plot_differences_between_methods(log_axes = TRUE)
```
![obraz](https://github.com/adrianolszewski/model-based-testing-hypotheses/assets/95669100/624e6ea0-56f8-4453-b067-2a731145d388)

Now it works! OK, there's some bias, but, typically - at the very low p-values.

OK, but what about the normal distribution? When I asked the program to compare two distributions with SAME unit-variance, only differing by location, it failed permanently.
Don't tell me it will work, when I scale it down (to a unit variance again)?

```r
simulate_wilcox_olr_distr(samples = 100, n_group = 50, 
                          arm_1_distr = "scale(rnorm(n_group, 0, 1), center=FALSE)",
                          arm_2_distr = "scale(rnorm(n_group, 5, 1), center=FALSE)",
                          title = "sc(N(0,1)) vs sc(N(5,1))") %>% 
  filter(complete.cases(.)) %>% 
  plot_differences_between_methods(log_axes = TRUE)
```
![obraz](https://github.com/adrianolszewski/model-based-testing-hypotheses/assets/95669100/dbd2db2b-c0a7-476a-a6ad-bf9fcbea9926)

You don't say! So we made a cirle and returned to the case that failed, with similar setting, but now it doesn't fail?

Well, yes I have to investigate it more.

The general conclusion is that:
1. Mann-Whitney (-Wilcoxon) and the proportional-odds model (ordinal logistic regression) perfectly agree for naturally ordinal data, like Likert items and you can use either, depending on purpose. If you are interested in a more complex comparison of several contrasts, adjusting for covariates - then you have the tool for it.
2. For numerical data - well, it really depends on various circumstances and it's difficult to say how your model will behave. In general - if it converges, then it's rather in agreement with Mann-Whitney (-Wilcoxon). If it fails to converge, you have bad luck and that's all. Either you will switch to simple Wilcoxon (and resign from covariate adjustments, comparisons, etc.) or admit it failed, or try another method, maybe quantile regression (with different objective and interpretation!).


   
