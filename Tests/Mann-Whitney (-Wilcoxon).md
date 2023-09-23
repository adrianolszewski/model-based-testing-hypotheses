# Mann-Whitney (-Wilcoxon)

## Table of contents
1. Description of the test
2. Model-based implementation via Ordinal Logistic Regression (aka Proportional-Odds model)
3. References
4. [Simulations - how well the model-based approach replicates the classic test](#simulations)

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
- n_group = 100 - typical situation in clinical trials, where I work
- n_group = 200 - above that, even small discrepancies between groups will be reported as statisticall significant, so testing at this data size is dangerous, especially without the practical significance defined!).

Actually, the type-1 error and power is totally off topic in this simulation, we only check how well the methods follow each other.

---
Definition of a helper function for printing the figures.
I use: dplyr, tidyr, ggplot2, ggrepel, and patchwork

``` r
plot_differences_between_methods <- function(results, sign_level = 0.05, log_axes=FALSE) {
  
  properties <- attr(results, "properties")
  samples <- properties$samples
  n_group <- properties$n_group

results %>% 
  mutate(diff = Test - Model,
         diff_sign = case_when(diff == 0 ~ "Test = Model", diff > 0 ~ "Test > Model", diff < 0 ~ "Test < Model"),
         ratio = Test / Model,
         discr = (Test <= sign_level & Model > sign_level) | (Test > sign_level & Model <= sign_level),
         discr_descr = ifelse(discr, sprintf("T: %.3f ; M :%.3f", Test, Model), NA)) -> results 

(p_rej_discrep <- results %>%
  ggplot(aes(x=iter, y = discr, label=ifelse(discr, discr_descr, NA), col=discr)) +
  geom_point(show.legend = FALSE) +
  ggrepel::geom_label_repel(min.segment.length = 0, seed = 100,
                           box.padding = .4,  segment.linetype = 6,   
                           direction = "both",
                           nudge_y = .1,
                           max.overlaps = 50, size=2.5, show.legend = FALSE) +
  scale_color_manual(values=c("TRUE"="red2", "FALSE"="green2")) +
  scale_y_discrete(labels = c('In agreement','Discrepant')) +
  theme_bw() +
  labs(title = "Discrepancy in rejection: Test vs Model",
      subtitle =  sprintf("N=%d samples, group size=%d observations; sign. level=%.3f", samples, n_group, sign_level)) +
  ylab("p-value Test - p-value Model")
)

(p_pval_ratios <- results %>% 
   ggplot() +
   geom_point(aes(x=iter, y = ratio)) +
   geom_hline(yintercept = 1, col="green") +
   theme_bw() +
   labs(title = "Ratio of p-values: Test vs Model",
        subtitle =  sprintf("N=%d samples, group size=%d observations", samples, n_group)) +
   ylab("p-value Test / p-value Model")
) 

(p_pval_diffs <- results %>% 
  ggplot() +
  geom_point(aes(x=iter, y = diff)) +
  geom_hline(yintercept = 0, col="green") +
  theme_bw() +
  labs(title = "Raw differences in p-values: Test vs Model",
      subtitle =  sprintf("N=%d samples, group size=%d observations", samples, n_group)) +
  ylab("p-value Test - p-value Model")
) 

(p_pval_rel <- results %>% 
   group_by(diff_sign) %>% 
   summarize(n=n(),
             q=list(setNames(quantile(diff), nm=c("Min", "Q1", "Med", "Q3", "Max"))),
             .groups = "drop_last") %>% 
   mutate(p=n/sum(n)) %>% 
   unnest_wider(q) %>% 
   
   ggplot() +
   geom_bar(aes(x = diff_sign, y=p), stat="identity") +
   scale_y_continuous(labels=scales::percent) +
   theme_bw() +
   xlab("Relationship") +
   labs(title = "Relationship between p-values: Test vs Model",
        subtitle =  sprintf("N=%d samples, group size=%d observations", samples, n_group)) +
   geom_label(aes(x = diff_sign, y=.5,
                  label=sprintf("Quartiles of differences\nMin=%.3f\nQ1=%.3f\nMed=%.3f\nQ3=%.3f\nMax=%.3f", Min, Q1, Med, Q3, Max)),
              size=3.5) +
   ylab(NULL)
)

(p_pval_vs <- results %>% 
  ggplot() +
  geom_point(aes(x=Test, y = Model)) +
    
  geom_abline(slope = 1, intercept = 0, col="green") +

  { if(properties$under_null)
       list(geom_hline(yintercept = sign_level, linetype="dashed", color="red"),
            geom_vline(xintercept = sign_level, linetype="dashed", color="red"),
            geom_label(aes(x=0.75, y=0.25, 
                           label = sprintf("Rejections of H0 at α=%.3f\nModel: %d\nTest: %d",
                                           sign_level,
                                           sum(Model <= sign_level),
                                           sum(Test <= sign_level)))))
  } +
  
  { if(log_axes)
      list(scale_y_continuous(trans='log10'),
           scale_x_continuous(trans='log10'),
           xlab("Test log(p-values)"),
           ylab("Model log(p-values)"))
    else
      list(xlab("Test p-values"),
           ylab("Model p-values"))
  } +
    
  theme_bw() +
  labs(title = "P-values: Test vs Model",
       subtitle = sprintf("N=%d samples, group size=%d obs., %s",
                          samples, n_group, 
                          ifelse(properties$under_null, "Under H0", "Under H1")))
)

((p_rej_discrep + p_pval_ratios +
  p_pval_diffs + p_pval_rel +
  patchwork::plot_layout(ncol = 2, nrow = 2)) | p_pval_vs) + 
  plot_layout(widths = c(2, 1))
}
```
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
  data.frame( 
    do.call( 
      rbind, 
      lapply(1:samples, function(i) {
        print(i)
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
                                     under_null = identical(arm_1_prob, arm_2_prob),
                                     samples = samples,
                                     n_group = n_group)
  return(result)
}
```
#### Results
##### 20 observations per group
###### Under H0 - same probabilities of obtating each score in both groups = same ordering of observations
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
4. When we look at the area below the 0.05 significance level, we notice 9 obervations. Remembering we "work under the null" it means, that both methods exceeded the 5% significance level, reaching 9% in this simulation. At this sample size it's not bad for exploratory purposes.
5. What's nice, the tests did not contradict each other. 

###### Under H1 - in one group the probabilities of obtating each score are reversed (as in the figure explaining the data above)
``` r
simulate_wilcox_olr(samples = 100, n_group = 20, set = 0:5, 
                    arm_1_prob =  rev(c(20, 10, 5, 2, 2, 2)),
                    arm_2_prob =  c(20, 10, 5, 2, 2, 2)) %>%  
plot_differences_between_methods()
```

Well well weel! Let's make some notes:
1. Now we went to small and very small p-values, so I log-transformed both axes.
1. There was a clear relationship between the methods, but it's far from a perfect agreement!
2. This time model was giving noticeably larger p-values
3. The fact, that practically all p-values coming a test were larger from the p-values obtained from the model is confirmed also by the bar plot (central-bottom).
4. When we look at the area below the 0.05 significance level, we notice 9 obervations. Remembering we "work under the null" it means, that both methods exceeded the 5% significance level, reaching 9% in this simulation. At this sample size it's not bad for exploratory purposes.
5. What's nice, the tests did not contradict each other. 

###### Again under H1 - this time less "aggressively"
OK, this was extreme. What about a more comparable groups? The closer to H0, the more "peacefully" the methods behave :-)
As if below 0.0001 something bad happened. Well, we will observe this pattern in the next simulations.

``` r
simulate_wilcox_olr(samples = 100, n_group = 20, set = 0:5, 
                    arm_1_prob =  rev(c(20, 10, 20, 5, 5, 5)),
                    arm_2_prob =  c(20, 10, 20, 5, 5, 5)) %>%  
  plot_differences_between_methods(log_axes = TRUE)
```

