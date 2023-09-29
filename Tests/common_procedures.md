# Function sampling Likert items {0-5} with predefined probability of each item occurrence and applying both Mann-Whitney (-Wilcoxon) test and the ordinal logistic regression (proportional-odds model)

I use the orm::olr() as it returns both LRT and Rao score. The Wald's is calculated by hand to save time.
```r
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

----

# Function sampling from specified distributions and applying both Mann-Whitney (-Wilcoxon) test and the ordinal logistic regression (proportional-odds model)
```r

simulate_wilcox_olr_distr <- function(samples, n_group, arm_1_distr, arm_2_distr, title, which_p) {
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
  
  attr(result, "properties") <- list(arm_1_prob = arm_1_distr, 
                                     arm_2_prob = arm_2_distr, 
                                     under_null = identical(arm_1_distr, arm_2_distr),
                                     samples = samples,
                                     n_group = n_group,
                                     title = title,
                                     method = which_p)  
  return(result)
}
```

----
# Comparing MWW test and PO-model for Likert under H0 and H1 (Halt)
```r
compare_likert_h0_h1_3methods <- function(n_group, log_transform=TRUE, 
                                   arm_prob =  c(20, 10, 5, 2, 2, 2)) {
  
  cbind(method="Wald", simulate_wilcox_olr(samples = 100, n_group = n_group, set = 0:5, 
                                           arm_1_prob = arm_prob,
                                           arm_2_prob = arm_prob, which_p = "Wald")) %>% 
    bind_rows(cbind(method="Wilk's LRT", simulate_wilcox_olr(samples = 100, n_group = n_group, set = 0:5, 
                                                             arm_1_prob =  arm_prob,
                                                             arm_2_prob =  arm_prob, which_p = "LRT")),
              cbind(method="Rao's score", simulate_wilcox_olr(samples = 100, n_group = n_group, set = 0:5, 
                                                              arm_1_prob =  arm_prob,
                                                              arm_2_prob =  arm_prob, which_p = "Rao"))) ->cmp_h0


  ph0 <- cmp_h0 %>% 
    ggplot() +
    geom_point(aes(x=Test, y = Model, col=method)) +
    
    geom_abline(slope = 1, intercept = 0, col="green") +
    
    { if(TRUE)
      list(geom_hline(yintercept = 0.05, linetype="dashed", color="red"),
           geom_vline(xintercept = 0.05, linetype="dashed", color="red"))
    } +
    
    { if(FALSE)
      list(scale_y_continuous(trans='log10'),
           scale_x_continuous(trans='log10'))
    } +
    
    xlab("Test p-values") +
    ylab("Model p-values") +
    theme_bw() +
    labs(title = "P-values from Test vs. Model. Comparison of LRT vs. Rao vs. Wald",
         subtitle = sprintf("N=%d samples, group size=%d obs., %s",
                            100, n_group, "under null hypothesis"),
         caption = ifelse(FALSE, "Both axes log10 transformed", "")) +
    theme(plot.caption = element_text(color = "red", face="italic"))


  cbind(method="Wald", simulate_wilcox_olr(samples = 100, n_group = n_group, set = 0:5, 
                                           arm_1_prob =  rev(arm_prob),
                                           arm_2_prob =  arm_prob, which_p = "Wald")) %>% 
    bind_rows(cbind(method="Wilk's LRT", simulate_wilcox_olr(samples = 100, n_group = n_group, set = 0:5, 
                                                             arm_1_prob =  rev(arm_prob),
                                                             arm_2_prob =  arm_prob, which_p = "LRT")),
              cbind(method="Rao's score", simulate_wilcox_olr(samples = 100, n_group = n_group, set = 0:5, 
                                                              arm_1_prob =  rev(arm_prob),
                                                              arm_2_prob = arm_prob, which_p = "Rao"))) ->cmp_h1


  ph1 <- cmp_h1 %>% 
    ggplot() +
    geom_point(aes(x=Test, y = Model, col=method)) +
    
    geom_abline(slope = 1, intercept = 0, col="green") +
    
    { if(FALSE)
      list(geom_hline(yintercept = 0.05, linetype="dashed", color="red"),
           geom_vline(xintercept = 0.05, linetype="dashed", color="red"))
    } +
    
    { if(log_transform)
      list(scale_y_continuous(trans='log10'),
           scale_x_continuous(trans='log10'))
    } +
    
    xlab("Test p-values") +
    ylab("Model p-values") +
    theme_bw() +
    labs(title = "P-values from Test vs. Model. Comparison of LRT vs. Rao vs. Wald",
         subtitle = sprintf("N=%d samples, group size=%d obs., %s",
                            100, n_group, "under alternative hypothesis"),
         caption = ifelse(log_transform, "Both axes log10 transformed", "")) +
    theme(plot.caption = element_text(color = "red", face="italic"))

  (ph0 + ph1 & theme(legend.position = "top")) + plot_layout(guides = "collect")
}
```

----

# Comparing test vs. model for the 3 methods
```r

compare_3methods_distr <- function(n_group, samples=100,
                                   distr_pair_spec = list("N(50,20) vs. N(50,20)" = 
                                                            c(arm_1_distr = "rnorm(n_group, mean=50, sd=20)", 
                                                              arm_2_distr = "rnorm(n_group, mean=50, sd=20)",
                                                              log = "TRUE"))) {
                                   
  lapply(names(distr_pair_spec), function(dps_name) {
    dps <- distr_pair_spec[[dps_name]]

    log_axes <- as.logical(dps["log"])
    log_axes <- ifelse(is.na(log_axes), TRUE, log_axes)
    log_axes <- (dps["arm_1_distr"] != dps["arm_2_distr"]) & log_axes
    
    bind_rows(
      bind_cols(method="Wald", simulate_wilcox_olr_distr(samples = samples, 
                                                         n_group = n_group,
                                                         arm_1_distr = dps["arm_1_distr"],
                                                         arm_2_distr = dps["arm_2_distr"],
                                                         title = dps_name, 
                                                         which_p = "Wald")), 
      bind_cols(method="Wilks' LRT", simulate_wilcox_olr_distr(samples = samples, 
                                                         n_group = n_group,
                                                         arm_1_distr = dps["arm_1_distr"],
                                                         arm_2_distr = dps["arm_2_distr"],
                                                         title = dps_name, 
                                                         which_p = "LRT")),
      bind_cols(method="Rao score", simulate_wilcox_olr_distr(samples = samples, 
                                                               n_group = n_group,
                                                               arm_1_distr = dps["arm_1_distr"],
                                                               arm_2_distr = dps["arm_2_distr"],
                                                               title = dps_name, 
                                                               which_p = "Rao"))) %>%  
      ggplot() +
      
      { if(dps["arm_1_distr"] == dps["arm_2_distr"])
        list(geom_hline(yintercept = 0.05, linetype="dashed", color="red"),
             geom_vline(xintercept = 0.05, linetype="dashed", color="red"))
      } +
      
      { if(log_axes)
        list(scale_y_continuous(trans='log10'),
             scale_x_continuous(trans='log10'))
      } +
      
      geom_point(aes(x=Test, y = Model, col=method)) +
      geom_abline(slope = 1, intercept = 0, col="green") +

      xlab("Test p-values") +
      ylab("Model p-values") +
      theme_bw() +
      labs(title = "P-values from Test vs. Model. Comparison of LRT vs. Rao vs. Wald",
           subtitle = sprintf("N=%d samples, group size=%d obs. | comparison: %s",
                              samples, n_group, dps_name),
           caption = ifelse(log_axes, "Both axes log10 transformed", "")) +
      theme(plot.caption = element_text(color = "red", face="italic"))
  }) -> plots
  
  (patchwork::wrap_plots(plot = plots, ncol = 2) & theme(legend.position = "top")) + plot_layout(guides = "collect")
}
```

----

# Comparing p-values from model only - 3 approaches
```r
compare_3methods_distr1 <- function(n_group, samples=100,
                                   distr_pair_spec = list("N(50,20) vs. N(50,20)" = 
                                                            c(arm_1_distr = "rnorm(n_group, mean=50, sd=20)", 
                                                              arm_2_distr = "rnorm(n_group, mean=50, sd=20)",
                                                              log = "TRUE"))) {
  
  lapply(names(distr_pair_spec), function(dps_name) {
    dps <- distr_pair_spec[[dps_name]]
    
    log_axes <- as.logical(dps["log"])
    log_axes <- ifelse(is.na(log_axes), TRUE, log_axes)
    log_axes <- (dps["arm_1_distr"] != dps["arm_2_distr"]) & log_axes
    
    bind_rows(
      bind_cols(method="Wald", simulate_wilcox_olr_distr(samples = samples, 
                                                         n_group = n_group,
                                                         arm_1_distr = dps["arm_1_distr"],
                                                         arm_2_distr = dps["arm_2_distr"],
                                                         title = dps_name, 
                                                         which_p = "Wald")), 
      bind_cols(method="Wilks' LRT", simulate_wilcox_olr_distr(samples = samples, 
                                                               n_group = n_group,
                                                               arm_1_distr = dps["arm_1_distr"],
                                                               arm_2_distr = dps["arm_2_distr"],
                                                               title = dps_name, 
                                                               which_p = "LRT")),
      bind_cols(method="Rao score", simulate_wilcox_olr_distr(samples = samples, 
                                                              n_group = n_group,
                                                              arm_1_distr = dps["arm_1_distr"],
                                                              arm_2_distr = dps["arm_2_distr"],
                                                              title = dps_name, 
                                                              which_p = "Rao"))) %>%  
      ggplot() +
      
      { if(dps["arm_1_distr"] == dps["arm_2_distr"])
        list(geom_hline(yintercept = 0.05, linetype="dashed", color="red"))
      } +
      
      { if(log_axes)
        list(scale_y_continuous(trans='log10'),
             scale_x_continuous(trans='log10'))
      } +
      
      geom_point(aes(x=iter, y = Model, col=method)) +

      xlab("Test p-values") +
      ylab("Model p-values") +
      theme_bw() +
      labs(title = "P-values: Wilks LRT vs. Rao vs. Wald | Ordinal Logistic Regression",
           subtitle = sprintf("N=%d samples, group size=%d obs. | comparison: %s",
                              samples, n_group, dps_name),
           caption = ifelse(log_axes, "Both axes log10 transformed", "")) +
      theme(plot.caption = element_text(color = "red", face="italic"))
  }) -> plots
  
  (patchwork::wrap_plots(plot = plots, ncol = 2) & theme(legend.position = "top")) + plot_layout(guides = "collect")
}
```

----
# Printing Test vs. Model comparison for 
``` r
plot_differences_between_methods <- function(results, sign_level = 0.05, log_axes=FALSE) {
  
  properties <- attr(results, "properties")
  samples <- properties$samples
  n_group <- properties$n_group
  under_null <- properties$under_null
  title <- properties$title
  hypothesis <- case_when(is.na(under_null) ~ "Mixed hypotheses",
                          under_null ~ "Under H0",
                          .default = "Under H1")
  method <- properties$method
  
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
  ylab(NULL)
)

(p_pval_ratios <- results %>% 
   ggplot() +
   geom_point(aes(x=iter, y = ratio)) +
   geom_hline(yintercept = 1, col="green") +
   theme_bw() +
   scale_y_continuous(trans='log10') +
   labs(title = "Ratio of p-values: Test vs Model",
        subtitle =  sprintf("N=%d samples, group size=%d observations", samples, n_group),
        caption = "Y axis log10 transformed") +
   ylab("p-value Test / p-value Model") +
   theme(plot.caption = element_text(color = "red", face="italic"))
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

  { if(under_null)
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
           scale_x_continuous(trans='log10'))
  } +

  xlab("Test p-values") +
  ylab("Model p-values") +
  theme_bw() +
  labs(title = sprintf("P-values [%s]: Test vs Model", method),
       subtitle = sprintf("N=%d samples, group size=%d obs., %s %s",
                          samples, n_group, 
                          hypothesis, title),
       caption = ifelse(log_axes, "Both axes log10 transformed", "")) +
    theme(plot.caption = element_text(color = "red", face="italic"))
)

(p_pval_vs | (p_rej_discrep + p_pval_ratios +
  p_pval_diffs + p_pval_rel +
  patchwork::plot_layout(ncol = 2, nrow = 2))) + 
  plot_layout(widths = c(1, 2))
}
```
