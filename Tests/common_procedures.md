Definition of a helper function for printing the figures.
I use: dplyr, tidyr, ggplot2, ggrepel, and patchwork

``` r
plot_differences_between_methods <- function(results, sign_level = 0.05, log_axes=FALSE) {
  
  properties <- attr(results, "properties")
  samples <- properties$samples
  n_group <- properties$n_group
  hypothesis <- case_when(is.na(properties$under_null) ~ "Mixed hypotheses",
                          properties$under_null ~ "Under H0",
                          !properties$under_null ~ "Under H1")
  under_null <- !is.na(properties$under_null) & properties$under_null

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
                          hypothesis))
)

((p_rej_discrep + p_pval_ratios +
  p_pval_diffs + p_pval_rel +
  patchwork::plot_layout(ncol = 2, nrow = 2)) | p_pval_vs) + 
  plot_layout(widths = c(2, 1))
}
```
