# Proving the equivalence between the 2-sample Wald’s z-statistic for comparing proportions with pooled variances and the Pearson’s χ2 (independence) test for a 2×2 contingency table.
#### Adrian Olszewski, 27.03.2025

The Wald’s z-statistic for difference in 2 proportions with unpooled variances is of the following form:

``` math
\begin{equation}
z=\ \frac{\widehat{p_B}-\widehat{p_A}}{\sqrt{p\left(1-p\right)\left(\frac{1}{n_A}+\frac{1}{n_B}\right)}};\ z~N(0,\ 1)
\hspace{2cm} (1)
\end{equation}
```
Where $\widehat{p_A}\$ stands for the estimated probability (sample proportion, %) in the 1st group, 
$\widehat{p_B}$ is the estimated probability in the 2nd group, $n_A$ and $n_B$ denote respective group sizes, and p is the pooled probability $p=\frac{x_A+x_B}{n_A+n_B}=\frac{p_An_A+p_Bn_B}{n_A+n_B}$.

Traditionally the Wald’s statistic is expressed in the squared form, becoming: $z^2=\chi_{df=1}^2$.
Both forms yield the same p-value. For convenience I will show that this  $\chi_{df=1}^2$ statistic is 1:1 equivalent to the $χ^2$ (independence) test for a 2×2 contingency table.

Let’s define the form of the 2×2 contingency table with the observed (O) counts:
| Status   <br>                      Group |  Success <br>(S) | Failure<br> (F) | Total per group<br>(ni·) |
|-----------------------------------------:|:----------------:|:---------------:|:------------------------:|
|                                    **A**    |     $n_{11} = O_{SA}$    |    $n_{12} = O_{FA}$    | $n_A = n_{1·}= n_{11} + n_{12}$       |
|                                     **B**   |     $n_{21} = O_{SB}$    |    $n_{22} = O_{FB}$    | $n_B = n_{2·}= n_{21} + n_{22}$       |
|                   **Total per status (n·i)** | $n_{·1} = n_{11} + n_{21}$  | $n_{·2} = n_{12} + n_{22}$ | $n (=n_A + n_B)$             |

Where $O_{SA}$ stands for “observed number of successes in group A”, $O_{FA}$ stands for “observed number of failures in group A”, and so on.
The test statistic is defined as follows:

``` math
\begin{equation}
X^2=\sum_{r=1}^{2}\sum_{c=1}^{2}\frac{\left(O_{rc}-E_{rc}\right)^2}{E_{rc}};\ X^2~\chi^2(df=1)
\hspace{2cm} (2)
\end{equation}
```

Let’s express the observed and expected number of events as products of totals and probabilities, remembering that the expected number (E) is calculated under $H_0$, i.e. using pooled probability p:

``` math
\begin{equation}
O = \begin{cases} 
O_{11}=O_{SA}=n_Ap_A \\ 
O_{12}=O_{FA}=n_A\left(1-p_A\right) \\
O_{21}=O_{SB}=n_Bp_B \\
{O_{22}=O}_{FB}=n_B\left(1-p_B\right) 
\end{cases}
\hspace{2cm} (3)
\end{equation}
```
and
``` math
\begin{equation}
E = \begin{cases} 
E_{11}=E_{SA}=n_Ap \\ 
E_{12}=E_{FA}=n_A\left(1-p\right) \\
E_{21}=E_{SB}=n_Bp \\
{E_{22}=E}_{FB}=n_B\left(1-p\right)
\end{cases}
\hspace{2cm} (4)
\end{equation}
```

Let’s substitute the O and E elements in the $X^2$ test statistic:

``` math
\begin{equation}
\begin{aligned}
X^2&=\frac{\left(n_Ap_A-n_Ap\right)^2}{n_Ap}+\frac{\left(n_A\left(1-p_A\right)-n_A\left(1-p\right)\right)^2}{n_A\left(1-p\right)}+\frac{\left(n_Bp_B-n_Bp\right)^2}{n_Bp}+\frac{\left(n_B\left(1-p_B\right)-n_B\left(1-p\right)\right)^2}{n_B\left(1-p\right)} \\
&=\frac{\left(n_A\left(p_A-p\right)\right)^2}{n_Ap}+\frac{\left(n_A\left(p-p_A\right)\right)^2}{n_A\left(1-p\right)}+\frac{\left(n_B\left(p_B-p\right)\right)^2}{n_Bp}+\frac{\left(n_B\left(p-p_B\right)\right)^2}{n_B\left(1-p\right)} \\
&=\frac{{n_A\left(p_A-p\right)}^2}{p}+\frac{{n_A\left(p-p_A\right)}^2}{\left(1-p\right)}+\frac{{n_B\left(p_B-p\right)}^2}{p}+\frac{{n_B\left(p-p_B\right)}^2}{\left(1-p\right)} \\
&=\frac{\left(1-p\right)n_A\left(p_A-p\right)^2+pn_A\left(p-p_A\right)^2}{p\left(1-p\right)}+\frac{{{\left(1-p\right)n}_B\left(p_B-p\right)}^2+pn_B\left(p-p_B\right)^2}{p\left(1-p\right)} \\
&=\frac{n_A\left(p_A-p\right)^2\left(1-p+p\right)}{p\left(1-p\right)}+\frac{n_B\left(p_B-p\right)^2\left(1-p+p\right)}{p\left(1-p\right)}=\frac{n_A\left(p_A-p\right)^2}{p\left(1-p\right)}+\frac{n_B\left(p_B-p\right)^2}{p\left(1-p\right)} \\
&=\frac{n_A\left(p_A-\frac{p_An_A+p_Bn_B}{n_A+n_B}\right)^2+n_B\left(p_B-\frac{p_An_A+p_Bn_B}{n_A+n_B}\right)^2}{p\left(1-p\right)} \\
&=\frac{n_A\left(\frac{p_An_B-p_Bn_B}{n_A+n_B}\right)^2+n_B\left(\frac{p_Bn_A-p_An_A}{n_A+n_B}\right)^2}{p\left(1-p\right)}=\frac{n_A\left(\frac{n_B\left(p_A-p_B\right)}{n_A+n_B}\right)^2+n_B\left(\frac{n_A\left(p_B-p_A\right)}{n_A+n_B}\right)^2}{p\left(1-p\right)} \\
&=\frac{\frac{n_An_B^2\left(p_A-p_B\right)^2}{\left(n_A+n_B\right)^2}+\frac{n_Bn_A^2\left(p_A-p_B\right)^2}{\left(n_A+n_B\right)^2}}{p\left(1-p\right)}=\frac{\frac{\left(p_A-p_B\right)^2\left(n_An_B^2+n_Bn_A^2\right)}{\left(n_A+n_B\right)^2}}{p\left(1-p\right)}=\frac{\frac{{{(p}_A-p_B)}^2\left(n_An_B\right)\left(n_B+n_A\right)}{\left(n_A+n_B\right)^2}}{p\left(1-p\right)} \\
\end{aligned}
\hspace{2cm} (5a)
\end{equation}
```

``` math
\begin{equation}
\begin{aligned}
&=\frac{\frac{{{n_An_B(p}_A-p_B)}^2}{n_A+n_B}}{p\left(1-p\right)}=\frac{{{(p}_A-p_B)}^2}{p\left(1-p\right)\frac{n_A+n_B}{n_An_B}}=\frac{{{(p}_A-p_B)}^2}{p\left(1-p\right)\left(\frac{1}{n_A}+\frac{1}{n_B}\right)}\blacksquare
\end{aligned}
\hspace{2cm} (5b)
\end{equation}
```
This way I have proven that the z2 statistic is equivalent to the Pearson’s χ2 statistic for 2×2 table.

-----
By the way, it is worth noticing, that in the χ2 test, the expected frequencies (E) in each column are based on the pooled proportion [p] which is just the weighted average proportion across both groups. 
So, in other words, this test evaluates how much each group’s proportion $p_A$ and $p_B$ deviates from this overall average (p). And this effectively comparing the two groups directly, which is what the z-test actually does.

-----

``` r
> (m <- matrix(c(16, 2, 12, 11),
+              nrow = 2, ncol = 2,
+              dimnames=list(c("A", "B"),
+                            c("Success", "Failure"))))
  Success Failure
A      16      12
B       2      11
> 
> prop.test(m, correct = FALSE)

	2-sample test for equality of proportions without continuity correction

data:  m
X-squared = 6.2859, df = 1, p-value = 0.01217
alternative hypothesis: two.sided
95 percent confidence interval:
 0.1491317 0.6860332
sample estimates:
   prop 1    prop 2 
0.5714286 0.1538462 

> chisq.test(m, correct = FALSE)

	Pearson's Chi-squared test

data:  m
X-squared = 6.2859, df = 1, p-value = 0.01217
```

![obraz](https://github.com/user-attachments/assets/96158adb-97c5-4299-98f2-16a132760912)
