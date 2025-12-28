# Proving the equivalence between the 2-sample Wald’s z-statistic for comparing proportions with pooled variances and the Rao score test over logistic regression with a single binary predictor
#### Adrian Olszewski, 21.03.2025

The Wald’s z-statistic for difference in 2 proportions with unpooled variances is of the following form:

``` math
\begin{equation}
z=\ \frac{\widehat{p_B}-\widehat{p_A}}{\sqrt{p\left(1-p\right)\left(\frac{1}{n_A}+\frac{1}{n_B}\right)}};\ z~N(0,\ 1)
\hspace{2cm} (1)
\end{equation}
```
Where $\widehat{p_A}\$ stands for the estimated probability (sample proportion, %) in the 1st group, 
$\widehat{p_B}$ is the estimated probability in the 2nd group, $n_A$ and $n_B$ denote respective group sizes, and p is the pooled probability $p=\frac{x_A+x_B}{n_A+n_B}$.

Traditionally the Wald’s statistic is expressed in the squared form, becoming: $z^2=\chi_{df=1}^2$.
Both forms yield the same p-value. For convenience I will show that this  $\chi_{df=1}^2$ statistic is 1:1 equivalent to the Rao score test over the logistic regression with a single binary predictor playing role of indicator for the compared samples.
To simplify calculations, I will derive the formula for pooled probability p and the overall statistic form separately.

**The equation of logistic regression**

Let’s start from the equation of the logistic regression with a single binary predictor:

``` math
\begin{equation}
logit\left(E\left(Y\middle| X\right)=logit\left(\hat{p}(Y=1|X\right)\right)=\ln{\left(\frac{\hat{p}\left(Y=1\middle| X\right)}{1-\hat{p}\left(Y=1\middle| X\right)}\right)}=\beta_0+X_1\beta_1
\hspace{2cm} (2)
\end{equation}
```

or equivalently, after applying the inverse-logit, i.e. sigmoid function (let’s also simplify X<sub>1</sub> to just X)

``` math
\begin{equation}
\hat{p}\left(Y=1\middle| X\right)=\frac{e^{\beta_0+X\beta_1}}{1+e^{\beta_0+X\beta_1}}=\frac{e^\eta}{1+e^\eta}
\hspace{2cm} (3)
\end{equation}
```
where Y<sub>i</sub> are independent Bernoulli random variables with probabilities p<sub>i</sub>. The second, simpler form with η=β<sub>0</sub>+Xβ<sub>1</sub> will facilitate later calculations.

Let’s also encode the two levels {A, B} using a single binary predictor X such that: A: X=0, B: X=1, let’s simplify notation and express $p_A$ and $p_B$ in terms of beta coefficients:

``` math
\begin{equation}
\begin{cases} 
\hat{p}_A = P(Y = 1, X = 0) = \frac{e^{\beta_0}}{1 + e^{\beta_0}} \\ 
\hat{p}_B = P(Y = 1, X = 1) = \frac{e^{\beta_0 + \beta_1}}{1 + e^{\beta_0 + \beta_1}} 
\end{cases}
\hspace{2cm} (4)
\end{equation}
```

Let’s also skip the hat notation for estimated p and use simply $p_A$ and $p_B$ until stated otherwise.

**Pooled probability**

Under $H_0$, i.e. $β_1=0$, $p_A = p_B$, and the best estimator of this common proportion is the pooled proportion p. Let’s find its form.
We need to assume that the data consists of two independent binomial samples:
``` math
\begin{equation}
\begin{cases} 
{group\ A=X}_A~Binomial\left(n_A,p_A\right) \\ 
{group\ B=X}_B~Binomial\left(n_B,p_B\right) 
\end{cases}
\hspace{2cm} (5)
\end{equation}
```

So the likelihood function is:

``` math
\begin{equation}
L(p_A,p_B)=\left(\begin{matrix}n_A\\x_A\\\end{matrix}\right)p_A^{x_A}\left(1-p_A\right)^{n_A-x_A}\bullet\left(\begin{matrix}n_B\\x_B\\\end{matrix}\right)p_B^{x_B}\left(1-p_B\right)^{n_B-x_B}
\hspace{2cm} (6)
\end{equation}
```
where $x_A$ and $x_B$ are the observed number of successes in groups A and B, respectively, and $n_A$ and $n_B$ are the total sample sizes in each group.

Under the null hypothesis $H0: p_A=p_B=p_{pooled}$. The pooled p is of some form yet unknown.
Knowing that:

``` math
\begin{equation}
\begin{aligned}
L(p_{pooled}) &= \binom{n_A}{x_A} \binom{n_B}{x_B} p_{pooled}^{x_A + x_B} (1 - p_{pooled})^{(n_A - x_A) + (n_B - x_B)} \\
&= \binom{n_A}{x_A} \binom{n_B}{x_B} p_{pooled}^{x_A + x_B} (1 - p_{pooled})^{(n_A + n_B - x_A - x_B)}
\end{aligned}
\hspace{2cm} (7)
\end{equation}
```
we obtain a single binomial likelihood where the total number of successes is $x_A+x_B$ and the total number of trials is $n_A+n_B$, 
i.e. 

``` math
\begin{equation}
$X_A+X_B~Binomial\left(n_A+n_B,p_{pooled}\right)
\hspace{2cm} (8)
\end{equation}
```

**Note**: The binomial coefficients are multiplied rather than pooled, as they originally come from two independent samples, the true $H_0$ is an assumption, not the property of the data. Actually, the choice doesn’t matter, as this term will be zeroed when taking the derivative.

Let’s simplify notation and replace $p_{pooled}$ with just p. Now, the log-likelihood, $log(L(p))=\ell\(p)$, is defined as:

``` math
\begin{equation}
\ell(p)=const\bullet(x_A+x_B)\bullet log(p)+(n_A+n_B-x_A-x_B)\bullet log(1-p)
\hspace{2cm} (9)
\end{equation}
```

(I wrote _const_ to highlight that this term will disappear after taking the derivative).
Now by taking ${\frac{d\ell\left(p\right)}{dp}}$ and setting it to 0, we obtain (provided that $p\notin0,1$):

``` math
\begin{equation}
\begin{aligned}
\frac{d\ell\left(p\right)}{dp}=\frac{x_A+x_B}{p}-\frac{n_A+n_B-x_A-x_B}{1-p}=0, \\
\left(x_A+x_B\right)\left(1-p\right)=\left(n_A+n_B-x_A-x_B\right)\bullet p \\
x_A+x_B={(n}_A+n_B)\bullet p\\
p\ =\ \frac{x_A+x_B}{n_A+n_B}
\end{aligned}
\hspace{2cm} (10)
\end{equation}
```

Or, alternatively, since $p_i={\frac{x_i}{n_i}}$:

``` math
\begin{equation}
p=\frac{p_An_A+p_Bn_B}{n_A+n_B}
\hspace{2cm} (11)
\end{equation}
```

**Log-likelihood of the logistic regression with a single binary predictor**
The log-likelihood function is of the form:

``` math
\begin{equation}
\ell(\beta)=log(L(\beta))=log\left(\prod_{i=1}^{n}p_i^{y_i}\bullet\left(1-p_i\right)^{1-y_i}\right)=\sum_{i=1}^{n}{{[y}_ilog}(p_i)+(1-y_i)log(1-p_i)]
\hspace{2cm} (12)
\end{equation}
```
where β is the vector of estimated coefficients, i.e. $(β_0, β_1)$.

Let’s express p and 1-p in terms of $η=β_0+Xβ_1$:

``` math
\begin{equation}
\begin{cases} 
p=\frac{e^\eta}{1+e^\eta} \\ 
1-p=\frac{1}{1+e^\eta} 
\end{cases}
\hspace{2cm} (13)
\end{equation}
```
Then:
``` math
\begin{equation}
\begin{aligned}
\ell\left(\beta\right)&=\sum_{i=1}^{n}\left[y_ilog\left(\frac{e^{\eta_i}}{1+e^{\eta_i}}\right)+\left(1-y_i\right)log\left(\frac{1}{1+e^{\eta_i}}\right)\right] \\
&=\sum_{i=1}^{n}\left[y_i\left(log\left(e^{\eta_i}\right)-log\left(1+e^{\eta_i}\right)\right)-\left(1-y_i\right)log\left(1+e^{\eta_i}\right)\right] \\
&=\sum_{i=1}^{n}\left[y_i\eta_i\ -y_ilog\left(1+e^{\eta_i}\right)-\left(1-y_i\right)log\left(1+e^{\eta_i}\right)\right] \\
&=\sum_{i=1}^{n}{y_i\eta_i-log\left(1+e^{\eta_i}\right)}
\end{aligned}
\hspace{2cm} (14)
\end{equation}
```

**Gradient, Score, Hessian, Fisher Information, Covariance…**

We will need both the gradient and Hessian of the log-likelihood. For future use, we will call the gradient as _Rao score function_ denoted by “U” and the Hessian as “H”.
First, let’s find the form of U(β):
``` math
\begin{equation}
U\left(\beta\right)=\left[\begin{matrix}\frac{\partial\ell\left(\beta\right)}{\partial\beta_0}\\\frac{\partial\ell\left(\beta\right)}{\partial\beta_1}\\\end{matrix}\right]
\hspace{2cm} (15)
\end{equation}
```

By noticing that ${\frac{\partial\eta_i}{\partial\beta_0}}=1 ,  {\frac{\partial\eta_i}{\partial\beta_1}}=x_i$ and remembering that $p=\frac{e^\eta}{1+e^\eta}$ we obtain:

``` math
\begin{equation}
\frac{\partial\ell\left(\beta\right)}{\partial\beta_0}=\sum_{i=1}^{n}{\left(y_i\frac{\partial\eta_i}{{\partial\beta}_0}-\frac{1}{1+e^{\eta_i}}e^{\eta_i}\frac{\partial\eta_i}{{\partial\beta}_0}\right)=\sum_{i=1}^{n}\left(y_{i\bullet}1-\frac{e^{\eta_i}}{1+e^{\eta_i}}\bullet1\right)=\sum_{i=1}^{n}\left(y_i-p_i\right)}
\hspace{2cm} (16)
\end{equation}
```
and

``` math
\begin{equation}
\frac{\partial\ell\left(\beta\right)}{\partial\beta_1}=\sum_{i=1}^{n}{\left(y_i\frac{\partial\eta_i}{{\partial\beta}_1}-\frac{1}{1+e^{\eta_i}}e^{\eta_i}\frac{\partial\eta_i}{{\partial\beta}_1}\right)=\sum_{i=1}^{n}\left(y_ix_i-\frac{e^{\eta_i}}{1+e^{\eta_i}}x_i\right)=\sum_{i=1}^{n}{x_i\left(y_i-p_i\right)}}
\hspace{2cm} (17)
\end{equation}
```
So finally:
``` math
\begin{equation}
U\left(\beta\right)=\left[\begin{matrix}\sum_{i=1}^{n}\left(y_i-p_i\right)\\\sum_{i=1}^{n}{x_i\left(y_i-p_i\right)}\\\end{matrix}\right]
\hspace{2cm} (18)
\end{equation}
```

Now, the Hessian:
``` math
\begin{equation}
H\left(\beta\right)=\left[\begin{matrix}\frac{\partial^2\ell\left(\beta\right)}{\partial\beta_0^2}&\frac{\partial^2\ell\left(\beta\right)}{\partial\beta_0\beta_1}\\\frac{\partial^2\ell\left(\beta\right)}{\partial\beta_1\beta_0}&\frac{\partial^2\ell\left(\beta\right)}{\partial\beta_1^2}\\\end{matrix}\right]
\hspace{2cm} (19)
\end{equation}
```
The partial derivatives are as follows:
``` math
\begin{equation}
\frac{\partial^2\ell\left(\beta\right)}{\partial\beta_0^2}=\frac{\partial}{\partial\beta_0}\sum_{i=1}^{n}\left(y_i-p_i\right)=\sum_{i=1}^{n}{-\frac{\partial p_i}{\partial\beta_0}}=-\sum_{i=1}^{n}{\frac{\partial p_i}{\partial\beta_0}=-}\sum_{i=1}^{n}{\frac{\partial p_i}{\partial\eta_0}\frac{\partial\eta_0}{\partial\beta_0}=}-\sum_{i=1}^{n}{p_i\left(1-p_i\right)}
\hspace{1cm} (20)
\end{equation}
```
``` math
\begin{equation}
\frac{\partial^2\ell\left(\beta\right)}{\partial\beta_1^2}=\frac{\partial}{\partial\beta_1}\sum_{i=1}^{n}{x_i\left(y_i-p_i\right)}=\sum_{i=1}^{n}{x_i\left(-\frac{\partial p_i}{\partial\beta_1}\right)}=-\sum_{i=1}^{n}{x_ip_i\left(1-p_i\right)}x_i=-\sum_{i=1}^{n}{x_i^2p_i\left(1-p_i\right)}
\hspace{1cm} (21)
\end{equation}
```
``` math
\begin{equation}
\frac{\partial^2\ell\left(\beta\right)}{\partial\beta_0\beta_1}=\frac{\partial^2\ell\left(\beta\right)}{\partial\beta_1\beta_0}=\frac{\partial}{\partial\beta_1}\sum_{i=1}^{n}\left(y_i-p_i\right)=-\sum_{i=1}^{n}\frac{\partial p_i}{\partial\beta_1}=-\sum_{i=1}^{n}{p_i\left(1-p_i\right)}x_i
\hspace{1cm} (22)
\end{equation}
```
Therefore:
``` math
\begin{equation}
H\left(\beta\right)=\left[\begin{matrix}-\sum_{i=1}^{n}{p_i\left(1-p_i\right)}&-\sum_{i=1}^{n}{p_i\left(1-p_i\right)}x_i\\-\sum_{i=1}^{n}{p_i\left(1-p_i\right)}x_i&-\sum_{i=1}^{n}{x_i^2p_i\left(1-p_i\right)}\\\end{matrix}\right]
\hspace{1cm} (23)
\end{equation}
```
Let’s also determine the Fisher Information matrix:
``` math
\begin{equation}
I\left(\beta\right)=-E(H(\beta))=\left[\begin{matrix}\sum_{i=1}^{n}{p_i\left(1-p_i\right)}&\sum_{i=1}^{n}{p_i\left(1-p_i\right)}x_i\\\sum_{i=1}^{n}{p_i\left(1-p_i\right)}x_i&\sum_{i=1}^{n}{x_i^2p_i\left(1-p_i\right)}\\\end{matrix}\right]=\left[\begin{matrix}\sum_{i=1}^{n}{p_iq_i}&\sum_{i=1}^{n}{p_iq_i}x_i\\\sum_{i=1}^{n}{p_iq_i}x_i&\sum_{i=1}^{n}{p_iq_i}x_i^2\\\end{matrix}\right]
\hspace{1cm} (24)
\end{equation}
```
where $q_i=1-p_i$.

This can be further expanded, by substituting sums with appropriate (per respective group, A and B) counts of elements, 
remembering that $n=n_A+n_B$, A: X=0, B: X=1, $p_A=p(Y=1|X=A)$ and $p_B=p(Y=1|X=B)$, and:

``` math
\begin{equation}
\sum_{i=1}^{n}{p_i=}\sum_{i:\ X_i=0} p_i+\sum_{i:\ X_i=1} p_i=\sum_{i=1}^{n_A}p_i+\sum_{i=n_A+1}^{{n_B+n}_B}p_i=n_Ap_A\ +\ n_Bp_B
\hspace{2cm} (25)
\end{equation}
```
So the final, useful form is:
``` math
\begin{equation}
I\left(\beta\right)=\left[\begin{matrix}\sum_{i:\ X_i=0}{p_iq_i}+\sum_{i:\ X_i=1}{p_iq_i}&\sum_{i:\ X_i=1}{1p_iq_i}\\\sum_{i:\ X_i=1}{1p_iq_i}&\sum_{i:\ X_i=1}{1p_iq_i}\\\end{matrix}\right]=\left[\begin{matrix}n_Ap_Aq_A+n_Bp_Bq_B&n_Bp_Bq_B\\n_Bp_Bq_B&n_Bp_Bq_B\\\end{matrix}\right]
\hspace{1cm} (26)
\end{equation}
```
This matrix will be used to find the covariance one:
``` math
\begin{equation}
{I\left(\beta\right)}^{-1}\ =\ \Sigma\left(\beta\right)=\left[\begin{matrix}var(\beta_0)&cov(\beta_0,\beta_1)\\cov(\beta_0,\beta_1&var(\beta_1)\\\end{matrix}\right]
\hspace{2cm} (27)
\end{equation}
```

-----
**Another way to obtain I(β)**
The covariance matrix can also be obtained from:
``` math
\begin{equation}
I\left(\beta\right)=X^TWX
\hspace{2cm} (28)
\end{equation}
```
where X is the design matrix with 2 columns ($β_0$ of 1s and $β_1$ indicating when X=1), with $n_A$ and $n_B$ number of rows corresponding to group A and B, respectively.

``` math
\begin{equation}
X=\left[\begin{matrix}1&A=0\\1&A=0\\\vdots&\vdots\\1&A=0\\1&B=1\\\vdots&\vdots\\1&B=1\\\end{matrix}\right]\ 
\hspace{2cm} (29)
\end{equation}
```
Now, W is the diagonal matrix of weights, of the block-diagonal form:
``` math
\begin{equation}
W=diag\left(p_i\times\left(1-p_i\right)\right)=\left[\begin{matrix}p_Aq_A&0&\ldots&0&0&0&\ldots&0\\0&p_Aq_A&\cdots&0&0&0&\ldots&0\\\vdots&\vdots&\ddots&\vdots&\vdots&\vdots&\ddots&\vdots\\0&0&\ldots&p_Aq_A&0&0&\ldots&0\\0&0&\ldots&0&p_Bq_B&0&\ldots&0\\0&0&\ldots&0&0&p_Bq_B&\ldots&0\\\vdots&\vdots&\ddots&\vdots&\vdots&\vdots&\ddots&\vdots\\0&0&\ldots&0&0&0&\ldots&p_Bq_B\\\end{matrix}\right]
\hspace{2cm} (30)
\end{equation}
```
which can be expressed in a simpler form as:

``` math
\begin{equation}
W=diag\left(p_i\times\left(1-p_i\right)\right)=\left[\begin{matrix}p_Aq_AI_{nA}&0\\0&p_Bq_BI_{nB}\\\end{matrix}\right]
\hspace{2cm} (31)
\end{equation}
```
where $I_nA$ and $I_nB$ are respective identity matrices.

The result of the matrix multiplication can be expressed as appropriate sums:

``` math
\begin{equation}
X^TWX=\left[\begin{matrix}\sum_{i=1}^{n}{1p_iq_i}&\sum_{i=1}^{n}{x_ip_iq_i}\\\sum_{i=1}^{n}{x_ip_iq_i}&\sum_{i=1}^{n}{x_i^2p_iq_i}\\\end{matrix}\right]=I\left(\beta\right)
\hspace{2cm} (32)
\end{equation}
```
where 1 is the result of multiplying 1 x 1 (the $β_0$ vector), and X refer to the other products of the $β_0$ and $β_1$ vectors. 
Notice, that this is **exactly the matrix #24**.

-----
**The Rao Score and Information Matrix under $H_0: β_1=0$**

_The Rao Score_

The Rao Score and Information Matrix under $H0: β_1=0$
Recall the formula #18:

``` math
\begin{equation}
U\left(\beta\right)=\left[\begin{matrix}\sum_{i=1}^{n}\left(y_i-p_i\right)\\\sum_{i=1}^{n}{x_i\left(y_i-p_i\right)}\\\end{matrix}\right]
\hspace{2cm} (18)
\end{equation}
```
Again, $n_A$ and $n_B$ depict the respective number of rows corresponding to group A and B, 
for group A: X=0 and for B: X =1 and $n=n_A+n_B$. Also, under $H_0$, $p_A=p_B=p$ (i.e. $\forall_i{\ p}_i=p$), where p is the pooled probability. 

``` math
\begin{equation}
U\left(\beta_0,\ \beta_1=0\right)=\left[\begin{matrix}\left(\sum_{i=1}^{n}y_i\right)-(n_{A\ }+n_B)\bullet p\\\left(\sum_{i:\ X_i=1} y_i\right){-\ n}_Bp\\\end{matrix}\right]
\hspace{2cm} (33)
\end{equation}
```

Here, $y_i$ is the response vector containing **1**s in respective group. Summing all those **1**s yields the total number of successes in this group, which can be expressed as:

``` math
\begin{equation}
\sum_{i=1}^{n} y_i = \sum_{i:\ X_i\in\{0, 1\}} y_i = n_A p_A + n_B p_B \overset{H_0}{\Rightarrow} (n_A + n_B) \cdot p
\hspace{2cm} (34)
\end{equation}
```
So, the first element of the vector becomes 0 out and finally:
``` math
\begin{equation}
U\left(\beta_0,\ \beta_1=0\right)=\left[\begin{matrix}0\\n_B(p_B-p)\\\end{matrix}\right]
\hspace{2cm} (35)
\end{equation}
```

_Information Matrix_

By recalling **matrix #26** and remembering that under $H_0$, $p_A=p_B=p$ we obtain:
``` math
\begin{equation}
I\left(\beta_0,\ \beta_1=0\right)=\left[\begin{matrix}pq(n_A+n_B)&n_Bpq\\n_Bpq&n_Bpq\\\end{matrix}\right]
\hspace{2cm} (36)
\end{equation}
```
Let’s also calculate $I^{-1}$

``` math
\begin{equation}
\mathrm{\Sigma}={I\left(\beta_0,\ \beta_1=0\right)}^{-1}=\frac{1}{pq(n_A+n_B)\bullet n_Bpq-\left(n_Bpq\right)^2}\left[\begin{matrix}n_Bpq&-n_Bpq\\-n_Bpq&pq(n_A+n_B)\\\end{matrix}\right]
\hspace{1cm} (37)
\end{equation}
```
After simplifying the denominator term:

``` math
\begin{equation}
pq(n_A+n_B)\bullet n_Bpq-\left(n_Bpq\right)^2=p^2q^2n_An_B\ +\ p^2q^2n_B^2-p^2q^2n_B^2\ =p^2q^2n_An_B
\hspace{1cm} (38)
\end{equation}
```
we finally obtain:

``` math
\begin{equation}
\mathrm{\Sigma}(\beta_0,\ \beta_1=0)={I\left(\beta_0,\ \beta_1=0\right)}^{-1}=\frac{1}{p^2q^2n_An_B}\left[\begin{matrix}n_Bpq&-n_Bpq\\-n_Bpq&pq(n_A+n_B)\\\end{matrix}\right]=\left[\begin{matrix}\frac{1}{n_Apq}&-\frac{1}{n_Apq}\\-\frac{1}{n_Apq}&\frac{n_A+n_B}{n_An_Bpq}\\\end{matrix}\right]
\hspace{1cm} (39)
\end{equation}
```

**Rao score test under $H_0: β_1=0$**

The Rao score test (called also Lagrange multiplier test) under $H_0$ is defined as the following quadratic form:

``` math
\begin{equation}
R={U(\beta_0,\ \beta_1=0)}^T\bullet {I(\beta_0,\ \beta_1=0)}^{-1}\bullet U(\beta_0,\ \beta_1=0)
\hspace{2cm} (40)
\end{equation}
```
But since the first element of U is 0, this reduces to just scalar operation: $U^2/I = U^2Σ$:

``` math
\begin{equation}
\begin{aligned}
R&=U\left(\beta_0,\beta_1=0\right)\mathrm{\Sigma}\left(\beta_0,\beta_1=0\right)= \\
&=\left[n_B\left(p_B-p\right)\right]^2\bullet\frac{n_A+n_B}{pqn_An_B}=\frac{\left[n_B\left(p_B-p\right)\right]^2}{pq\frac{n_An_B}{n_A+n_B}} \\
&=\frac{n_B^2\left(p_B-\frac{n_Ap_A+n_Bp_B}{n_A+n_B}\right)^2}{pq\frac{n_An_B}{n_A+n_B}}=\frac{{n_B^2\left(\frac{p_Bn_A+p_Bn_B-n_Ap_A-n_Bp_B}{n_A+n_B}\right)}^2}{pq\frac{n_An_B}{n_A+n_B}} \\
&=\frac{\frac{n_B^2\left(p_Bn_A-n_Ap_A\right)^2}{\left(n_A+n_B\right)^2}}{pq\frac{n_An_B}{n_A+n_B}}=\frac{\frac{n_B^2n_A^2\left(p_B-p_A\right)^2}{\left(n_A+n_B\right)^2}}{pq\frac{n_An_B}{n_A+n_B}} \\
&=\frac{\left(p_B-p_A\right)^2}{pq\frac{n_An_B}{n_A+n_B}\frac{\left(n_A+n_B\right)^2}{n_B^2n_A^2}}=\frac{\left(p_B-p_A\right)^2}{pq\frac{n_A+n_B}{n_An_B}}=\frac{{{(p}_B-p_A)}^2}{pq\left(\frac{1}{n_B}+\frac{1}{n_A}\right)} \\
&=\frac{{{(p}_B-p_A)}^2}{p\left(1-p\right)\left(\frac{1}{n_A}+\frac{1}{n_B}\right)}=z^2\blacksquare
\end{aligned}
\hspace{2cm} (41)
\end{equation}
```

This way I have proven the equivalence of the 2-sample Wald’s z-statistic for comparing proportions with pooled variances and the Rao score test over the logistic regression with a single binary predictor distinguishing the compared samples.

``` r
> wald_z_test_pooled(x1 = 6, n1 = 20, x2 = 10, n2 = 20)
  diff         z     chi2        se   p.value p.value_1        LCI       HCI
1 -0.2 -1.290994 1.666667 0.1549193 0.1967056 0.0983528 -0.5036363 0.1036363

> prop.test(c(6,10), c(20,20), correct = FALSE)

	2-sample test for equality of proportions without continuity correction

data:  c(6, 10) out of c(20, 20)
X-squared = 1.6667, df = 1, p-value = 0.1967
alternative hypothesis: two.sided
95 percent confidence interval:
 -0.49724326  0.09724326
sample estimates:
prop 1 prop 2 
   0.3    0.5 

> data <- data.frame(response = factor(c(rep("Success", 6), rep("Failure", 20-6),
+                                        rep("Success", 10), rep("Failure", 20-10))),
+                    grp    = factor(rep(c("B", "A"), each=20)))

> m <- glm(response ~ grp, data = data, family = binomial(link = "logit"))

> anova(m, test = "Rao")
Analysis of Deviance Table

Model: binomial, link: logit
Response: response
Terms added sequentially (first to last)

     Df Deviance Resid. Df Resid. Dev    Rao Pr(>Chi)
NULL                    39     53.841                
grp   1   1.6805        38     52.160 1.6667   0.1967
```

![obraz](https://github.com/user-attachments/assets/183ecccc-14d1-4596-bfb8-0aa174621fee)

The implementation of the z statistic:
``` r
wald_z_test_pooled <- function(x1, n1, x2, n2, conf.level=0.95) {
      p1 <- x1/n1
  p2 <- x2/n2
  
  p_pool <- (p1*n1 + p2*n2) / (n1+n2)
  se_p1 <- sqrt(p_pool * (1 - p_pool) / n1)
  se_p2 <- sqrt(p_pool * (1 - p_pool) / n2)
  
  se_diff <- sqrt(se_p1^2 + se_p2^2)
  
  z <- (p1 - p2) / se_diff
  p <- 2 * (1 - pnorm(abs(z)))
  hCI <- abs(qnorm((1 - conf.level)/2)) * se_diff
  
  return(data.frame(diff=p1-p2, 
                    z = z, 
                    chi2 = z^2,
                    se = se_diff, p.value = p, p.value_1 =p/2,
                    LCI = (p1-p2) - hCI,
                    HCI = (p1-p2) + hCI,
                    row.names = NULL))
}
```
