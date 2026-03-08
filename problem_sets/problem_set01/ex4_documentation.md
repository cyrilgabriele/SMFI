# Exercise 4 — Custom Derivative Product under GBM and Risk-Neutral Pricing

## Prompt used

Create an innovative derivative product beyond vanilla options, but still economically sensible and consistent with the GBM / Black–Scholes framework. Then provide an intuitive description of the product, its formal payoff function, and its price under the risk-neutral measure $\mathbb{Q}$. Pick your own parameters and implement the pricing in R.

------------------------------------------------------------------------

## Product idea: Protected Target-Range Bull-Spread Note

I propose the following structured equity-linked derivative:

A **Protected Target-Range Bull-Spread Note** is a maturity payoff product with three components:

1.  **Principal protection:** the investor receives a guaranteed amount at maturity.
2.  **Capped upside participation:** the investor participates in the stock upside, but only up to a cap.
3.  **Target-range bonus:** the investor receives an additional fixed bonus if the stock finishes in a predefined “good” range.

This is more innovative than a vanilla call or put, but still fully sensible from a financial engineering perspective.

------------------------------------------------------------------------

## Intuitive description

This product is suitable for an investor who is **moderately bullish** on the stock but does not want to bear downside equity risk.

-   If the stock performs poorly, the investor still receives the guaranteed maturity amount.
-   If the stock rises moderately, the investor participates in the upside.
-   If the stock ends in a favorable range, the investor receives an extra fixed bonus.
-   If the stock increases too strongly above the cap, the investor no longer benefits from further gains. This capped upside helps finance the guarantee and the bonus.

Hence, the product rewards **moderate positive performance** rather than either very poor or extremely high outcomes.

------------------------------------------------------------------------

## Model assumptions

Assume a non-dividend-paying stock under the risk-neutral measure $\mathbb{Q}$ follows a geometric Brownian motion:

$$
dS_t = r S_t\,dt + \sigma S_t\,dW_t^{\mathbb{Q}},
$$

with solution

$$
S_T = S_0 \exp\left(\left(r - \frac{1}{2}\sigma^2\right)T + \sigma \sqrt{T} Z\right),
\qquad Z \sim N(0,1).
$$

The derivative is priced by risk-neutral valuation:

$$
V_0 = e^{-rT}\mathbb{E}^{\mathbb{Q}}[X_T].
$$

------------------------------------------------------------------------

## Choice of parameters

I choose the following parameters:

$$
S_0 = 100, \qquad r = 0.02, \qquad \sigma = 0.22, \qquad T = 2.
$$

For the product design:

-   guaranteed amount: $G = 100$
-   participation rate: $\lambda = 0.70$
-   lower strike of upside participation: $K = 100$
-   upper strike / cap: $U = 130$
-   target range lower bound: $L = 95$
-   target range upper bound: $U = 130$
-   bonus amount: $B = 6$

------------------------------------------------------------------------

## Formal payoff function

Let $X_T$ denote the payoff at maturity. Then

$$
X_T
=
G
+
\lambda\left[(S_T-K)^+-(S_T-U)^+\right]
+
B\,\mathbf{1}_{\{L \le S_T \le U\}}.
$$

With the chosen parameters, this becomes

$$
X_T
=
100
+
0.70\left[(S_T-100)^+-(S_T-130)^+\right]
+
6\,\mathbf{1}_{\{95 \le S_T \le 130\}}.
$$

------------------------------------------------------------------------

## Interpretation of the payoff

The payoff can be decomposed into:

1.  **Zero-coupon bond component:**\
    $$
    100
    $$ at maturity.

2.  **Bull call spread with participation factor** $0.70$:\
    $$
    0.70\left[(S_T-100)^+-(S_T-130)^+\right]
    $$ which gives upside exposure between 100 and 130.

3.  **Range digital bonus:**\
    $$
    6\,\mathbf{1}_{\{95 \le S_T \le 130\}}
    $$ which pays an additional fixed amount if the terminal stock price lies in the target range.

------------------------------------------------------------------------

## Pricing formula

By linearity of expectation,

$$
V_0
=
100e^{-rT}
+
0.70\Big(C_0(100)-C_0(130)\Big)
+
6e^{-rT}\mathbb{Q}(95 \le S_T \le 130),
$$

where $C_0(K)$ is the Black–Scholes call price with strike $K$:

$$
C_0(K)=S_0\Phi(d_1)-Ke^{-rT}\Phi(d_2),
$$

with

$$
d_1=\frac{\ln(S_0/K)+(r+\tfrac12\sigma^2)T}{\sigma\sqrt{T}},
\qquad
d_2=d_1-\sigma\sqrt{T}.
$$

For the corridor probability,

$$
\mathbb{Q}(95 \le S_T \le 130)
=
\Phi(d_2(95))-\Phi(d_2(130)),
$$

where

$$
d_2(x)=\frac{\ln(S_0/x)+(r-\tfrac12\sigma^2)T}{\sigma\sqrt{T}}.
$$

Therefore the analytical price is

$$
\boxed{
V_0
=
100e^{-rT}
+
0.70\Big(C_0(100)-C_0(130)\Big)
+
6e^{-rT}\Big[\Phi(d_2(95))-\Phi(d_2(130))\Big]
}
$$

------------------------------------------------------------------------

## Numerical price

Using

$$
S_0 = 100, \quad r = 0.02, \quad \sigma = 0.22, \quad T = 2,
$$

the numerical value is approximately

$$
\boxed{V_0 \approx 104.73.}
$$

A Monte Carlo simulation under $\mathbb{Q}$ confirms the analytical result.

------------------------------------------------------------------------

## Conclusion

The **Protected Target-Range Bull-Spread Note** is a non-vanilla but economically meaningful derivative product. It combines:

-   capital protection,
-   capped equity upside participation,
-   and a conditional range bonus.

Under the GBM model and risk-neutral pricing, the fair price of the product is approximately

$$
\boxed{104.73.}
$$

This product is a good example of financial engineering beyond standard vanilla options while remaining fully consistent with the Black–Scholes framework.
