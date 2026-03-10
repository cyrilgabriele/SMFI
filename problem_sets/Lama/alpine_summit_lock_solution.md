# Alpine Summit-Lock Asian Knock-Out Call

## Stochastic Modelling — Exotic Derivative Design & Monte Carlo Pricing

---

## 1. The Prompt

This solution was generated in response to a detailed Master's-level assignment prompt requesting: (i) the design of an innovative path-dependent exotic equity option combining at least two distinct exotic features, (ii) a formal payoff specification under the risk-neutral measure $\mathbb{Q}$, (iii) an intuitive investor-facing description, and (iv) a fully annotated Monte Carlo pricing implementation in R using exact GBM simulation with $N = 100{,}000$ paths.

## Full Prompt:
Role: Act as a Senior Quantitative Analyst and Financial Engineer at a top-tier Swiss investment bank. You possess expert-level knowledge in stochastic calculus, risk-neutral pricing, and exotic derivative structuring.
Task: I am taking a Master's level course in Stochastic Modelling. I need to design an innovative equity derivative product (strictly beyond standard vanilla options) and price it using Monte Carlo simulation under the risk-neutral measure $\mathbb{Q}$.
Product Requirements:
1. Innovation: Design a realistic but highly creative path-dependent exotic option. It should combine at least two distinct exotic features (e.g., an Asian-style time-weighted average combined with a Knock-Out barrier, or a Lookback feature with a conditional step-up guarantee). Give it a catchy, professional market name (e.g., "Yield-Enhanced Asian Knock-Out" or "Memory-Ratchet Barrier Call").
2. Underlying Dynamics: The underlying stock follows a Geometric Brownian Motion (GBM). Switch to the risk-neutral measure $\mathbb{Q}$, where the drift is the continuous risk-free rate $r$. The SDE is:
$$d\tilde{S}_{t} = r\tilde{S}_{t}dt + \sigma\tilde{S}_{t}d\tilde{W}_{t}^{\mathbb{Q}}$$
3. Parameters: Define realistic market parameters. Pick a starting price $S_0$, a continuous risk-free rate $r$ (e.g., 2%), a volatility $\sigma$ (e.g., 20%), and a maturity $T$ (e.g., 1 or 2 years). Define specific strike(s) or barrier levels required for your product.
Required Output Format (Please strictly follow these sections):
1. The Prompt(s): Briefly state that this entire instruction was the prompt used to generate the solution.
2. Intuitive Description: Write a 2-3 paragraph pitch for this structured product. Explain to a potential investor how the product behaves in different market scenarios (bull, bear, sideways) and what specific risk/reward profile it caters to.
3. Formal Payoff Function: Provide the exact mathematical definition of the maturity payoff $\tilde{P}_T$. Use proper stochastic calculus notation (e.g., integrals for continuous averages or summations for discrete time steps, indicator functions for barriers).
4. Pricing and R Code: Provide highly readable, fully annotated, and executable R code.
* The code must simulate $N = 100,000$ paths using the exact closed-form GBM solution (not Euler discretization) under $\mathbb{Q}$.
* Calculate the path-dependent payoff for each simulation.
* Discount the expected payoff back to $t=0$ to find the option price: $V_0 = \exp(-rT) \cdot \mathbb{E}^{\mathbb{Q}}[\tilde{P}_T]$.
* Include a set.seed(42) for reproducibility.
* Print the final calculated price clearly.

---

## 2. Intuitive Description — Investor Pitch

**The Alpine Summit-Lock Asian Knock-Out Call** is a 2-year equity-linked structured note designed for investors who are moderately bullish on the underlying stock but want built-in downside discipline and protection against late-cycle reversals. The product innovates by fusing three distinct exotic mechanisms into a single coherent payoff: a *memory ratchet lock-in*, an *Asian tail settlement*, and a *discrete down-and-out barrier*.

During the first 21 months — the **accumulation phase** — the product continuously monitors the stock's monthly closing price. Each time the stock trades above the strike, the investor automatically "locks in" a guaranteed floor equal to 50% of the in-the-money amount at that date. The product *remembers* the best lock-in achieved across all accumulation dates. This is the memory ratchet: if the stock rallies to 130 in month 8 but subsequently retreats, the investor retains a guaranteed floor payoff of $0.50 \times (130 - 100) = 15$, regardless of where the stock finishes. The final three months — the **settlement phase** — use an Asian-style arithmetic average of the stock price, which determines the conventional call component. At maturity the investor receives the *greater* of the Asian call payoff and the locked-in floor — provided the barrier has never been breached.

**Bull market**: The product excels. Both the Asian average and the lock-in floor are deep in the money. The payoff closely tracks a vanilla call, with the added cushion that a late pullback is mitigated by the lock-in memory. **Sideways market**: The lock-in activates opportunistically during transient rallies, potentially delivering a modest payoff even if the stock finishes near its starting level. **Bear market**: If the decline is gradual, the payoff is zero but capital at risk is limited to the premium. If the decline is severe (stock falls through the 70-level barrier at any monitoring date), the option is knocked out entirely with zero payoff — this is the trade-off that cheapens the structure and funds the ratchet feature. The product is unsuitable for investors expecting a deep correction or a crash.

---

## 3. Formal Payoff Function

### 3.1 Setup and Notation

Let $\tilde{S}_t$ denote the stock price process under the risk-neutral measure $\mathbb{Q}$, following Geometric Brownian Motion:

$$d\tilde{S}_t = r\,\tilde{S}_t\,dt + \sigma\,\tilde{S}_t\,d\tilde{W}_t^{\mathbb{Q}}$$

with exact closed-form solution:

$$\tilde{S}_t = S_0 \exp\!\Bigl[\Bigl(r - \tfrac{\sigma^2}{2}\Bigr)t + \sigma\,\tilde{W}_t^{\mathbb{Q}}\Bigr]$$

Define $M = 24$ equally-spaced discrete monitoring dates:

$$t_i = \frac{i \cdot T}{M}, \qquad i = 1, 2, \ldots, 24$$

where $T = 2$ years is the maturity. The monitoring dates are partitioned into:

- **Accumulation phase**: $\mathcal{A} = \{t_1, t_2, \ldots, t_{21}\}$ (months 1--21)
- **Settlement phase**: $\mathcal{S} = \{t_{22}, t_{23}, t_{24}\}$ (months 22--24)

### 3.2 Product Parameters

| Symbol | Value | Description |
|--------|-------|-------------|
| $S_0$ | 100 | Initial stock price |
| $K$ | 100 | Strike price (ATM) |
| $B$ | 70 | Down-and-out barrier |
| $\beta$ | 0.50 | Lock-in fraction |
| $r$ | 0.02 | Continuous risk-free rate |
| $\sigma$ | 0.20 | Annualised volatility |
| $T$ | 2 | Maturity (years) |

### 3.3 Payoff Components

**Component 1 — Down-and-Out Barrier Indicator (discrete monitoring):**

$$\mathbb{1}_{\text{survived}} = \prod_{i=1}^{24} \mathbb{1}\!\bigl\{\tilde{S}_{t_i} > B\bigr\} = \mathbb{1}\!\Bigl\{\min_{1 \le i \le 24}\, \tilde{S}_{t_i} > B\Bigr\}$$

**Component 2 — Summit Lock-In Floor (memory ratchet over accumulation phase):**

$$F = \max_{i \in \{1, \ldots, 21\}}\; \beta \cdot \bigl(\tilde{S}_{t_i} - K\bigr)^+$$

where $(x)^+ = \max(x, 0)$. The floor $F$ captures and retains the best partial intrinsic value achieved during the accumulation phase.

**Component 3 — Asian Tail Average (settlement phase):**

$$A = \frac{1}{3} \sum_{j=22}^{24} \tilde{S}_{t_j}$$

**Component 4 — Asian Call Payoff:**

$$C = (A - K)^+$$

### 3.4 Maturity Payoff

The maturity payoff of the Alpine Summit-Lock Asian Knock-Out Call is:

$$\boxed{\tilde{P}_T = \mathbb{1}_{\text{survived}} \cdot \max\!\bigl(C,\; F\bigr) = \mathbb{1}\!\Bigl\{\min_{1 \le i \le 24}\, \tilde{S}_{t_i} > B\Bigr\} \cdot \max\!\Bigl(\bigl(A - K\bigr)^+,\; \max_{1 \le i \le 21}\; \beta\bigl(\tilde{S}_{t_i} - K\bigr)^+\Bigr)}$$

### 3.5 Present Value (Risk-Neutral Price)

The time-0 fair value is the discounted risk-neutral expectation:

$$V_0 = e^{-rT}\;\mathbb{E}^{\mathbb{Q}}\!\bigl[\tilde{P}_T\bigr]$$

estimated via Monte Carlo as:

$$\hat{V}_0 = e^{-rT} \cdot \frac{1}{N}\sum_{n=1}^{N} \tilde{P}_T^{(n)}$$

with standard error $\;\hat{\sigma}_{\hat{V}} = \dfrac{e^{-rT}\,\hat{\sigma}_P}{\sqrt{N}}$.

---

## 4. Monte Carlo Pricing Results

With `set.seed(42)` and $N = 100{,}000$ simulated paths:

| Metric | Value |
|--------|-------|
| **Monte Carlo Price $\hat{V}_0$** | **15.1865** |
| Standard Error | 0.0605 |
| 95% Confidence Interval | [15.0678, 15.3051] |
| Barrier breach rate | 16.9% |
| Lock-in floor dominates payoff | 43.1% of paths |
| Asian payoff dominates | 35.5% of paths |
| Black-Scholes vanilla call benchmark | 13.0957 |
| Exotic / Vanilla price ratio | 115.97% |

The exotic trades at roughly a 16% premium over the vanilla European call, reflecting the value of the embedded memory ratchet. This premium is partially offset by the knock-out barrier, which extinguishes the option in approximately 17% of scenarios.

---

## 5. Diagnostic Visualisations

The R script generates six diagnostic plots that illuminate the structure's behaviour and pricing characteristics:

### Plot 1 — Simulated Price Paths

A sample of 50 Monte Carlo paths is drawn over the 2-year horizon, overlaid with the barrier at $B = 70$ and strike at $K = 100$. Paths that breach the barrier are highlighted in red; surviving paths appear in blue. A shaded region marks the 3-month Asian settlement phase (months 22--24), and a vertical marker separates the accumulation and settlement regimes. This plot gives an immediate visual sense of path dispersion, barrier proximity, and the frequency of knock-out events.

### Plot 2 — Payoff Distribution

A histogram with kernel density overlay shows the distribution of positive payoffs at maturity (conditional on the payoff being non-zero). The expected payoff is marked, and a legend reports the fraction of paths that produce positive versus zero payoffs. The pronounced right skew reflects the optionality embedded in the max operator over the Asian and lock-in components.

### Plot 3 — Monte Carlo Convergence

The running estimate of $\hat{V}_0$ is plotted against the number of simulated paths on a logarithmic x-axis, together with a 95% confidence band. This demonstrates how the estimate stabilises as the sample size grows and verifies that $N = 100{,}000$ paths yield adequate precision (standard error below 0.07).

### Plot 4 — Payoff Source Decomposition

A bar chart decomposes all $N$ paths into four mutually exclusive categories: paths where the lock-in floor dominates, paths where the Asian payoff dominates, surviving paths with zero payoff (out-of-the-money), and knocked-out paths. This visualisation quantifies the relative importance of each product feature and shows that the lock-in ratchet is the dominant payoff source in a plurality of scenarios.

### Plot 5 — Volatility Sensitivity

The exotic option price is re-estimated across a grid of volatilities from 10% to 45% (using 20,000 paths per point), plotted alongside the analytical Black-Scholes vanilla call price for the same grid. The calibrated volatility $\sigma = 20\%$ is marked. This reveals the exotic's non-monotonic sensitivity to volatility: higher volatility increases both the upside optionality and the probability of barrier breach, creating a flattening effect relative to the vanilla call at elevated vol levels.

### Plot 6 — Asian Payoff vs Lock-In Floor Scatter

A scatter plot of the Asian tail payoff against the best lock-in floor for a subsample of surviving paths. Points are coloured by which component dominates the final payoff, with a 45-degree reference line. This reveals the substitution dynamics between the two payoff mechanisms: the lock-in floor provides significant value precisely in scenarios where the Asian average underperforms (late-cycle reversals), confirming the product's design rationale.

---

## 6. R Code

The complete, executable R script is provided in the accompanying file **`alpine_summit_lock_pricer.R`**. Key implementation highlights:

- Uses the **exact closed-form GBM** transition $\tilde{S}_{t_{i}} = \tilde{S}_{t_{i-1}} \exp\bigl((r - \tfrac{\sigma^2}{2})\Delta t + \sigma\sqrt{\Delta t}\,Z_i\bigr)$ — not Euler discretisation.
- Fully vectorised matrix operations (no per-path loops) for computational efficiency.
- Includes a payoff decomposition analysis showing how often the lock-in vs. Asian component dominates.
- Benchmarks the exotic price against the analytical Black-Scholes European call.
- Generates **six diagnostic plots**: sample paths, payoff distribution, MC convergence, payoff decomposition bar chart, volatility sensitivity curve, and Asian-vs-lock-in scatter.
- All Unicode box-drawing characters replaced with ASCII equivalents for full compatibility with `pdflatex` when compiling via Quarto/RMarkdown.
