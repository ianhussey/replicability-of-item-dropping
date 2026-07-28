# Estimands for a Monte Carlo study of item removal by "alpha-if-item-removed"

Scope: replicability of dropping items from multi-item
scales using the Cronbach's-alpha-if-item-removed rule. 



---

## 4. Replicability of the selection

### 4a. Selection distribution

$$
q_j = \Pr(J^{*} = j)\quad(\text{plus } q_\varnothing = \Pr(\text{no drop}) \text{ for the conditional rule}),
\qquad \hat q_j = \frac{1}{R}\sum_{r=1}^{R}\mathbf 1\{J^{*}_r = j\}.
$$

Estimated from a **single-sample stream** (one sample per replicate). The whole replicability story follows from this vector.

### 4b. Same-item agreement (two independent samples)

$$
A = \Pr\big(J^{{*}(1)} = J^{{*}(2)}\big) = \sum_{j} q_j^2
$$

(count "both no-drop" separately: report $A$ with and without adding $q_\varnothing^2$).

- **Does this need a second simulation / two samples?** **No, not for $A$.** Because the two samples are i.i.d. from the same DGP, $A$ is a deterministic function of the marginal $q$ vector: $A = \sum_j q_j^2$. You compute it from the single-sample selection distribution — no pairing required. Pairing is only needed for the out-of-sample estimand in §5; if you build that paired design anyway, it *also* yields a direct empirical estimate of $A$ (the fraction of replicates where discovery and replication drop the same item), which should match $\sum_j \hat q_j^2$ and serves as an internal consistency check. Either way it is the **same DGP and same rule** — not a distinct simulation, just (optionally) two draws per replicate at 2$\times$ sampling cost.

### 4c. Noise floor and normalization

Under an **exchangeable DGP** (all items identical: equal $\lambda$ **and** equal $\psi$, i.e. compound symmetry), symmetry forces $q_j = 1/k$ and hence
$$A_{\text{floor}} = \sum_j (1/k)^2 = 1/k.$$

Report $A$ **relative to $1/k$**, never against 1. Chance-corrected index against the exchangeable floor:

$$
A^{*} = \frac{A - 1/k}{1 - 1/k} \in [0,1], \qquad A^{*} = 0 \text{ (uniform)}, \quad A^{*} = 1 \text{ (deterministic)}.
$$

### 4d. Why NOT Cohen's kappa

Kappa corrects observed agreement against the **marginal-product** chance model. For two i.i.d. samples from one fixed DGP, observed agreement and that chance model are the *same number*:
$$
p_o = \sum_j q_j^2, \qquad p_e^{(\text{marginals})} = \sum_j q_j^2 \ \Rightarrow\ \kappa = \frac{p_o - p_e}{1 - p_e} \approx 0
$$
regardless of how concentrated selection is. Kappa deliberately subtracts off the concentration that *is* the replicability signal, so it is structurally near zero and uninformative here. The correct chance baseline is the **uniform/exchangeable floor** $1/k$ (giving $A^{*}$ above), not the empirical marginals. Use raw agreement vs $1/k$, or $A^{*}$.

### 4e. Exchangeable validation cell (mandatory calibration)

Include a condition with **equal loadings and equal unique variances** (compound symmetry). It must return $\hat q_j \approx 1/k$ and $A \approx 1/k$. If it does not, there is a bug. Note that equal loadings with *unequal* errors is **not** exchangeable and will not sit at the floor. Then increase the loading gap / Beta spread and watch $A$ (and $A^{*}$) climb from the floor toward 1 as a function of $n$ and heterogeneity.

### 4f. Set-valued rule (drop-all-that-improve)

Decision is a binary vector $D \in \{0,1\}^k$. Per-item drop probability and agreement:
$$
p_j = \Pr(D_j = 1) = \mathbb E[D_j], \qquad a_j = \Pr\big(D_j^{(1)} = D_j^{(2)}\big) = p_j^2 + (1-p_j)^2.
$$
The same kappa degeneracy applies per item, so do **not** use kappa. Recommendations:

- Report the **per-item drop-probability vector** $p = (p_1,\dots,p_k)$ directly (raw agreement $a_j$ is inflated by "both keep it" for rarely-dropped items).
- For a scalar summary, use agreement **conditional on both samples dropping at least one item**, and/or compare each $a_j$ to its value in the exchangeable cell.
- In the exchangeable cell, drops under the conditional rule are noise-driven and rare; report conditional-on-dropping agreement there too (it should sit at the floor).

---

## 5. Does the improvement replicate out of sample?

Requires **paired** discovery + replication samples per replicate, both i.i.d. from the same DGP.

1. Select $J^{*}$ in the **discovery** sample.
2. Discovery gain (positive by selection): $G_D = \hat\alpha_D(-J^{*}) - \hat\alpha_D(\text{full})$.
3. **Hold $J^{*}$ fixed** and recompute in the **replication** sample:
   $G_R = \hat\alpha_R(-J^{*}) - \hat\alpha_R(\text{full})$.

Estimands:

$$
\text{Shrinkage / selection optimism (headline)}:\quad \mathbb E\big[\,G_D - G_R \mid \text{drop in } D\,\big] > 0,
$$
$$
\mathbb E[G_R] \approx \mathbb E[\Delta_{\text{pop}}]\ \text{(the "honest" gain; $J^{*}$ fixed} \Rightarrow \hat\alpha_R \text{ unbiased for } \alpha_{\text{pop}} \text{ of that fixed set)}.
$$

- **Answer to "is it bias between in and out of sample?"** Yes — the primary quantity is the **in-vs-out difference** $\mathbb E[G_D - G_R]$, i.e. the winner's-curse portion of the observed gain, and it is regression-to-the-mean made concrete. Report $G_D$, $G_R$, and their difference. $G_R$ can be $\le 0$ on average even when $G_D > 0$: the drop looked good, replicates to nothing or to harm. Optionally compare $G_R$ to the true-reliability gain $\Delta_\rho$ to show that even the honest alpha gain need not reflect a reliability gain.
- **Keep distinct from §4:** here $J^{*}$ is *held fixed* from discovery and only re-*evaluated* in replication. Letting the replication sample *re-select* its own item is the agreement question (§4), a different estimand. Do not conflate them; one paired design supports both (agreement uses each sample's own selection; replication uses the discovery selection evaluated in the replication sample).

---

## Summary table

| # | Estimand | Symbol | Level | Sign / role |
|---|----------|--------|-------|-------------|
|  |                         |                                                              |       |                                            |
|      |                         |                                                              |       |                                            |
|      |                         |                                                              |       |                                            |
|      |                         |                                                              |       |                                            |
|      |                         |                                                              |       |                                            |
| 4 | Same-item agreement | $A = \sum_j q_j^2$; floor $1/k$; $A^{*} = \frac{A-1/k}{1-1/k}$ | — | replicability vs uniform floor (not kappa) |
| 5 | Out-of-sample shrinkage | $\mathbb E[G_D - G_R \mid \text{drop}]$ | — | $> 0$; regression to the mean |

---

