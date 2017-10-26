## Frailty model simulation strategy

### Inducing positive correlation among ($R_0$, $R_1$)

A subject-specific frailty parameter $\gamma_i$ captures literal *frailty* that leads to both earlier death and the need for rehospitalization, regardless of treatment status. When $\gamma$ is set to zero, there is no induced correlation (i.e., $\gamma_i = 1 \ \forall \ i$). Larger $\gamma$ means that there is greater correlation among $(R_0, R_1, T_0, T_1)$.

\[ \gamma_i \sim \mathrm{Gamma}(\theta^{-1}, \theta^{-1}) \]

### Rehospitalization ($R$ model)

$R_{z}$, the potential rehospitalization time under $Z=z$, is a mixed distribution. When rehospitalization time is truncated by death (i.e., death occurs without rehospitalization), $R_{z}$ is set to a non-real number, $\bar{\mathbb{R}}$. (For the purposes of operators such as $\mathrm{min}(\cdot)$ and $\mathrm{max}(\cdot)$, one can think of $\bar{\mathbb{R}}$ as $\infty$.)

A $\mathrm{Weibull}(\alpha_1, \kappa_1)$ model is adopted for the hazard of rehospitalization:

\[ h_{1}(r_{zi} | \gamma_i, z) = \gamma_i \alpha_1 \kappa_1 r_{zi}^{\alpha_1 - 1} \mathrm{exp} \left\lbrace \beta_1 z \right\rbrace \text{ for } r_{zi} > 0 \]

\[ h_{1}(r_{zi} | \gamma_i, z) = \gamma_i \alpha_1 r_{zi}^{\alpha_1 - 1} \mathrm{exp} \left\lbrace \beta_{01} + \beta_1 z \right\rbrace \text{ for } r_{zi} > 0 \]


### Death ($T$ model), without prior rehospitalization

$T_{z}$ is the potential death time under $Z=z$. We adopt a $\mathrm{Weibull}(\alpha_2, \kappa_2)$ model for the hazard of death:

\[ h_{2}(t_{zi} | \gamma_i, z) = \gamma_i \alpha_2 \kappa_2 t_{zi}^{\alpha_2 - 1} \mathrm{exp} \left\lbrace \beta_2 z \right\rbrace \text{ for } t_{zi} > 0 \]

\[ h_{2}(t_{zi} | \gamma_i, z) = \gamma_i \alpha_2 t_{zi}^{\alpha_2 - 1} \mathrm{exp} \left\lbrace \beta_{02} + \beta_2 z \right\rbrace \text{ for } t_{zi} > 0 \]

### Death ($T$ model), after rehospitalization

Rehospitalization may change the hazard of death. We use a $\mathrm{Weibull}(\alpha_3, \kappa_3)$ model with a Semi-Markov formulation; that is, the hazard of death at time $t > r_{zi}$ depends only on $(t-r_{zi})$. Concretely, two patients alive at 89 days that were rehospitalized at 15 and 30 days, respectively, will have different baseline hazards of death at Day 90. However, two patients with identical covariate patterns will experience the same baseline hazard of death 10 days after rehospitalization, regardless of whether the rehospitalization occurs at Day 15 or Day 30.

\[ h_{3}(t_{zi} | \gamma_i, r_{zi}, z) = \gamma_i \alpha_3 \kappa_3 (t_{zi}-r_{zi})^{\alpha_3 - 1} \mathrm{exp} \left\lbrace \beta_3 z \right\rbrace \text{ for } t_{zi} > r_{zi} > 0 \]

\[ h_{3}(t_{zi} | \gamma_i, r_{zi}, z) = \gamma_i \alpha_3 (t_{zi}-r_{zi})^{\alpha_3 - 1} \mathrm{exp} \left\lbrace \beta_{03} + \beta_3 z \right\rbrace \text{ for } t_{zi} > r_{zi} > 0 \]
