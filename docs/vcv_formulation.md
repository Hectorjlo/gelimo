# VCV Formulation and PGLS Mathematics in `gelimo`

## 1. Phylogenetic Variance-Covariance Matrix (VCV)
Given a phylogenetic tree $T$ with $N$ species (tips), the VCV matrix $V \in \mathbb{R}^{N \times N}$ represents shared evolutionary history under Brownian Motion:
* $V_{ii} = d(\text{root}, \text{tip}_i)$ (height/depth of tip $i$ from root)
* $V_{ij} = d(\text{root}, \text{MRCA}(i, j))$ (distance from root to Most Recent Common Ancestor of tip $i$ and tip $j$)

## 2. Pagel's Transformations

### Pagel's $\lambda$ (Lambda)
Rescales off-diagonal elements by $\lambda \in [0, 1]$:
$$V(\lambda)_{ij} = \begin{cases} V_{ii} & \text{if } i = j \\ \lambda \cdot V_{ij} & \text{if } i \neq j \end{cases}$$

### Pagel's $\kappa$ (Kappa)
Rescales branch lengths $b$ by $b^\kappa$.

### Pagel's $\delta$ (Delta)
Rescales node heights $h$ by $h^\delta$.

## 3. GLS / PGLS Solver via Cholesky Decomposition
Given response vector $Y \in \mathbb{R}^{N \times 1}$, design matrix $X \in \mathbb{R}^{N \times p}$, and covariance matrix $V$:

1. Compute lower triangular Cholesky factor $L$ such that $V = L L^T$.
2. Transform $X^* = L^{-1} X$ and $Y^* = L^{-1} Y$ using forward substitution.
3. Compute OLS coefficients on transformed system:
   $$\hat{\beta} = (X^{*T} X^*)^{-1} X^{*T} Y^*$$
4. Compute residual variance $\hat{\sigma}^2$ and Log-Likelihood:
   $$\ln L = -\frac{N}{2} \ln(2\pi \hat{\sigma}^2) - \frac{1}{2} \ln |V| - \frac{N}{2}$$
