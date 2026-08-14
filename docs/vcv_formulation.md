# VCV Formulation, PGLS Mathematics and High-Performance Optimizations in `gelimo`

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

## 3. Classic GLS / PGLS Solver via Cholesky Decomposition ($O(N^3)$)
Given response vector $Y \in \mathbb{R}^{N \times 1}$, design matrix $X \in \mathbb{R}^{N \times p}$, and covariance matrix $V$:

1. Compute lower triangular Cholesky factor $L$ such that $V = L L^T$.
2. Transform $X^* = L^{-1} X$ and $Y^* = L^{-1} Y$ using forward substitution.
3. Compute OLS coefficients on transformed system:
   $$\hat{\beta} = (X^{*T} X^*)^{-1} X^{*T} Y^*$$
4. Compute residual variance $\hat{\sigma}^2$ and Log-Likelihood:
   $$\ln L = -\frac{N}{2} \ln(2\pi \hat{\sigma}^2) - \frac{1}{2} \ln |V| - \frac{N}{2}$$

## 4. Linear-Time $O(N)$ Algorithm (Ho & Ané, 2014)
For high-dimensional traits or large trees ($N > 1000$), `gelimo` incorporates the $O(N)$ tree-traversal algorithm described by **Ho & Ané (2014)** (*Systematic Biology*, [doi:10.1093/sysbio/syu005](https://doi.org/10.1093/sysbio/syu005)):

* **Key Advantage**: Computes Log-Likelihood and GLS estimates $\hat{\beta}$ directly via 3 post-order/pre-order tree passes without ever explicitly forming or inverting the $N \times N$ dense covariance matrix $V$.
* **Complexity**: Reduces computation time from $O(N^3)$ to $O(N)$ and memory overhead from $O(N^2)$ to $O(N)$.

---

## 5. Architectural Optimization Principles (Ponytail Strategy)

### 5.1. Parallel Tree & Parameter Optimization (SIMD + Rayon Zero-Copy)
In empirical phylogenetic studies, analyses evaluate thousands of posterior trees (to account for phylogenetic uncertainty) or fit multiple evolutionary models ($\lambda, \kappa, \delta, \text{OU}$) simultaneously.
* **R / Python Limitation**: Parallelizing duplicates tree structures across process boundaries (`multiprocessing`), causing severe memory inflation.
* **Rust / Rayon Solution**: Leveraging thread safety (`Send + Sync`) and zero-copy shared tree references, `gelimo` parallelizes likelihood evaluation across thousands of trees using native thread pools without copying data. Combined with SIMD auto-vectorization, CPU multi-core utilization reaches 100%.

### 5.2. Zero-Allocation In Optimization Loop
During Maximum Likelihood estimation of parameters like Pagel's $\lambda$, numerical optimizers evaluate dozens of candidate parameters.
* **R / caper Limitation**: Re-allocates matrices and vectors in R memory on every optimization step.
* **`gelimo` Solution**: Work buffers for node states, vector operations, and workspace structures are pre-allocated once before entering the optimization loop. Optimizers mutate these pre-allocated buffers in-place safely, achieving 0 heap allocations during parameter optimization loops.
