# Stable Inversions

## Overview

In DQMC, to obtain the equal-times Green's function we commonly perform the inversion

$$G = \left[1 + B\right]^{-1}.$$

Similarily, we compute

$$G = \left[A + B\right]^{-1}$$

for the time-displaced Green's function. The following methods are exported to facilitate these tasks.

- `inv_one_plus`, `inv_one_plus!`
- `inv_sum`, `inv_sum!`

When function names are suffixed with `_loh`, i.e. `inv_one_plus_loh`, a more sophisticated method is used for numerical stabilization (see the paper linked above for more details).

## Details

### `inv_one_plus`

```math
\begin{aligned}
	G &= [\mathbb{1} + UDX]^{-1} \\
	&= [U\underbrace{(U^\dagger X^{-1} + D)}_{udx}X]^{-1}\\
	&= [(Uu)d(xX)]^{-1}\\
	&= U_r D_r X_r,
\end{aligned}
```
with $U_r = (xX)^{-1}$, $D_r = d^{-1}$, $X_r = (Uu)^{-1}$.

### `inv_one_plus_loh`

```math
\begin{aligned}
	G &= [\mathbb{1} + UDX]^{-1}\\
	&= [\mathbb{1} + UD_mD_pX]^{-1}\\
	&= [(X^{-1} D_p^{-1} + U D_m) D_p X]^{-1}\\
	&= X^{-1} \underbrace{[D_p^{-1} (\underbrace{X^{-1} D_p^{-1} + UD_m}_{udx})^{-1}]}_{udx} \\
	&= U_r D_r X_r,
\end{aligned}
```
with $D_m = \min(D, 1)$, $D_p = \max(D, 1)$, $U_r = X^{-1}u$, $D_r = d$, and $X_r = x$. [Loh2005, Loh1989]


### `inv_sum`

```math
\begin{aligned}
	G(\tau_1, \tau_2) &= [U_L D_L X_L + U_R D_R X_R]^{-1}\\
	&= [U_L \underbrace{(D_L X_L X_R^{-1} + U_L^\dagger U_R D_R)}_{udx} X_R ]^{-1}\\
	&= [(U_L u) d^{-1} (x X_R)]^{-1}\\
	&= U_r D_r X_r,
\end{aligned}
```
where $U_r = (x X_R)^{-1}$, $D_r = d^{-1}$, and $X_r = (U_L u)^{-1}$.


### `inv_sum_loh`

```math
\begin{aligned}
	G(\tau_1, \tau_2) &= [U_L D_L X_L + U_R D_R X_R]^{-1} \\
	&= [U_L D_{Lm} D_{Lp} X_L + U_R D_{Rm} D_{Rp} X_R]^{-1} \\
	&= \left[U_L D_{Lp} \underbrace{\left( \dfrac{D_{Lm}}{D_{Rp}} X_L X_R^{-1} + U_L^\dagger U_R \dfrac{D_{Rm}}{D_{Lp}} \right)}_{udx} X_R D_{Rp} \right]^{-1} \\
	&= X_R^{-1} \underbrace{\dfrac{1}{D_{Rp}} [udx]^{-1} \dfrac{1}{D_{Lp}}}_{udx} U_L^\dagger\\
	&= U_r D_r X_r,
\end{aligned}
```
with $D_{Rm} = \min(D_R, 1)$, $D_{Rp} = \max(D_R, 1)$, $D_{Lm} = \min(D_L, 1)$, $D_{Lp} = \max(D_L, 1)$, $U_r = X_R^{-1} u$, $D_r = d$, and $X_r = x U_L^\dagger$. [Loh2005]


## Resources

#### [[Loh2005]](http://www.worldscientific.com/doi/abs/10.1142/S0129183105007911)
```
@article{Loh2005,
	author = {Loh, E. Y. and Gubernatis, J. E. and Scalettar, R. T. and White, S. R. and Scalapino, D. J. and Sugar, R. L.},
	title = {{Numerical Stability and the Sign Problem in the Determinant Quantum Monte Carlo Method}},
	journal = {International Journal of Modern Physics C},
	volume = {16},
	number = {08},
	pages = {1319--1327},
	month = {aug},
	year = {2005}
	issn = {0129-1831},
	doi = {10.1142/S0129183105007911},
	url = {http://www.worldscientific.com/doi/abs/10.1142/S0129183105007911},
}
```

#### [[Loh1989]](http://link.springer.com/10.1007/978-1-4613-0565-1{\_}8)
```
@incollection{Loh1989,
	author = {Loh, E. Y. and Gubernatis, J. E. and Scalettar, R. T. and Sugar, R. L. and White, S. R.},
	title = {{Stable Matrix-Multiplication Algorithms for Low-Temperature Numerical Simulations of Fermions}},
	year = {1989}
	pages = {55--60},
	doi = {10.1007/978-1-4613-0565-1_8},
	url = {http://link.springer.com/10.1007/978-1-4613-0565-1{\_}8},
}
```
