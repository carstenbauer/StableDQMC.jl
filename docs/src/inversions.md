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
```math
with $D_m = \min(D, 1)$, $D_p = \max(D, 1)$, $U_r = X^{-1}u$, $D_r = d$, and $X_r = x$.


### `inv_sum`

```math
\begin{aligned}
	G(\tau_1, \tau_2) &= [U_L D_L X_L + U_R D_R X_R]^{-1}\\
	&= [U_L \underbrace{(D_L X_L X_R^{-1} + U_L^\dagger U_R D_R)}_{udx} X_R ]^{-1}\\
	&= [(U_L u) d^{-1} (x X_R)]^{-1}\\
	&= U_r D_r X_r,
\end{aligned}
```math
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
with $D_{Rm} = \min(D_R, 1)$, $D_{Rp} = \max(D_R, 1)$, $D_{Lm} = \min(D_L, 1)$, $D_{Lp} = \max(D_L, 1)$, $U_r = X_R^{-1} u$, $D_r = d$, and $X_r = x U_L^\dagger$.
