# Stable Inversions

In DQMC, to obtain the equal-times Green's function we commonly perform the inversion

$$G = \left[1 + B\right]^{-1}.$$

Similarily, we compute

$$G = \left[A + B\right]^{-1}$$

for the time-displaced Green's function. The following methods are exported to facilitate these tasks.

- `inv_one_plus`, `inv_one_plus!`
- `inv_sum`, `inv_sum!`

When function names are suffixed with `_loh`, i.e. `inv_one_plus_loh`, a more sophisticated method is used for numerical stabilization (see the paper linked above for more details).
