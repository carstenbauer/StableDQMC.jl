# Stable Inversions

In DQMC, we commonly perform inversions like `G = [1 + B]^-1` to obtain the equal-times Green's function and `G = [A + B]^-1` for the time-displaced pendant. The following methods are exported to facilitate these tasks.

- `inv_one_plus`, `inv_one_plus!`
- `inv_sum`, `inv_sum!`

When function names are suffixed with `_loh`, i.e. `inv_one_plus_loh`, a more sophisticated method is used for numerical stabilization (see the paper linked above for more details).
