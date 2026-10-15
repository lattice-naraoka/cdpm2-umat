# Norm-definition isolation tests

These five comparisons use a **diagnostic copy**, not the public UMAT. Their
results do not change the original PASS/FAIL records or the 342-run total.
All tests retain viscosity `1e-3`, original material inputs and the same criteria.

Only the plastic-increment norm in `update_damage` and its derivative in
`damage_sensitivity` differ from the public source:

```fortran
! Public tensor norm:
pnorm = norm2([dep(1:3), dep(4:6)/sqrt(2.0_real64)])
! Diagnostic OOFEM-compatible vector norm:
pnorm = norm2(dep)
```

The derivative uses the corresponding unweighted six-component dot product
instead of a half-weight on engineering-shear components. No yield function,
hardening parameters, softening parameters or convergence criteria are changed.
This is an explanation of the experiment, not a recommendation to alter the
public UMAT. The diagnostic source and library hashes are recorded separately.

`norm_isolation.csv` records the five diagnostic comparisons. The other JSON
files record the unchanged OOFEM coordinate-representation test and the unchanged
public UMAT comparison against the principal-axis OOFEM reference. These controls
use `Hp=0.5`, 960 increments, period 1000 and common characteristic length sqrt(2).
They test this specific path and do not establish arbitrary-rotation objectivity.
