# One-element verification

[日本語PDF](../docs/CDPM2_Verification_JA.pdf) | [English PDF](../docs/CDPM2_Verification_EN.pdf)

The bilingual Verification Guide reports native OOFEM / Abaqus UMAT comparisons
for 24 prescribed strain paths and 13 material sets. Every UMAT run uses viscosity
`1e-3`; OOFEM is rate-independent. Analytical elasticity provides an additional
independent reference. Completion and numerical acceptance are reported separately.

The catalogue includes monotonic, reversed, cyclic, held, multiaxial and
fixed-principal-axis nonproportional histories. Side strains are prescribed;
these are not traction-free uniaxial-stress tests. Each specimen is one C3D8
element, with all eight integration points assessed. Disconnected specimens with
identical properties share an Abaqus batch only to reduce launch overhead.

## Reproduce

Install NumPy, native OOFEM and Abaqus/Standard with its supported Fortran compiler.
Keep the repository's `src/cdpm2.f90` unchanged.

```text
python generate.py
python run.py suite --oofem /path/to/oofem --reference-source /path/to/concretedpm2.C --abaqus abaqus
python analyze.py suite
```

The generator uses `NLGEOM=YES` and exponential nodal stretches so that Abaqus
logarithmic strains match the native OOFEM small strains at the comparison times.
Both use the same loading duration and increment grid. OOFEM receives the
independently derived geometric characteristic length of each cubic specimen:
`L` for axis-aligned paths and `sqrt(2)*L` for the 45-degree shear paths.
The UMAT computes this length from its initial geometry.

Generate the recorded refinement tests separately:

```text
python generate.py --group Hp_05 --n 960 --out refine_Hp_960
python generate.py --group E_high --case compression_cycles --n 960 --out refine_Ehigh_960
python generate.py --group As_5 --group Gf_half --group length_2 --case shear_reversal --n 960 --out refine_shear_960
python generate.py --group Hp_05 --case pure_shear_xy --case shear_reversal --n 3840 --out refine_Hp_3840
```

Run and analyze each generated directory as above. The runner uses one CPU and
does not overwrite existing results. Site-specific license scheduling is external
to the reproducible tests. The OOFEM executable and source hashes are recorded.
The reference executable used for the reported runs identifies itself as OOFEM
2.6; its embedded Git identifier is unavailable, so its exact executable hash is
provided rather than an inferred binary commit identifier.
The reference material source is pinned to OOFEM commit
[`1c1dbd2`](https://github.com/oofem/oofem/tree/1c1dbd2b4b8c95dcc26878bf5009d44ed535e1a6).
The recorded CRLF source matches the official LF source after line-ending
normalization; both source hashes are supplied in `provenance.json`.

## Results

The guides and per-run CSV files include unsuccessful and incomplete comparisons,
including the original coarse-increment run. Refinement runs do not erase those
results. A batch interruption prevents complete assessment of its other specimens;
it does not establish a material failure in every specimen of that batch.

Five original OOFEM comparisons exceed the criteria because the two implementations
use different shear weights in the plastic-strain-increment norm entering damage
history. The public UMAT uses the tensor Frobenius norm; the fixed OOFEM version
uses an unweighted Euclidean norm of engineering-Voigt components in these routines.
The guide documents controlled norm-only isolation tests and a coordinate-
representation check. The five diagnostic-copy comparisons meet the original
criteria, but are reported separately: they do not change `results.csv`,
`latest.csv`, the 342-run count or the public UMAT source.

`diagnostics/` contains the numerical isolation results and reference-identity
checks. Its results are not release-UMAT acceptance results. The finding concerns
the fixed reference version and the tested paths, not all OOFEM versions.

These tests assess the explicitly defined constitutive paths. They do not establish
global Newton convergence, arbitrary-rotation objectivity, general mesh objectivity,
experimental validation or application-specific qualification.
