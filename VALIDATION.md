# Validation and limitations — 2026-11-06

Tested source SHA-256:
`d373a8bee58490c07bccce8a60c16cf3c09b1986c1f484f1d84e6c187192fbd0`.
The packaged source is byte-identical; its older header modification-date line
is supplemented by the dated change notice below and in THIRD_PARTY_NOTICES.
No packaged binaries are needed to reproduce the source build.

## Material-point evidence

| Check | Result | Interpretation |
|---|---|---|
| Archived case2 failed attempts, four batches | 20+28+23+20 = 91 converged | Local replay, not full-model completion |
| Broader archived trials | 724/1094 converged, 370 rejected | Original r4 had 696 successes; none lost |
| Independent spectral divided differences | 160/160; max relative error 3.1148e-15 | Includes generalized slope at exact zero |
| Tangent width sweep | 700/724 meet three consecutive widths below 1e-4 | 15 widths, 1e-5 to 1e-12; 24 do not qualify |

On the same original 11-width grid, only 685 points qualify both before and
after the spectral change. Do not attribute the 685-to-700 increase solely to
the code change. At fixed perturbation 1e-9, maximum error increases to about
0.832 from 0.00851 in the previous candidate; perturbations can cross sign
boundaries. There is no uniform finite-difference error reduction claim.
Trial 947 has a protected multiple-root rejection under perturbation.
The private replay corpus is not included in this source distribution; these
counts are archived internal evidence, not an independently reproducible public
benchmark. The independent spectral test and small FE decks are included.

## Actual Abaqus/Standard 2026 HF0 tests

Each run used a separate directory, source/library hashes, an exclusive claim,
successful datacheck (DAT marker and no errors), then a single continue launch.
S, LE, SDV, U and RF were extracted from all saved frames.

| Comparison | Result |
|---|---|
| Small prescribed-traction C3D8, previous/new | 40 increments, 40 iterations, no cutbacks; saved values identical |
| Mixed-boundary axial reversal C3D8 | 41 increments, one cutback each; iterations 47 to 46; max stress difference 1.9001e-5 MPa |
| Rotated two-C3D8 shear reversal | 40 increments, 43 iterations, no cutbacks; saved values identical |
| Two C3D6 wedge shear tests, previous/new | 40 increments, 44 iterations, no cutbacks; saved values identical |
| Serial versus requested 2-thread shear | Saved values identical |
| Repeated requested 2-thread shear | Saved values identical |

The comparator is the preceding characteristic-length-tangent candidate,
not the old public release. SMP results are a small smoke test, not a race-free
proof or MPI validation. These historical verification materials use As=15;
they are not the case2 calibration (As=5).

Machine-readable comparison summaries are in `evidence/`. Restart results are
reported separately there. The initial restart harness omitted the old ODB and
failed datacheck; that failed attempt is retained locally, not counted as a pass.

The corrected 89-SDV automatic-FIL run and restart both completed. Restart used
copied restart/ODB files without the original FIL or external geometry table.
Across 22 common frames, all 27 geometry entries match exactly; maximum S and
material SDV differences are 3.58e-13 and 3.56e-13. U differs by at most 2.63e-18.
RF differs by 5.666e-9, exceeding the initial absolute 1e-9 comparison threshold.
That strict all-field check is therefore **not passed**; the threshold was not
silently relaxed. These are observed small differences, not exact equivalence.

Follow-up frame-level RF audit: the 5.666e-9 maximum occurs only in the
restart-initial frame (phase3, time zero, instance A node 3, component 3).
All 21 later saved frames, including the next step's initial frame, have RF
differences at most 4.247e-15. Thus the observed discrepancy does not persist
in the subsequent computed response. Casting both values to float32 does not
make all initial-frame values equal; simple output rounding is not established
as the cause. Restart initialization/output reconstruction is a hypothesis,
not a proven implementation diagnosis. The original all-frame failed criterion
remains recorded. See `evidence/RESTART_RF_AUDIT.json`.

## Not established

Production qualification, arbitrary loading and mesh objectivity, broad global
Newton robustness, general parallel safety, completion of stopped case2, and
experimental validation are not established. Use an explicitly marked
research prerelease; do not silently replace stable-release validation claims.

## Full 312-condition old/new regression completed

The final adopted conditions from the old public `latest.csv` were used: 284
at 240 increments, 26 at 960, and 2 at 3840. Eighteen original batch inputs
per source were run with one CPU (36 analyses). All completed, with zero MSG
errors and successful extraction. All 312 selected conditions had complete,
identical time grids; all eight integration points were compared. Additional
specimens in shared batches were not counted as extra distinct conditions.

Old source SHA-256: `78cc9a681a3526cd1d04e7232ada14ecb18cb116ecd4426b186029ef5536b2de`.
New source SHA-256 is given above. Both rebuilt libraries matched their earlier
recorded hashes. Materials, viscosity and original increment grids were retained.

Maximum absolute differences across selected conditions:

| Output | Maximum difference |
|---|---:|
| Stress S | 3.81149e-9 MPa |
| LE | 2.97912e-57 |
| SDV14 | 8.49859e-8 |
| SDV37 / SDV38 | 1.92469e-10 / 7.59535e-11 |
| SDV39 (characteristic length) | 0 |
| SDV59 / SDV60 | 1.92423e-10 / 7.58806e-11 |

Maximum stress difference scaled by max(1 MPa, old history maximum) is
3.37009e-11. Solver wall-time sums were 164.40 s old and 158.93 s new; this
single serial campaign is not a statistically controlled speed benchmark.
The differences show no appreciable response regression in these conditions.
No post-hoc numerical PASS threshold was assigned: machine-readable status is
COMPARED for 312 conditions, separate from solver completion.

This is an old/new UMAT comparison, not a fresh native OOFEM run. The historical
five OOFEM discrepancies remain distinct and are not claimed resolved. The
prescribed strain paths do not establish global Newton convergence on general
structural models. See `evidence/regression312_summary.json` and the full
per-condition `evidence/regression312_comparison.json`.

## Implementation changes (2026-11-06)

Normalized/log-domain plastic flow; scaled plastic multiplier; balanced linear
solve fallback; root-checked algebraic continuation; characteristic-length
direction sensitivity; spectral divided-difference tangent. Material constants,
viscosity and acceptance tolerances were not relaxed. Persistent material-state
layout was not changed. Numerical algorithms do change and require project-
specific verification before adopting this preview.

Documentation basis: Abaqus 2026, User Subroutines > Abaqus/Standard User
Subroutines > UMAT (`SIMACAESUBRefMap/simasub-c-umat.htm`, initial single-element
traction testing); Execution > File Extension Definitions
(`SIMACAEEXCRefMap/simaexc-c-fileextensions.htm`, restart files). The manuals
are not distributed and do not certify this implementation.
