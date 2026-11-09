# CDPM2 UMAT

[日本語](README_JA.md)

Fortran implementation of the existing CDPM2 model of Grassl and coauthors for implicit Abaqus/Standard, based on the OOFEM reference implementation.
**Research preview: not production-qualified and not an official SIMULIA or OOFEM product.**

## Documentation

- [Manual](docs/CDPM2_Manual_EN.pdf): current material model, numerical update, inputs and usage.
- [Verification guide](docs/CDPM2_Verification_EN.pdf): tested source, methods, results and limitations.
- [Analysis examples](examples/pullout/CDPM2_Examples_EN.pdf): contents, Example 1 and Example 2.
- [Change history (bilingual)](docs/CDPM2_Change_History.pdf): one consolidated history PDF.

Documentation edition: 2026-11-09. Examples contents are on page 1, Example 1 on pages 2-6, and Example 2 on pages 7-11.

## Source and tested scope

- [UMAT source](src/cdpm2.f90)
- [Validation details](VALIDATION.md) and [machine-readable evidence](evidence/)
- [Example 1 inputs](examples/pullout/), [Example 2 inputs](examples/pullout/case2/) and [Example 2 notes](CASE2.md)
- [Notched-beam inputs](examples/notched_beam/) and [one-element inputs/reference data](verification/)

Regression against the reference source covers 312 conditions: all 36 runs completed, with maximum stress difference 3.81149e-9 MPa.
This is not a fresh native OOFEM campaign or experimental validation. Archived local trials give 724 converged and 370 rejected out of 1094;
the specified 15-width tangent criterion qualifies 700/724 points. Same-source restart exceeds its initial-frame reaction criterion,
so all-field equivalence is not established. See the verification guide for scope and unresolved issues.
Existing saved examples and reference datasets are not reruns of the currently distributed source.

## Usage

Use licensed Abaqus and a supported Fortran compiler, following the manual. Initial geometry is read from ASCII FIL output;
basic input uses 24 properties and 89 state entries. Cross-version restart compatibility is not guaranteed.
Download using **Code → Download ZIP** or git clone.

Independent spectral test: `python tests/spectral_reference.py --compiler gfortran`
Never link the standalone test utility stubs into Abaqus.

LGPL-2.1-or-later, without warranty. Retain [LICENSE](LICENSE) and [third-party notices](THIRD_PARTY_NOTICES.md).
