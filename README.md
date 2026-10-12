# CDPM2 UMAT

[日本語](README_JA.md)

Fortran implementation of CDPM2 for implicit analysis with Abaqus/Standard.
CDPM2 was developed by Grassl and co-authors; this project implements the existing constitutive model as a UMAT, based on the OOFEM reference implementation.

## Start here

- [English manual](docs/CDPM2_Manual_EN.pdf): theory, verification, inputs and usage.
- [Japanese manual](docs/CDPM2_Manual_JA.pdf).
- [UMAT source](src/cdpm2.f90).
- [Examples](examples/).
- [Notched-beam comparison](examples/notched_beam/): six paired CDP/CDPM2 cases and result plotting.

No geometry-table preparation is required. The Fortran source reads initial mesh geometry from Abaqus-generated ASCII `.fil` output. The examples include the required output requests; see the manual when using your own model. Automatic restart is not yet supported.

Download the repository using **Code → Download ZIP**, or clone it. The manual is the single detailed guide. The repository contains the source, runnable examples, PDF manuals and license notices needed by users.

This is a research implementation, not an official SIMULIA or OOFEM product. Review the manual's limitations before use.

LGPL-2.1-or-later, without warranty. See [LICENSE](LICENSE) and [third-party notices](THIRD_PARTY_NOTICES.md). Third-party experimental material retains its owners' rights.
