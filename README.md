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

No geometry-table preparation is required. The Fortran source reads initial mesh geometry from Abaqus-generated ASCII `.fil` output. Use `*DEPVAR` with **89** variables: the original 62 material/work variables are unchanged; the additional 27 retain reference geometry in Abaqus restart data. The examples include the required settings. Restart jobs do not need the original FIL or a separate geometry file. Existing 62-variable inputs still run, but their old restart files do not contain the added geometry. See the manual for details.

Download the repository using **Code → Download ZIP**, or clone it. The manual is the single detailed guide. The repository contains the source, runnable examples, PDF manuals and license notices needed by users.

This is a research implementation, not an official SIMULIA or OOFEM product. Review the manual's limitations before use.

LGPL-2.1-or-later, without warranty. See [LICENSE](LICENSE) and [third-party notices](THIRD_PARTY_NOTICES.md). Third-party experimental material retains its owners' rights.
