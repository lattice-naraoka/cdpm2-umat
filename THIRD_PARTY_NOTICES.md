# Copyright and modification notices

[日本語](THIRD_PARTY_NOTICES_JA.md)

Source, added tools, tests and project-authored explanations are GNU Lesser General Public License version 2.1 or later (LGPL-2.1-or-later). Preserve the complete LICENSE and original notices. No warranty. Source-only distribution does not waive obligations for binary redistribution. LICENSE contains the authoritative legal text; these explanations are not substitute license terms.

## OOFEM

Copyright (C) 1993 - 2013 Borek Patzak

Reference commit: 1c1dbd2b4b8c95dcc26878bf5009d44ed535e1a6. The original constitutive sources (concretedpm2.C, concretedpm2.h), supporting sources (floatmatrixf.h, structuralmaterial.C) and license (COPYING.LGPLv2.1) are available from the upstream reference below.

Reference: https://github.com/oofem/oofem/tree/1c1dbd2b4b8c95dcc26878bf5009d44ed535e1a6

Constitutive source: Grassl et al., CDPM2: A damage-plasticity approach to modelling the failure of concrete, International Journal of Solids and Structures 50(24), 2013, 3805-3816. No substitution with a different damage model is claimed.

## Modifications

Source modifications were made from 2026-09-30 to 2026-10-12. The OOFEM reference model was ported to Fortran with the Abaqus stress, engineering-shear, STATEV, rotation and mesh-geometry interface; direct endpoint update; hybrid tangent; viscous relaxation; and failure-state protection/external cutback. Changes include an additional Newton initial guess, comparison of converged solutions and acceptance control on damage increments. Mesh geometry is now read automatically from Abaqus ASCII .fil files without a prepared geometry table. C3D8I eight-node geometry is also supported. Added tests use the same LGPL conditions. The automatic input adapter does not change the constitutive update.

OOFEM, Abaqus and HYPLAS are separate projects. This is not an OOFEM- or SIMULIA-provided/endorsed UMAT. The nonlinear-FEM book and HYPLAS informed general numerical methods; their pages/programs are not included. Abaqus software, manuals and libraries are not included. Third-party experimental material retains its owners' rights. Comparison plots are independently generated from simulation results and numerical values read from published plots; original third-party images and measured datasets are not distributed.
