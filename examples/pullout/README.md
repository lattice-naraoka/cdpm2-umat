# Reinforcing-bar pullout

[日本語](README_JA.md)

A quarter-model Abaqus/Standard example comparing the CDPM2 UMAT and built-in CDP. Both cases use the same mesh, steel, support, contact and loading conditions. CDP tensile and compressive damage are capped at 0.99. Viscosity is 1e-4 for CDP and 1e-3 for CDPM2.

- [Example guide and results](Pullout_Example_EN.pdf)
- [Japanese guide](Pullout_Example_JA.pdf)

## Run

Use licensed Abaqus and a supported Fortran compiler. From this folder, run:

```text
abaqus job=cdp interactive
abaqus job=cdpm2 user=../../src/cdpm2.f90 interactive
```

The CDPM2 source is shared with the repository's other examples. Keep the supplied ASCII FIL requests and `*DEPVAR` value 89. The mesh is in `mesh.inp`; CDP material tables are in `material.inp`. The [main manual](../../docs/CDPM2_Manual_EN.pdf) describes the UMAT installation and material model.

## Results

Both cases completed to a loaded-end displacement of 2 mm using Abaqus/Standard 2026 HF0. Full-specimen peak loads were 178.85 kN (CDP) and 139.95 kN (CDPM2). The guide includes the load-slip comparison and maximum principal logarithmic strain contours.

Extract the completed ODB files, then plot with Python, NumPy and Matplotlib:

```text
abaqus python extract_curve.py cdp
abaqus python extract_curve.py cdpm2
python plot_results.py
```

Extraction writes `cdp_curve.csv` and `cdpm2_curve.csv`. `comparison.png` plots the newly extracted analysis curves. Load is four times the quarter-model loaded-end reaction. Slip is the mean unloaded-end bar displacement relative to the mean concrete displacement at (50,0,0) and (0,50,0) mm.

## Model and sources

Units: N, mm, MPa. The model contains 86,152 C3D8R elements, a 200 mm concrete cube represented by quarter symmetry, and a nominal 20 mm bar with 100 mm bond length. Ribs are rings; the entire concrete top face is restrained vertically. Linear penalty contact uses a 0.005 stiffness multiplier and friction coefficient 0.5. The guide gives the material properties, geometry and sources.

This example was constructed with reference to S. Seok, *Finite Element Modeling of Bond-Zone Behavior in Reinforced Concrete*, Purdue University, 2019, Chapter 4 ([record](https://docs.lib.purdue.edu/dissertations/AAI30502919/)), and G. Metelli and G. A. Plizzari, *Influence of the relative rib area on bond behaviour*, 2014 ([DOI](https://doi.org/10.1680/macr.13.00198)). Ring ribs and the adopted support conditions are modelling choices, not the original test fixture. Experimental curves in the guide are approximate readings of Seok Fig.4.8; original source figures and measured datasets are not included.
