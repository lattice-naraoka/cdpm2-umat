# Notched beam: CDP and CDPM2

[日本語](README_JA.md)

Six paired Abaqus/Standard examples: 70 (`coarse`), 280 (`medium`) and 1120 (`fine`) C3D8R elements, Riks loading, no viscosity. Each pair has the same mesh, boundary conditions and step controls. CDPM2 replaces the material, adds its required hourglass stiffness and changes material output to SDV. Its ASCII FIL output requests supply initial mesh geometry to the UMAT.

## Run

Build `../../src/cdpm2.f90` in the repository's `build/` folder as described in the [manual](../../docs/CDPM2_Manual_EN.pdf). From this folder, run one job at a time:

```text
abaqus job=cdp_coarse cpus=1 mp_mode=threads interactive output_precision=full
abaqus job=cdpm2_coarse cpus=1 mp_mode=threads interactive output_precision=full
```

Replace `coarse` with `medium` or `fine` for the other meshes. No geometry table or mesh-preparation tool is needed. Keep the included ASCII FIL and mesh-wide file-output requests when adapting a CDPM2 input. Do not run the same job name twice concurrently. Automatic restart is not supported, although restart output is written by the retained step settings.

Extract each completed ODB with Abaqus Python, then plot with ordinary Python (NumPy and Matplotlib):

```text
abaqus python extract_curve.py cdp_coarse.odb
abaqus python extract_curve.py cdpm2_coarse.odb
python plot_comparison.py
```

Repeat extraction for all six jobs before plotting. The plot uses your newly extracted CDP/CDPM2 curves and is saved to `results/comparison.png`. The manual includes the existing verification results and experimental comparison.

## Conditions and source

These example inputs were reconstructed with reference to Abaqus 2026 Benchmarks, *Notched unreinforced concrete beam under 3-point bending*, official `notchedconcbeam_3d_{coarse,medium,fine}_std.inp`. CDP tensile softening and tensile damage are replaced with the official 2D TS-2 tables. This is not a reproduction of the original 3D linear-softening case. Mesh-generator keywords were expanded to explicit nodes, elements and sets, without changing physical geometry or mesh connectivity. The UMAT reads Abaqus-generated FIL geometry, not these input keywords.

Units: N, m, Pa. E=30 GPa, nu=0.2, fc=30 MPa, ft=3.33 MPa. TS-2 opening/stress points: (0 mm, 3.33 MPa), (0.025 mm, 0.833 MPa), (0.198 mm, 0.0333 MPa). The finite table area is 0.12697 N/mm; CDPM2 has the same bilinear-softening area, but is not identical to CDP's 1% residual tail. Model widths are 20/10/5 mm. Curves are converted to full-beam load at 50 mm width. Riks step time is not physical loading time. The manual contains the remaining material details and limitations.

Experiment: Petersson (1981), *Crack Growth and Development of Fracture Zones in Plain Concrete and Similar Materials*, Report TVBM-1006, Lund Institute of Technology. Direct graph-reading source: the above Abaqus benchmark, Fig.11. Points are approximate readings, not raw experimental records. Gf=115/137 N/m identifies experimental series, not the present analysis input. Third-party material retains its original rights; the project license does not relicense it.
