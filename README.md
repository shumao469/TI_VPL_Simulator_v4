# Mouse TI Simulator

[![Open in Streamlit](https://static.streamlit.io/badges/streamlit_badge_black_white.svg)](https://tivplsimulatorv4-dp7r8p99mgktdkfyjkmhyf.streamlit.app/)
[![CI](https://github.com/shumao469/TI_VPL_Simulator_v4/actions/workflows/ci.yml/badge.svg?branch=main)](https://github.com/shumao469/TI_VPL_Simulator_v4/actions/workflows/ci.yml)

A Xu Lab research prototype for exploring temporal interference (TI) electric fields around the mouse ventral posterolateral thalamic nucleus (VPL). It combines **Allen Mouse CCFv3 2017 anatomical annotations**, a two-channel, four-electrode montage, and TI envelope calculations. The web application provides a fast analytical preview; a separate local workflow supports finite element modeling (FEM) using an EFMouse mesh.

**[Open the online simulator](https://tivplsimulatorv4-dp7r8p99mgktdkfyjkmhyf.streamlit.app/)** — no Python installation required. If the application is asleep, follow the hosting page's wake-up instructions and wait for it to load. Availability depends on Streamlit hosting.

[Chinese documentation](README_CN.md)

<img width="1424" alt="Mouse TI Simulator interface example" src="https://github.com/user-attachments/assets/09c81e7b-06b4-4ce6-a3a8-5f96415368ac" />

## 1. Which Allen Institute atlas is used?

The web asset records its provenance in [assets/allen_ccf_100um_web.json](assets/allen_ccf_100um_web.json). Its generation pipeline is implemented in [scripts/build_atlas_assets.py](scripts/build_atlas_assets.py).

| Item | Current implementation |
|---|---|
| Atlas | Allen Mouse Brain Common Coordinate Framework, version 3 (**CCFv3**) |
| Annotation release | **2017 annotation (`ccf_2017`)** |
| Resolution used by the web application | **100 µm isotropic** |
| Annotation volume | [annotation_100.nrrd](https://download.alleninstitute.org/informatics-archive/current-release/mouse_ccf/annotation/ccf_2017/annotation_100.nrrd) |
| Corresponding average template | [average_template_100.nrrd](https://download.alleninstitute.org/informatics-archive/current-release/mouse_ccf/average_template/average_template_100.nrrd) |
| Target structure | VPL, Allen structure ID **718**; the web application extracts the left hemisphere |
| Anatomical uses | Whole-brain display mask, transparent brain surface, left VPL surface, and slice outlines |

CCFv3 is a three-dimensional reference space whose average template was constructed from 1,675 specimen brains. See the [Allen documentation](https://brain-map.org/support/documentation) and [Wang et al., 2020](https://doi.org/10.1016/j.cell.2020.04.007). The atlas supplies anatomical references; electrical conductivities and electric fields are not measurements supplied by the atlas.

The configuration retains an `atlas_name: allen_mouse_25um` field. However, the current web assets and advanced script actually read **100 µm** annotations. That field does not establish that the current calculation uses 25 µm data.

### Coordinate conventions

Interface coordinates are in millimeters relative to bregma: **AP positive anteriorly, ML positive to the right, and DV positive dorsally**. The default requested target is `AP=-1.94, ML=-1.75, DV=-3.10 mm`.

The web application uses a generic QC mapping: `AP=5.4−CCF_AP`, `ML=CCF_ML−5.7`, and `DV=−CCF_DV`, with all distances in mm and CCF array axes ordered AP/DV/ML. This approximation does not replace subject-specific or validated template registration.

The asset records an approximately **0.64 mm** distance from the default target to the nearest left VPL voxel. The sidebar reads this precomputed value; **it is not recalculated when the target inputs change**.

## 2. Simulation methods

| Mode | Method | Intended use |
|---|---|---|
| Online or local web preview | Analytical point source–sink fields in a homogeneous conductor, displayed with Allen anatomy, followed by TI envelope calculation | Qualitative exploration of electrode placement, coordinates, current balance, and parameters |
| Local advanced FEM | Two-channel electroquasistatic FEM on an EFMouse multi-tissue tetrahedral mesh | Template electric fields and numerical QC under specified material, mesh, and registration assumptions |

### Analytical preview

Implemented in [ti_vpl/analytic.py](ti_vpl/analytic.py) and [app_advanced.py](app_advanced.py). Each source–sink electrode pair produces:

```text
E(r) = k I / (4πσ) · [(r−rs)/|r−rs|³ − (r−rt)/|r−rt|³]
```

The web application uses `σ=0.275 S/m` and a fixed scaling factor `k=0.025`. The distance denominator has a **0.45 mm** lower bound to limit near-electrode singularities. The underlying `PreviewParameters` default factor is 0.028; the web application explicitly overrides it with 0.025. The repository does not provide experimental calibration of this factor, so absolute preview amplitudes should not be used to determine animal stimulation doses.

**Allen masks constrain the display region; they do not enter a multi-tissue conductivity boundary-value calculation.** Atlas-derived brain anatomy does not mean that the web heatmap has undergone realistic head-model FEM.

### Finite element modeling

The advanced entry point is [scripts/run_advanced_fem.py](scripts/run_advanced_fem.py), with the Python solver in [ti_vpl/fem.py](ti_vpl/fem.py). It converts the `Alekseichuk_mesh_clean.mat` mesh from the [official EFMouse repository](https://github.com/klabhub/EFMouse) and applies linear tetrahedral FEM. It does not directly call the EFMouse MATLAB toolbox, and equivalence to that toolbox's results has not been established.

```text
∇·(σ∇φ) = 0
E = −∇φ
```

The E1–E2 and E3–E4 channels are solved separately. Electrode contact surfaces receive balanced Neumann inflow and outflow currents; the remaining boundary has no imposed flux, and one node's potential is fixed. Small systems use a sparse direct solver. Larger systems use conjugate gradients, preferably with an algebraic multigrid (AMG) preconditioner, falling back to diagonal preconditioning when AMG is unavailable.

Configured tissue conductivities, in S/m, are: gray matter 0.275, cerebrospinal fluid 1.654, bone 0.010, skin 0.465, and eye 0.500. These are model assumptions that require assessment for the intended experiment. The advanced entry point uses real, isotropic, frequency-independent conductivities. It does not simulate neuronal responses, heating, or frequency dispersion.

### TI envelopes and frequency

[ti_vpl/envelope.py](ti_vpl/envelope.py) calculates the following for two real vector fields:

```text
Directional envelope: m(n) = 2 min(|n·E1|, |n·E2|), |n|=1
Maximum-orientation envelope: mTI = max_n m(n)
Beat frequency: Δf = |f2−f1|
```

The maximum-orientation envelope selects a direction independently at every location. It is an upper-bound metric, not the stimulation strength along an actual nerve fiber or a probability of neural activation.

Defaults are `f1=2000 Hz` and `f2=2010 Hz`, giving a 10 Hz beat frequency. In the current model, changing a nonzero frequency difference changes the temporal beat rate but not the spatial amplitude. For equal frequencies, the web/analytical preview sets the TI beat envelope to zero. The advanced FEM entry point currently lacks this equal-frequency branch: use distinct carriers for TI analysis, and do not interpret its equal-frequency spatial formula output as a beat envelope.

## 3. Using the online application

1. Open the online link and wait for the interface and brain model to load.
2. In the sidebar, enter target AP/ML/DV coordinates and adjust both channel currents, f1/f2, cross half-spacing, and montage rotation.
3. Click the parameter submission button, labeled **“运行并提交新参数”** (“Run and submit new parameters”). Check the Run number and electrode coordinate table. Editing inputs alone does not submit a new Run.
4. Drag the 3D view to rotate and scroll to zoom. Inspect the coronal, sagittal, and horizontal slices alongside the 3D field.
5. Review `target mTI` (the analytical envelope at the target point), `pair current balance` (smaller channel current divided by larger channel current), and `beat frequency`. The current balance percentage is not a numerical current-conservation residual.

**The red/pink object is the left VPL anatomical label, not an electric-field focus.** Colored nodes and slice heatmaps show the calculated field. The color scale is fixed at 0–0.35 V/m; higher values may saturate, so the color-scale limit does not establish the maximum field.

Defaults are 0.50 mA per channel, a 2.08 mm half-spacing, 0° rotation, and electrode DV=+0.50 mm:

| Electrode | AP (mm) | ML (mm) | DV (mm) | Channel |
|---|---:|---:|---:|---|
| E1 | +0.14 | −1.75 | +0.50 | E1–E2 |
| E2 | −4.02 | −1.75 | +0.50 | E1–E2 |
| E3 | −1.94 | +0.33 | +0.50 | E3–E4 |
| E4 | −1.94 | −3.83 | +0.50 | E3–E4 |

E3 crosses the midline. This montage is an adjustable starting geometry. Equal currents and geometric symmetry do not guarantee that the VPL contains the whole-brain field maximum.

## 4. Downloading and launching locally

### Web preview only: Windows, Linux, or macOS

Install Git and **Python 3.11**, then run:

```bash
git clone https://github.com/shumao469/TI_VPL_Simulator_v4.git
cd TI_VPL_Simulator_v4
python -m venv .venv
```

Activate the environment:

- Windows PowerShell: `.\.venv\Scripts\Activate.ps1`
- Linux/macOS/WSL: `source .venv/bin/activate`

Then install the preview dependencies and launch:

```bash
python -m pip install --upgrade pip
python -m pip install -r requirements.txt
python -m streamlit run streamlit_app.py --server.port 8504
```

Open **[http://localhost:8504](http://localhost:8504)** in a browser. Keep the terminal running; press Ctrl+C to stop. In an activated environment, you can also launch `run_web.bat` on Windows or `bash run_web.sh` on Linux/WSL. If the port is occupied, use 8505 and open the corresponding address.

The preview uses the bundled anatomical assets and requires neither AllenSDK nor the large EFMouse mesh.

### Full FEM environment: Linux/WSL with Conda or Miniforge recommended

From the repository root:

```bash
bash create_environment.sh
conda activate ti_vpl_fem311
bash run_web.sh
```

The script creates or updates `ti_vpl_fem311` using `environment.yml` and installs the project. If gmsh/VTK reports missing OpenGL libraries on Ubuntu/WSL, install `libglu1-mesa` and `libgl1`. The full environment is larger than the preview environment; setup time depends on downloads and dependency resolution.

## 5. Running local FEM

Run the following from the repository root with the full environment activated.

**Step 1: Download and convert the mesh.** The large MAT/VTU files are not bundled. The download script records a SHA-256 digest.

```bash
python scripts/download_efmouse_assets.py
python scripts/convert_efmouse_mesh.py
```

**Step 2: Prepare real landmarks and register the model.** `registration/landmarks_example.csv` demonstrates the format only; it is not an experimental registration. Create `registration/landmarks.csv` with these columns:

```text
bregma_AP,bregma_ML,bregma_DV,fem_x,fem_y,fem_z
```

Both coordinate sets use mm. At least four noncoplanar corresponding points are required. Prefer more well-distributed points and independent validation; a low residual from a four-point affine fit alone does not demonstrate generalization.

```bash
python scripts/fit_registration.py --landmarks registration/landmarks.csv --output-affine registration/bregma_to_fem_affine.txt --output-qc registration/registration_qc.json --maximum-rms-mm 0.20
```

**Step 3: Review the configuration and solve.** Check the coordinates, currents, frequencies, conductivities, contacts, and output directory in `config_advanced.yaml`.

```bash
python scripts/run_advanced_fem.py --config config_advanced.yaml
```

When mesh and registration files are detected, the web application also enables **“运行完整注册 FEM”** (“Run full registered FEM”). This writes the current interface parameters to `outputs_advanced/web_run_config.yaml` and launches the solver. A `ready` indicator means files exist; their contents are checked during calculation.

Main outputs:

- `outputs_advanced/ti_vpl_efmouse_registered.vtu`: both channel fields, TI envelopes, and ROI cell data; view with ParaView.
- `outputs_advanced/ti_vpl_efmouse_registered_metrics.json`: registration, electrode projections, ROI statistics, current balance, and linear residuals. The web interface can download the latest result JSON.

The standard entry point stops for missing registration, excessive registration RMS or electrode projection distance, or a VPL ROI without gray-matter elements. `--demo-unregistered` is for pipeline debugging only, not formal inference.

**Current ROI limitation:** the advanced script still samples the 100 µm annotation through a generic bregma-to-CCF mapping. Its `Allen_VPL` selection uses structure ID 718 without explicitly restricting the left hemisphere as the web application does. Passing bregma-to-FEM registration does not establish subject-specific Allen atlas registration. Verify hemisphere, axis order, and ROI alignment before interpreting left VPL statistics.

## 6. Tests, troubleshooting, and scope

```bash
python -m pip install -e . pytest
python -m pytest
```

Tests cover default coordinates, electrode geometry, envelope properties, linear current scaling, a single-tetrahedron solution, and affine recovery. CI also checks Python syntax. **Passing CI establishes installation and existing test success, not animal or realistic head-model validation.**

| Issue | Action |
|---|---|
| Plots do not update after editing parameters | Submit the parameter form and check the Run number |
| FEM button is disabled | Download/convert the mesh and provide real registration files |
| `No module named nrrd` | Install `pynrrd` in the active environment |
| `Multiple top-level packages discovered` | Update the repository; `pyproject.toml` explicitly limits package discovery to `ti_vpl` |
| No FEM mesh on the hosted app | Use the online preview; configure the full workflow locally |

Research use still requires mesh convergence, contact/tissue parameter sensitivity analyses, registration/ROI checks, and independent software or measurement comparisons. The program does not establish neural activation, efficacy, or safe stimulation doses. [Chinese supplementary instructions](README_CN.md) provide additional operational background; use this page and the current source code for implementation details.

## 7. Data and method references

- [Allen CCF documentation](https://brain-map.org/support/documentation): CCFv3, 2017 annotations, and templates.
- [Wang et al. (2020), Cell](https://doi.org/10.1016/j.cell.2020.04.007): *The Allen Mouse Brain Common Coordinate Framework: A 3D Reference Atlas*.
- [Official EFMouse repository](https://github.com/klabhub/EFMouse) and [Sanchez-Romero, Akyuz & Krekelberg (2025), PLOS Computational Biology](https://doi.org/10.1371/journal.pcbi.1013471): mesh and original toolbox background.
- This repository implements its analytical preview, Python FEM, and envelopes in `ti_vpl/analytic.py`, `ti_vpl/fem.py`, and `ti_vpl/envelope.py`, respectively. Follow the licenses and citation requirements of external resources.
