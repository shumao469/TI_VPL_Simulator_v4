# Computational and rendering methods

## Scope

All six figures depict the existing v4 homogeneous-conductor analytic preview sampled on Allen CCF anatomy. They are deterministic model calculations, not animal measurements, neural simulations or registered multi-tissue FEM results. No empirical outcome, confidence interval or statistical significance is inferred from voxel counts.

## Geometry and anatomy

Coordinates are ordered (AP, ML, DV), in millimetres, with anterior, right and dorsal positive. The generic atlas mapping is

\[
(AP,ML,DV)=(5.4-0.1i,\;0.1k-5.7,\;-0.1j),
\]

where (i,j,k) are voxel indices in the bundled array ordered AP–DV–ML. The constants and voxel size are in mm. The mask contains 505,359 labelled brain voxels and 441 left VPL voxels at 100 µm spacing. This generic mapping does not establish registration to the experimental animal. The requested point is (-1.94,-1.75,-3.10) mm, and the mean of VPL voxel centres is (-1.4222,-1.9739,-4.2642) mm. VPL-centred sections are at the rounded coordinates AP=-1.40, ML=-2.00 and DV=-4.30 mm.

E1/E2 form the AP-oriented current pair and E3/E4 the ML-oriented pair. All electrode centres have DV=+0.50 mm. Their AP/ML coordinates follow the symmetric 2.08 mm half-spacing montage in README.md. Positions are ideal coordinates, not projected contacts on an anatomically reconstructed skull. A marker on a rendering does not imply a modelled contact area.

## Analytic fields

For channel c, the field at position r is

\[
\mathbf E_c(\mathbf r)=\frac{g I_c}{4\pi\sigma}
\left[
\frac{\mathbf r-\mathbf r_{c,+}}{\max(\|\mathbf r-\mathbf r_{c,+}\|,\rho)^3}
-\frac{\mathbf r-\mathbf r_{c,-}}{\max(\|\mathbf r-\mathbf r_{c,-}\|,\rho)^3}
\right].
\]

Here r and the source/return positions r(c,+/-) are in metres; I(c) is the per-channel current amplitude in amperes; sigma=0.275 S/m is uniform conductivity; rho=0.45 mm is the inherited near-source regularization radius; and g=0.028 is v4's dimensionless empirical gain. E(c) is returned in V/m. The gain is not a validated skull model. The reference currents are 0.50 mA per channel; the carrier frequencies are 2000 and 2010 Hz. Anatomical masks determine sampled/reporting locations only: no brain or skull boundary condition is imposed by this equation.

## Temporal-interference metric

For nonzero frequency difference the reported modulation depth is

\[
m(\mathbf r)=\max_{\|\mathbf n\|=1}
2\min\{|\mathbf n\cdot\mathbf E_1(\mathbf r)|,
       |\mathbf n\cdot\mathbf E_2(\mathbf r)|\}.
\]

The unit vector n represents a field projection direction. This maximization supplies an orientation-wise upper bound on projected envelope modulation depth. It does not specify axon orientation, neuronal polarization or recruitment. The closed-form implementation was checked against a dense independent search over directions. With one channel disabled or equal carrier frequencies, the preview reports zero difference-frequency modulation. No conclusion about high-frequency carrier stimulation follows from that convention. The formula's role follows the TI modelling framework of [Grossman et al. (2017)](https://pmc.ncbi.nlm.nih.gov/articles/PMC5520675/).

For Figure 05 the direction is the equal-current target-field bisector. With a=n·E1 and b=n·E2, the projected carrier signal is a cos(2πf1t)+b cos(2πf2t). Its envelope is sqrt(a²+b²+2ab cos(2π(f2-f1)t)), with time t in seconds and f in Hz. The curves are divided by a+b for this positive-projection configuration. Figure 05's temporal normalization is explicitly separate from the spatial-map normalization. It is not a demonstration of neuronal demodulation.

## Sampling and colour normalization

Every spatial map is divided by one baseline whole-brain maximum, 3.0793786491321113 preview V/m. It is evaluated over the fixed 100 µm brain mask using equal currents. The denominator is not recomputed for each panel or scenario. The logarithmic colour scale is 0.01–1 for mTI and 0.01–3 for the channel comparison; colour-bar extension triangles denote values beyond the display range. Raw, unnormalized values remain available in NPZ files. A plot does not imply that the VPL contains the global maximum.

Plane fields are evaluated directly from the analytic equations; brain and VPL masks are sampled by nearest atlas voxel. The plotting grid is denser than 100 µm for appearance, which does not increase the atlas's anatomical resolution. Three-dimensional cutaway views show a field section and a translucent atlas surface; they are not maximum-intensity projections or full volume renderings. VPL gold surfaces represent anatomy, not a selected field isosurface.

For appearance only, mesh vertices undergo eight neighbour-averaging steps with weight 0.15. All numerical field samples, voxel masks, ROI statistics and the unsmoothed source meshes remain unchanged. Orthographic projection is used. A 2 mm scale bar is defined in the projected image plane. All six figure PDFs/SVGs use vector polygons and curves with embedded/editable text; the PNGs are raster exports at 600 dpi.

## Summaries and tests

ROI values are simple voxel means, medians and empirical percentiles on equal-volume voxels; the stated voxel volume is 0.001 mm³. The figure 06 empirical CDFs describe these voxel distributions, not sampling uncertainty. The current sweep varies I1 from 0.05 to 0.95 mA and sets I2=1-I1 mA. It does not optimize electrode positions or estimate treatment efficacy.

Eight numerical tests cover masks/coordinates, linear current scaling, nonzero-frequency invariance, zero modulation limits, direction-search agreement, exact vector limits, preservation of the requested point and common normalization. These checks establish internal computation consistency only. They do not constitute external FEM or experimental validation.

Anatomy derives from the [Allen Mouse Brain CCF](https://pubmed.ncbi.nlm.nih.gov/32386544/). Multi-tissue follow-up can use [EFMouse](https://github.com/klabhub/EFMouse), with a separately verified registration and field calculation. The present figure package makes no claim that such a solve was performed.
