# Frontal area density of mangrove pneumatophores from hand-held lidar: method, synthetic validation and known biases

Technical report on the vertical-cylinder density estimate of raycloudtools (`--density_model vertical_cylinder` in
`rayrender` and `raydensity`). It describes the measurement model, the corrections, the cell outputs, a synthetic
validation, the known biases and what has not been tested, and gives a protocol for calibration and validation
with real data. Tiling, the Zarr output and large domains are described in
[`large_scale_density.md`](large_scale_density.md).

## Abstract

Hydrodynamic models of flow through mangrove pneumatophores need the frontal area per unit volume, $`a = N d`$ (stems
per m² of ground times stem diameter, m⁻¹), in cells of about 1 × 1 × 0.1 m over forests of about 300 × 300 m. We
estimate $`a`$ from hand-held SLAM lidar ray clouds by treating each ray as a probe of a medium of vertical cylinders:
the interception rate per metre of ray is $`a \sin\theta`$, which gives a $`G`$-weighted Beer–Lambert maximum-likelihood
estimate. Two corrections make it usable for 1 cm stems: (i) a Gaussian beam footprint model, which widens thin
stems to an effective width $`d_\mathrm{eff}(r)`$ of 2–3 cm and would otherwise overestimate $`a`$ by 2.1–2.6 times; and (ii) an
occlusion (space-occupancy) correction for non-overlapping stems, from a deterministic Monte Carlo response table
$`R(\phi, \rho, g)`$, applied to hits and path pooled over a horizontal neighbourhood of voxels. Voxel estimates are averaged
into cells (`density`, the main output) and also pooled (`pooled_density`), and their ratio flags patchiness.

On synthetic scans of randomly placed stems (30–1000 stems/m², 5 cm voxels), `density` in the lowest 0.1 m reads
4–8% high. This is self-shadowing: a stem removes the free path behind it within its own pooling window but keeps
its own hits. To first order the bias is $`c\,d_\mathrm{eff}/W`$ with $`c = 0.4732`$ ($`W`$ the window width), whatever the number of rays
or the stem density, and the optional `--self_shadow_correction` (here "ss") divides each window by this factor. A
wider validation matrix (Section 7.9: scan geometry, range and pose noise, diameter error and distribution, stem
heights, return convention, spacing and window size; 199 configurations, 4–8 seeds) tested a ±2% criterion in the
lowest 0.2 m, over cells at least 0.9 observed. The default `density` failed every configuration (mostly 4–9%
high); ss passed 144 of 199. Its failures were returns ranged to the widened outline rather than the stem surface
(up to +6% at 0.1–0.2 m), 5 cm trajectory error (−3 to −11%), overlapping stems (which the occlusion correction does
not model), 2 cm stems and a bimodal diameter mixture (+2.2 to +2.7% at 0.1–0.2 m), and near-misses at 1000 stems/m²;
returns from stem tops cause most of the remaining excess above 0.1 m.

In clumped stands ss reads low (patches 0.92–0.98, stems along cable roots 0.76–0.92), because there the
self-shadowing had been cancelling an opposite bias of uneven windows: the default's accuracy in patches, and its
13% low reading along cable roots, rest on that cancellation. A single-voxel configuration ("sw": a derived
second-order self-shadow factor $`1/(1-x)`$ plus three small-sample fixes; Section 8.1.2) resolves the clustering
instead (cable roots 0.98, against 0.87 with the default and 0.81 with ss; patches 0.97–1.02), but reads low for dense
stems under sparse scans (1000 stems/m²: 0.75–0.93), because without the neighbour top-up the densest voxels go
unobserved.

**Recommended configuration: the hybrid** (`--self_shadow_correction --hybrid`; Sections 8.1.3 and 10). In one pass
it computes ss and sw, and each cell takes sw where sw observes at least 0.9 of the cell's voxels, and ss otherwise.
Its errors in simulation (`density` / truth in the lowest 0.1 m, 4–8 seeds):

* uniform stands: within ±2% wherever ss passes; the worst case is 1.019 ± 0.004 (0.5 cm stems with a quarter of the
  standard rays);
* patches: 0.99–1.01;
* cable roots: 0.94–0.98 at 150 stems/m², and 0.43–0.71 at 500 stems/m² (0.47–0.75 over the cells at least 0.9
  observed);
* dense, poorly penetrated stands: the same as ss, including ss's shortfall where few cells are well observed;
* noise: it inherits sw's sensitivity, reading 0.990 with 1 cm range noise, 0.966–0.970 with 2 cm range noise and
  0.93–0.94 with 2 cm pose noise, where ss reads 0.998–1.008 and passes.

So where range or pose error exceeds about 1–2 cm, use ss, the robust fallback. In every configuration use only
cells with `observed_fraction` ≥ 0.9 (`--min_observed_fraction 0.9` blanks the others). `pooled_density` is within 3%
for uniform stems but about half the truth where stems are patchy within the cell.

The estimate is about as sensitive to the assumed stem diameter (−10% / +8% with ss for ±20% at 300 stems/m²) as to
the beam detection threshold (−11 to −12% / +22 to +27% for 0.4 / 3 times the true threshold), and to the minimum
gap between stems at high density (up to +34%). A beam footprint other than the one assumed acts only through the
error it makes in the path-weighted mean effective width at the ranges scanned (the density error is about −1.0
times that error with ss; Section 7.10). Common datasheet misreadings cost +18% / −21% (exit diameter) and +10% /
−15% (divergence), so without a rod calibration of the beam (Section 9.2) the footprint is the largest error source
for uniform stands. Section 9.3 ranks what a field campaign must measure. Layers above 0.1 m depend on how stem tops
and partial returns are treated by 5–20%, and
poorly penetrated cells (observed fraction below about 0.9) read low. All results are synthetic, from scenes
generated with the same idealised beam and stem models that the corrections assume; the method has not been checked
against field data, and Section 9 gives a procedure for doing so.

## Contents

1. [Introduction and target quantity](#1-introduction-and-target-quantity)
2. [Measurement model](#2-measurement-model)
3. [Beam footprint and effective width](#3-beam-footprint-and-effective-width)
4. [Estimators](#4-estimators)
5. [Occlusion and space-occupancy correction](#5-occlusion-and-space-occupancy-correction)
6. [Cell aggregation](#6-cell-aggregation)
7. [Synthetic validation](#7-synthetic-validation)
8. [Known biases and limitations](#8-known-biases-and-limitations)
9. [What was not simulated, and a calibration and validation protocol](#9-what-was-not-simulated-and-a-calibration-and-validation-protocol)
10. [Recommended settings and workflow](#10-recommended-settings-and-workflow)
11. [Equations and code](#11-equations-and-code)
12. [Reproducibility](#12-reproducibility)

## 1. Introduction and target quantity

Pneumatophores are the vertical aerial roots of mangroves such as *Avicennia*: thin (of order 1 cm), short (of order
0.1–0.5 m) and dense (hundreds to thousands per m²). For flow and wave models they are an array of vertical
cylinders, characterised per unit volume by

```math
a(z) = N(z)\,d, \qquad \phi(z) = N(z)\,\frac{\pi d^2}{4} = \frac{\pi d}{4}\,a(z) \qquad (1)
```

where $`N(z)`$ is the number of stems per m² of ground that reach height $`z`$, $`d`$ the stem diameter, $`a`$ the frontal area per
unit volume (m⁻¹) and $`\phi`$ the solid (volume) fraction. The drag force per unit volume of a canopy model is
proportional to $`a`$, and the frontal area index is $`\int a\,dz`$. For stems placed independently at random (a Poisson
process, so stems may overlap), the probability that a horizontal line of length $`L`$ misses every stem is
$`\exp(-a L)`$. Real stems don't overlap; for a non-overlapping (hard-core) arrangement the interception rate per unit
of free path is larger by $`R(\phi, 1)`$ (Section 5), so the gap fraction is approximately $`\exp(-a R(\phi, 1) L)`$: the exponent
is 3% larger at 300 and 10% larger at 1000 stems/m² of 1 cm stems.

The target is $`a`$ in cells of the hydrodynamic grid, e.g. 1 × 1 × 0.1 m (300 × 300 × 6 cells for a 300 m domain and
the lowest 0.6 m), from hand-held SLAM lidar. A ray cloud stores each return (or the end of a ray with no return)
with the sensor position at that time, so both where the beam stopped and the free space it crossed are known. The
estimator works on a voxel grid (5 cm by default) and aggregates the voxels into cells.

## 2. Measurement model

**Rays.** Each ray runs from the scanner position to its end point. A ray with nonzero alpha ended on a return
(a hit); a ray with zero alpha either had no return (its end is at the maximum range) or its return was removed
(e.g. classified as ground). Rays are clipped to the grid (`raylib/raydensitygrid.cpp`,
`DensityGrid::calculateDensities`), and the voxels they cross are walked (`DensityGrid::operator()` in
`raylib/raydensitygrid.h`).

**Interception by vertical cylinders.** A ray at zenith angle $`\theta`$ (measured from vertical) moves $`\sin\theta`$ horizontally
per metre of ray. It intercepts a vertical cylinder of diameter $`d`$ when the stem's centre lies within $`d/2`$ of the
ray's horizontal track, so for $`N`$ stems per m², placed at random, the interception rate per metre of ray is

```math
\lambda(\theta) = G(\theta)\,a, \qquad G(\theta) = \sin\theta. \qquad (2)
```

This replaces the default $`G = 0.5`$ of raycloudtools, which applies to leaves with a spherical (isotropic) angle
distribution and gives leaf area per volume. For vertical stems it would bias the estimate by a factor $`2 \sin\theta`$,
about 2 for the near-horizontal rays of a hand-held scanner (test `RayDensityVerticalCylinders`: at $`\theta = 90^\circ`$ the
estimate was 1.003 of the truth with $`G = \sin\theta`$, and 2.005 with $`G = 0.5`$). The end caps (stem tops) are not part of
the model (Section 8.2).

**Likelihood.** If hits along a ray form a Poisson process with rate $`G a`$, a ray $`i`$ that crosses a voxel with path
length $`L_i`$ inside it contributes the log-likelihood $`h_i \ln(G_i a) - a G_i L_i`$, where $`h_i`$ is 1 if it ends with a hit
in the voxel and 0 otherwise. Maximising over $`a`$ gives the $`G`$-weighted Beer–Lambert estimate

```math
\hat a = \frac{H}{\sum_i G_i L_i} \qquad (3)
```

with $`H`$ the number of hits in the voxel (or neighbourhood, Section 4). With the beam footprint (Section 3), $`G_i`$ is
replaced by $`G_i \rho_i`$.

**Path and hit rules.**

* *Path to the end.* A ray that ends in a voxel (hit or not) counts its path up to its end point there; a ray that
  passes through counts the full chord.
* *Misses are censored observations.* Rays with no return, and ground returns given zero alpha, add path but no
  hit. Ground returns must therefore be removed (set to zero alpha) before the estimate, or they count as stem hits.
* *Clipped rays.* A ray whose end point lies outside the grid has passed through the grid without a return inside
  it, so it is a miss there. Its start is also clipped, but the range from its original start is kept, as the
  footprint correction needs it (`range_offset_`).
* *Zero-length passes.* A ray that only grazes a voxel edge, corner or face adds no path and is not counted as one
  of the voxel's rays (this keeps the ray counts identical between tiled and single runs).
* *Ray angle selection.* Rays with $`\sin\theta`$ below `--min_sin_zenith` (default 0.1) or above `--max_sin_zenith`
  (default 1) are ignored. Near-vertical rays carry almost no information on vertical stems, and selecting rays on
  direction alone does not bias the estimate.

## 3. Beam footprint and effective width

A lidar pulse has a finite, range-dependent spot. We model it as a Gaussian beam of 1/e² diameter

```math
w(r) = w_0 + \gamma\,r \qquad (4)
```

with exit diameter $`w_0`$ and full divergence $`\gamma`$ (both at the 1/e² level, from the scanner's specification), so its
standard deviation is $`\sigma = w/4`$. A vertical cylinder is seen from any horizontal direction as a strip of width $`d`$. If
the beam axis passes at horizontal offset $`x`$ from the stem's centre, the fraction of the beam's power that the
stem intercepts is

```math
F(x) = \tfrac12\left[\operatorname{erf}\!\left(\frac{x + d/2}{\sqrt2\,\sigma}\right) - \operatorname{erf}\!\left(\frac{x - d/2}{\sqrt2\,\sigma}\right)\right] \qquad (5)
```

A return is triggered when $`F \ge \tau`$, the detection threshold (a fraction of the beam power). The stem is therefore
detected over an effective width

```math
d_\mathrm{eff}(r) = 2x^*, \quad F(x^*) = \tau, \qquad \rho(r) = d_\mathrm{eff}(r)/d \qquad (6)
```

found by bisection (`raylib/raystemmodel.cpp`, `stripInterceptFraction`, `effectiveStemWidth`), and zero when
$`F(0) < \tau`$ (a stem too thin ever to be detected). For 1 cm stems, a 10 mm, 3 mrad beam and $`\tau = 0.05`$, $`d_\mathrm{eff}`$ is 2.07 cm
at 1 m and 3.05 cm at 5 m; the path-weighted mean $`\rho`$ over the voxels of the synthetic scans was 2.13–2.42. The
interception rate becomes $`N d_\mathrm{eff} \sin\theta = a \rho G`$, so each ray segment is weighted by $`G \rho(r)`$, with $`r`$ the range of the
segment's midpoint from the ray's original start (`StemWidthModel`, tabulated every 1 cm up to 100 m and on a
logarithmic grid beyond), and the estimate stays $`a = N d`$ for the true mean diameter. Without this correction the
synthetic 1 cm stems were overestimated 2.1–2.6 times.

**Inputs.** $`w_0`$ and $`\gamma`$ come from the scanner's beam specification. The threshold $`\tau`$ is not a published parameter: it
depends on the detector, the return mode, the target's reflectance and the range, and has to be calibrated
(Section 9). The diameter $`d`$ is an input: it should be the mean diameter at the heights analysed (stems taper). The
result is about as sensitive to $`d`$ as to $`\tau`$ (Section 7.5), because with a 13–25 mm footprint the hit rate is roughly
proportional to $`N (d + 1\text{–}2\,\mathrm{cm})`$, so the data mainly measure the stem count, and $`N d`$ inherits about half of any
relative error in $`d`$.

**Return range convention.** When the beam axis misses a stem but enough power hits it, where along the beam is the
return reported? `--return_range` selects the model, in the generator and in the correction:

* `surface` (default): at the range of the physical stem surface (the nearest illuminated part of the stem), or of
  the beam axis's closest approach to the stem when the axis misses it. A time-of-flight scanner ranges to the
  illuminated surface, so this is the physically reasonable default.
* `outline`: where the beam axis enters the footprint-widened outline (diameter $`d_\mathrm{eff}`$).

In the generator the two conventions return the same rays (the side hit is tested where the beam axis enters the
widened outline), and with `surface` each partial return is reported 6.4–6.9 mm further along the ray
(horizontally), i.e. 5–13 mm lower for descending rays, and that ray's path is longer by as much. The occlusion
response is nearly the same for both (Section 5), and estimating a scene with the other convention changed the
estimates by at most 0.2%. But the generated data differ where stems end, by 1–15% depending on height
(Section 8.2). Which convention a real scanner follows has not been verified; Section 9 describes a check.

## 4. Estimators

Voxels of 5 cm are comparable to the stem spacing (about 6 cm at 300 stems/m²), so a voxel holds zero, one or a few
discrete stems, against the assumption of the standard per-voxel estimate (many small, randomly placed elements).
Three estimators are available (`--estimator`); the results quoted are plot totals in the lowest 0.1 m, over 30–1500
random stems/m² (Table 2).

**`beer_lambert`** (the raycloudtools default): Eq. (3) per voxel, multiplied by the small-sample factor $`(n-1)/n`$
of the original estimator, for the voxel's $`n`$ rays (`DensityGrid::Voxel::density`). It read 1.24–1.31 of the truth
in the lowest 0.1 m and 1.31–1.39 at 0.1–0.2 m (cell means). Most of this is the self-shadowing bias of
Section 8.1 at a window of one voxel ($`W`$ = 5 cm: +21% to first order, +28–29% measured with `--estimator pooled
--pool_radius 0`).

**`crossing`**: hit rays count their path to where they would have left the voxel, plus the mean chord of the
widened stem that remains beyond the return, divided by $`\sin\theta`$:

```math
\Delta L = \frac{1}{\sin\theta}\times\begin{cases}\pi d_\mathrm{eff}/4 & \text{outline}\\ \pi (d_\mathrm{eff}^2 + d^2)/(8 d_\mathrm{eff}) & \text{surface}\end{cases} \qquad (7)
```

(`DensityGrid::calculateDensities`). Without a footprint, density = hits / crossing length does not assume random
placement within the voxel. But it has no occlusion correction, and with a footprint the entries into a widened
stem that lie inside another widened stem are never observed. It read 1.05 at 30, 0.97 at 300, 0.82 at 1000 and
0.74 at 1500 stems/m², lower than with a near-perfect beam (1 mm, 0.1 mrad: 1.00 at 300 and 0.93 at 1000). It is not
recommended with a footprint.

**`pooled`** (recommended): the Beer–Lambert hits and weighted path are summed over a horizontal neighbourhood
$`\mathcal N(v)`$ of each voxel $`v`$ before taking the ratio,

```math
m_v = \frac{\sum_{u\in\mathcal N(v)} H_u}{\sum_{u\in\mathcal N(v)} \sum_i G_i \rho_i L_i}, \qquad \bar\rho_v = \frac{\sum_{u}\sum_i \rho_i\,(G_i \rho_i L_i)}{\sum_u \sum_i G_i \rho_i L_i} \qquad (8)
```

and $`m_v`$ is corrected for occlusion (Section 5) using the neighbourhood's mean widening $`\bar\rho_v`$. The voxel takes the
corrected value $`a_v`$ and keeps its own path length (so `pooled_density`, Section 6, weights it by its own path).
Voxels that no ray reached are not assigned a value. $`\mathcal N(v)`$ is the $`(2k+1) \times (2k+1)`$ voxels centred on $`v`$ in its
horizontal layer, with $`k`$ = round(`--pool_radius` / voxel width), so it spans $`\pm(k + \tfrac12)`$ voxels and its width is
$`W = (2k+1) w`$: 0.15 m (±0.075 m) for the default 0.05 m at 5 cm voxels, 0.25 m for 0.1 m, and 0.30 m for 0.05 m at
10 cm voxels (the run prints the window). It changes in steps with the voxel size. The window is computed with
summed-area tables (`DensityGrid::applyOverlapCorrection`).

Why pool: the Beer–Lambert estimate is unbiased for randomly placed stems however few there are in a voxel, provided
hits and path are counted over a region, not averaged as ratios of single voxels; pooling moves the ratio from one
voxel's few stems to a neighbourhood's, and the occlusion correction handles non-overlapping stems. Pooling trades
resolution for bias: it smooths over the window, so patch edges spread by about the window, stems clustered more
tightly than the window read low (Section 8.6), and the remaining self-shadowing bias falls as the window grows
(Section 8.1). Pooled corrected values were 1.01–1.08 of the truth (cell means; 0.98–1.01 path-weighted).

With `beer_lambert` and the occlusion correction on (a stem diameter given), the neighbourhood ratio is used only to
compute the correction factor, which is applied to the voxel's own estimate. With `pooled`, `--self_shadow_correction`
also divides $`m_v`$ by a self-shadow factor before the occlusion correction (Section 8.1, Eq. 12), and four further
options, all off by default, adapt the estimator to single-voxel windows (`--pool_radius 0`): the second-order
self-shadow factor (`--self_shadow_second_order`), the occlusion response inverted over a larger neighbourhood
(`--occlusion_radius`), a finite-ray correction of each window's ratio (`--ratio_bias_correction`) and no
neighbour top-up (`--no_neighbour_priors`); Section 8.1 describes them and the configuration "sw" that combines them.
`--hybrid` (with `--self_shadow_correction`, and not with the four single-voxel options, which it sets itself)
computes the default-window and the single-voxel estimates in the same pass and chooses between them per cell
(Section 8.1.3).

## 5. Occlusion and space-occupancy correction

The pooled estimate measures hits per unit of free path. For randomly placed, overlapping discs (a Boolean model)
the free path to the next hit is exponential, and the hit rate per unit free path is exactly $`N d_\mathrm{eff}`$, however much
the widened discs overlap. Real stems don't overlap. In a non-overlapping (hard-core) arrangement the
footprint-widened stems leave less free space than overlapping ones would, so the measured rate exceeds $`N d_\mathrm{eff}`$.
With a perfect beam ($`d_\mathrm{eff} = d`$) it is at least $`a/(1-\phi)`$, because the stems' own cross-sections are excluded from
the free path, and a little more because rays enter the stem layer at random points, so sample the gaps between
stems in proportion to their length, and hard-core gaps are more regular than exponential.

**Response.** We define the response as the measured rate relative to $`N d_\mathrm{eff}`$,

```math
R(\phi, \rho, g) = \frac{\text{hits per unit free path}}{N\,d_\mathrm{eff}}, \qquad m = a\,R\!\left(\tfrac{\pi d}{4} a,\ \bar\rho,\ g\right) \qquad (9)
```

where $`g`$ is the minimum gap between stem surfaces in diameters (stem centres at least $`d(1+g)`$ apart; `--min_gap`).
$`R`$ is computed by a 2D Monte Carlo simulation (`raylib/raystemmodel.cpp`, `freePathResponse`): in units of $`d`$,
1500 non-overlapping discs are placed by random sequential addition in a square sized for $`\phi`$; 10 000 lines start
from uniformly random points of the free space in the central 80% of the square, in random directions, and run to
their first reported return or the edge of the square; $`R = \sum\text{hits} / (\sum\text{path}\cdot N \rho)`$, pooled over 64 independent
realisations. The free space follows the return convention: under `outline`, outside the widened outlines; under
`surface`, outside the physical stems, except inside a widened outline beyond its reported return (a ray there
would already have returned). A response is saturated (infinite) when under 1% of candidate start points are free.

**Table.** For $`g = 0`$, $`R`$ is tabulated for both conventions over 13 solid fractions (0, 0.01, 0.02, 0.04, 0.07, 0.1,
0.13, 0.16, 0.2, 0.25, 0.3, 0.35, 0.4) and 17 width ratios (1 to 5 in steps of 0.25), embedded in the library
(`raylib/raystemresponsetable.h`, generated by `scripts/generate_stem_response_table.cpp`) and interpolated
bilinearly (`StemResponseTable::response`). For $`g > 0`$ the table is over the solid fraction of the exclusion discs,
$`\phi(1+g)^2`$, and its nodes are simulated when first needed. Examples: $`R(0.024, 2.3, 0) = 1.046`$ (300 stems/m² of 1 cm
stems), $`R(0.079, 2.3, 0) = 1.157`$ (1000 stems/m²) and $`R(0.2, 2.3, 0) = 1.51`$. With a perfect beam $`R(\phi, 1, 0)`$ is 1.03
and 1.10 at 300 and 1000 stems/m², so about a third of the excess above 1 comes from the footprint widening.

**Accuracy of the table.** Rerunning the simulation with 6 other seeds varied $`R`$ by 0.2–0.5% (s.d.) at $`\phi \le 0.2`$, and
by 0.6% (surface) and 0.9% (outline) at $`\phi = 0.25`$. Interpolation differed from direct simulation by at most 0.2% for
$`\phi \le 0.2`$, and 0.4% at $`\phi = 0.3`$. The two conventions agree within 0.2% at $`\phi \le 0.1`$ and 0.7% at $`\phi = 0.2`$, and diverge
above: 1.6% at $`\phi = 0.25`$ ($`\rho = 3`$) and 3.8% at $`\phi = 0.3`$ ($`\rho = 2.1`$), surface lower.

**Inversion.** For each voxel, Eq. (9) is solved for $`a`$ given the neighbourhood's $`m`$ and $`\bar\rho`$ by bracketing upwards
from $`m`$ (doubling, so $`R < 1`$ would also be handled) and bisecting (`invertDensityResponse`, `invertStemResponse` in
`raylib/raystemmodel.h`/`.cpp`).

**Limits.** The exclusion-disc solid fraction is capped at 0.4 (below the 0.547 random sequential addition limit
for discs), so $`N d \le 1.6/(\pi d (1+g)^2)`$: 51 m⁻¹ (5100 stems/m²) for 1 cm stems with no gap. Voxels beyond this are
set to the limit, and their number is printed. In the synthetic scans the correction lowered the estimates: the
uncorrected cell means were 1–6% higher at 30–300 stems/m², 11% at 600 and 18–30% at 1000–1500.

**Minimum gap.** With stems at least one diameter apart, assuming no gap overestimated the cell mean by 34% at 1000
stems/m², 10% at 300 and 4% at 100 (half a diameter apart: 15%, 5% and 2%), relative to assuming the right gap
(Table 6). $`g`$ should be set from field observations of how closely stems pack.

**Determinism.** The simulation uses fixed seeds and its own uniform random numbers from `std::mt19937`, so $`R`$, the
table and every estimate are the same on any platform, and a tiled run equals a single run (Section 7.7).

## 6. Cell aggregation

**Assignment.** Each voxel belongs to the cell that contains its centre (`aggregateDensityCells` in
`raylib/raydensitygrid.cpp`). With `rayrender`, the voxel and cell grids share the min bound, so cells that are
multiples of the voxel align exactly. With `raydensity --ground`, a voxel's height is that of its centre above the
ground mesh beneath the centre of its column (the mesh sampled at voxel resolution), so layers are misregistered by
up to half a voxel (2.5 cm) where the ground height isn't a multiple of the voxel: keep voxels no larger than half
the layer height. Don't flatten the cloud instead (`raytranslate subtract`): it shifts each whole ray by the ground
height at its end, misplacing the rest of the ray by the slope times its horizontal length (Section 7.7).

**Observed voxels.** A voxel with fewer than 10 rays is topped up from its 3 × 3 × 3 neighbourhood (faces, then
edges, then corners, until it has 10 rays; `DensityGrid::addNeighbourPriors`), except that for vertical cylinders a
voxel that no ray reached is not filled, and with `--no_neighbour_priors` (pooled estimator) no voxel is topped up. A voxel that
still has 2 or fewer rays after the top-up is unobserved. In
practice the unobserved voxels are those no ray reached (below the ground, beyond the scan, or completely hidden),
or whose neighbours are also nearly empty.

**Outputs.** Over the cell's $`n_\mathrm{obs}`$ observed voxels, with corrected voxel densities $`a_v`$ and weighted path lengths
$`P_v = \sum_i G_i \rho_i L_i`$:

```math
\texttt{density} = \frac{1}{n_\mathrm{obs}}\sum_v a_v, \qquad \texttt{pooled\_density} = \frac{\sum_v a_v P_v}{\sum_v P_v}, \qquad \texttt{patchiness\_ratio} = \frac{\texttt{density}}{\texttt{pooled\_density}} \qquad (10)
```

* `density` (main output): the mean of the observed voxel densities, each voxel weighted equally (by volume), so
  the volume average that a hydrodynamic cell needs. It follows patchy stems, but carries the self-shadowing bias
  of Section 8.1 (+4–8% for random stems in the lowest 0.1 m at the default settings) unless
  `--self_shadow_correction` is given.
* `pooled_density`: the cell's corrected hits over its weighted path, i.e. the path-weighted mean of the voxel
  densities. Rays penetrate dense parts less, so they get less weight, horizontally and vertically: it is within 3%
  only where the stems are uniform within the cell in both directions, about half the truth for patches within the
  cell, and low in layers where the density falls with height (Section 8.5).
* `patchiness_ratio`: in the lowest 0.2 m, 1.03–1.09 for random stems at the default pool (mostly the bias of
  `density`; 0.97–1.02 with `--self_shadow_correction`), 1.6–2.4 for patches and cable roots. Higher up, the
  vertical profile within the cell also raises it (up to 1.17 for random stems at 0.2–0.3 m).
* `hits` ($`\sum_v a_v P_v`$, the corrected hits), `weighted_path` ($`\sum_v P_v`$) and `observed_fraction` ($`n_\mathrm{obs}`$ over the cell's
  voxels). The profile (`--profile`) gives the same per horizontal layer of voxels, with `observed_voxels`,
  `total_voxels` and `mean_hits_per_voxel`.
* `hybrid_weight` (with `--hybrid` only; a column of `--cells` and `--profile`, and an array of the Zarr store): the
  fraction of the cell's observed voxels whose density is the single-voxel estimate. It is 0 where the cell takes
  the default-window estimate (ss), and at least 0.9 where it takes the single-voxel one (Section 8.1.3).

A cell without observed voxels is NaN (with `observed_fraction` 0 in the Zarr store). Both densities use only the
observed voxels, so they describe the observed part of the cell: where much of the cell is unobserved, they read
low (Section 8.3).

**Observed-fraction mask.** `--min_observed_fraction f` (in `rayrender`, for `--cells` and `--profile`, and in
`raydensity`; default 0, off) sets `density`, `pooled_density` and `patchiness_ratio` (and `solid_fraction` in the
store) to NaN in cells or layers whose `observed_fraction` is below f, and keeps `observed_fraction`, `hits`,
`weighted_path` and `hybrid_weight`, so the masked cells can still be inspected or aggregated. 0.9 is recommended
(Section 8.3). `raydensity` records the value in the store's `.zattrs` (`"min_observed_fraction": f`, only when set).

**Coarser cells.** Cells can be aggregated further from the stored values: `pooled_density` from the summed `hits`
and `weighted_path`, and `density` as the mean weighted by the number of observed voxels (`observed_fraction` times
the voxels per cell). When the cell dimensions are multiples of the voxel width, every cell holds the same number
of voxels (also with `--ground`, as the grid spans the full relief), so weighting by `observed_fraction` is exact;
otherwise the voxel count per cell varies with the ground height, and isn't stored.

## 7. Synthetic validation

### 7.1 Generator and scenes

`raycreate pneumatophores <seed>` (`raylib/raypneumatophoregen.cpp`) generates vertical cylinders, by default of one
diameter (1 cm, `--diameter`; `--diameter_sd sd` gives log-normal diameters of that mean, and `--diameter_mixture
d1,d2,f` two diameters with a fraction $`f`$ of the second), on flat ground at $`z = 0`$, and scans them with a hand-held
lidar walked along a lawnmower path. The defaults, used for every result unless stated, are:

* Heights: normal, mean 0.25 m and s.d. 0.08 m (`--height mean,sd`; heights below 2 cm are redrawn; 3%, 24%, 47%
  and 24% of stems end in successive 0.1 m layers). `--height 0.6,0` gives stems of one height, and
  `--uniform_height min,max` heights uniform in a range.
* Arrangement: random, non-overlapping (random sequential addition with the exact number of stems; `--min_gap`
  for a minimum gap, `--overlapping` for Poisson positions); `--patchy radius,cover,matrix_ratio` (stems only in
  random circular patches); or `--clustered` (stems along cable roots radiating from trees, 0.04 trees/m², 12 roots
  of 6 m per tree, stems offset from the root lines with s.d. `--cluster_spread`).
* Scan: a 6 × 6 m plot (`--plot_width`; stems generated 8 m beyond it), scanner at 1.2 m height
  (`--scanner_height`) walking lines 1 m apart (`--walk_spacing`), 2 000 000 rays (`--rays`; about 56 000 per m² of
  plot) uniform over zenith angles 0–101.5° (`--zenith`), 30 m range for rays with no return; ground returns are
  given zero alpha (perfectly classified; `--ground_returns` keeps them as returns).
* Beam: Gaussian, 10 mm exit diameter, 3 mrad divergence, $`\tau = 0.05`$, `surface` return range (`--footprint`,
  `--return_range`); stem tops are solid discs of the widened outline. No intensity.
* Measurement errors: none by default. `--range_noise sd` adds normal noise to the range of each return;
  `--pose_noise sd,interval` offsets each ray by a smooth trajectory error, independent normal values per axis every
  `interval` metres of walk, interpolated linearly. They use their own random sequence, so the stems and the ideal
  rays are those of the noise-free scene.
* Truth: $`N d`$ in each 1 × 1 × 0.1 m cell (`--truth_cells`) and each 1 cm layer, from the stems' positions and
  heights.

Scenes depend only on the seed and options, on any platform.

**Stage-1 self-shadowing scenes.** The results of Sections 8.1 and 8.3 that compare the estimator with and without
`--self_shadow_correction` used smaller scenes: a 4 × 4 m plot (`--plot_width 4`) scanned with 56 250 rays per
metre of walk (0.9 M rays at 1 m line spacing, the ray density of the standard scene), 16 cells per layer, seeds
101–104, reported as the mean over seeds (standard errors of `density` 0.001–0.01 for well observed random stems,
up to 0.015 in poorly observed layers and 0.02 in patches and along cable roots).

**Stage-2 scenes.** The validation matrix (Section 7.9) and the small-window study (Section 8.1) used the same 4 m
plots and ray density, seeds 101–104 (the matrix re-ran its failing configurations with 105–108 as well). The
small-window study also scaled the rays to ×0.25 and ×0.5 of that density (`--rays_per_m`), to test sparse scans.

**Metric.** Each table gives the plot total Σ estimate / Σ truth over the 36 cells of a 0.1 m layer of the 6 × 6 m
plot, as the mean ± s.d. over independent scenes (seeds), estimated with `--estimator pooled`, 5 cm voxels, the
default pool radius (0.05 m) and the footprint and stem diameter of the scene, unless stated. "Cell RMSE" is the RMSE
of the cells relative to the mean true cell value.

### 7.2 Random stems

**Table 1.** Random stems, 8 scenes each (seeds 101–108).

| Stems/m² | `density` 0–0.1 m | 0.1–0.2 m | 0.2–0.3 m | `density`, `--pool_radius 0.1`, 0–0.1 m | `pooled_density` 0–0.1 m | 0.1–0.2 m | 0.2–0.3 m | No occlusion correction, `density` 0–0.1 m |
|---|---|---|---|---|---|---|---|---|
| 100 | 1.073 ± 0.006 | 1.085 ± 0.006 | 1.123 ± 0.006 | 1.042 ± 0.007 | 0.983 ± 0.007 | 0.995 ± 0.005 | 1.028 ± 0.007 | 1.099 |
| 300 | 1.065 ± 0.005 | 1.077 ± 0.003 | 1.113 ± 0.003 | 1.033 ± 0.006 | 0.979 ± 0.006 | 0.988 ± 0.004 | 0.999 ± 0.004 | 1.126 |
| 1000 | 1.040 ± 0.005 | 1.075 ± 0.003 | 1.109 ± 0.002 | 1.021 ± 0.006 | 1.010 ± 0.005 | 0.999 ± 0.004 | 0.952 ± 0.002 | 1.231 |

The s.d. between scenes is under 1%, so the differences between rows and columns are systematic, not sampling
noise. The cell RMSE of `density` was 5–8% in the lowest 0.1 m and 8–9% at 0.1–0.2 m, mostly its bias;
`pooled_density` had 2.5–4%.

**Table 2.** Random stems over a wider range of densities, and the other estimators, 4 scenes each (seeds 101–104),
lowest 0.1 m (0.1–0.2 m in brackets for `density` and `pooled_density`).

| Stems/m² | Observed fraction | `density` | `pooled_density` | Uncorrected `density` | `crossing` | `beer_lambert` (uncorrected) | Mean $`\rho`$ |
|---|---|---|---|---|---|---|---|
| 30 | 1.00 | 1.080 (1.089) | 0.987 (0.999) | 1.093 | 1.055 | 1.299 | 2.42 |
| 100 | 1.00 | 1.075 (1.087) | 0.984 (0.996) | 1.101 | 1.032 | 1.299 | 2.35 |
| 300 | 1.00 | 1.063 (1.076) | 0.977 (0.987) | 1.123 | 0.974 | 1.313 | 2.25 |
| 600 | 1.00 | 1.061 (1.079) | 0.995 (0.989) | 1.175 | 0.898 | 1.296 | 2.19 |
| 1000 | 0.93 | 1.039 (1.076) | 1.008 (0.999) | 1.230 | 0.815 | 1.240 | 2.15 |
| 1500 | 0.68 | 1.011 (1.060) | 0.981 (1.014) | 1.317 | 0.739 | 1.246 | 2.13 |

`density` reads 6–8% high wherever the cells are fully observed, independently of the density: this is the
self-shadowing bias, which `--self_shadow_correction` removes for random stems (Section 8.1). The lower values at
1000 and 1500 stems/m² are not better accuracy: there 7% and 32% of the lowest layer's voxels are unobserved, and
the low bias of poorly observed cells (Section 8.3) cancels part of the high bias.

### 7.3 Patches and cable roots

Pneumatophores are often patchy and are normally arranged in rows along the cable roots from which they rise.

**Table 3.** Non-random arrangements, lowest 0.1 m (0.1–0.2 m in brackets). Patches of 0.3 m: 8 scenes (seeds
101–108); others: 4 scenes (101–104). s.d. under 1% except `pooled_density` in patches (2–3%).

| Arrangement | `density` | `density`, `--pool_radius 0.1` | `pooled_density` | `patchiness_ratio` | Uncorrected `density` |
|---|---|---|---|---|---|
| Cable roots, 150 stems/m², 3 cm lateral spread | **0.873 (0.904)** | 0.767 (0.801) | 0.514 (0.564) | 1.76 (1.64) | 0.924 |
| Patches of 0.5 m radius, 30% cover, 300 stems/m² within | 1.028 (1.050) | 0.951 (0.986) | 0.471 (0.602) | 2.39 (1.85) | 1.082 |
| Patches of 0.3 m radius, 40% cover, 300 stems/m² within | 1.015 (1.040) | 0.925 (0.966) | 0.547 (0.655) | 2.05 (1.71) | 1.065 |
| Patches of 0.2 m radius, 30% cover, 300 stems/m² within | 0.991 (1.014) | 0.882 (0.923) | 0.551 (0.640) | 1.88 (1.64) | 1.036 |

Stems along cable roots read 13% low with the default pool and 23% low with `--pool_radius 0.1`, because they are
clustered more tightly (3 cm) than the window (Section 8.6). `pooled_density` is about half the truth in every
non-random arrangement, and lower still over a whole 6 m profile layer (0.30 for 0.5 m patches and 0.34 for cable
roots, one scene each). The near-unity `density` in patches, and the cable roots' −13% rather than more, are partly
a cancellation with the self-shadowing bias: with `--self_shadow_correction` these arrangements read 4–8% and 19%
low (Section 8.1). Single-voxel windows (sw, Section 8.1.2) resolve them instead: 0.98 along the cable roots and
0.99–1.00 in patches (stage-2 scenes).

### 7.4 Ray angles and ray count

**Ray angles** (4 scenes, all rays versus only near-horizontal rays, `--min_sin_zenith 0.9`):

| | 0–0.1 m | 0.1–0.2 m | 0.2–0.3 m | 0.3–0.4 m |
|---|---|---|---|---|
| 300/m², all rays, `density` | 1.063 | 1.076 | 1.114 | 1.173 |
| 300/m², all rays, `pooled_density` | 0.977 | 0.987 | 1.000 | 1.075 |
| 300/m², shallow rays, `density` | 0.959 | 1.023 | 1.080 | 1.105 |
| 1000/m², all rays, `density` | 1.039 | 1.075 | 1.110 | 1.156 |
| 1000/m², all rays, `pooled_density` | 1.008 | 0.999 | 0.952 | 1.006 |
| 1000/m², shallow rays, `density` | 0.353 | 0.674 | 0.875 | 1.066 |

Using only shallow rays lowered every layer at 300 stems/m², not only those where stems end, and at 1000 stems/m² so
few shallow rays reach the lowest layer that it read 0.35: angle selection is not a remedy for the stem-top bias
(Section 8.2).

**Ray count.** The same scene (300 stems/m², seed 101) with 0.5, 2 and 8 million rays gave `density` 1.064, 1.056
and 1.062 at 0–0.1 m and 1.073, 1.072 and 1.073 at 0.1–0.2 m: more scanning does not reduce the high bias.

### 7.5 Sensitivity to the beam threshold and the stem diameter

**Table 4.** Scenes generated with $`\tau = 0.05`$, estimated assuming $`\tau`$ (4 scenes, `density`, lowest 0.1 m).

| Stems/m² | Assume $`\tau = 0.02`$ | $`\tau = 0.05`$ (true) | $`\tau = 0.15`$ |
|---|---|---|---|
| 100 | 0.943 (−12%) | 1.075 | 1.358 (+26%) |
| 300 | 0.939 (−12%) | 1.063 | 1.325 (+25%) |
| 1000 | 0.929 (−11%) | 1.039 | 1.265 (+22%) |

Relative changes are against the true-$`\tau`$ estimate; at 0.1–0.2 m and for `pooled_density` they are similar (up to
+27%).

**Table 5.** Scenes of 1 cm stems estimated assuming another diameter (one scene, seed 101, `density`, lowest
0.1 m).

| Stems/m² | Assume $`d = 0.8`$ cm | $`d = 1`$ cm (true) | $`d = 1.2`$ cm |
|---|---|---|---|
| 300 | 0.941 (−11%) | 1.056 | 1.147 (+9%) |
| 1000 | 0.954 (−8%) | 1.038 | 1.094 (+5%) |

The effect is smaller at high density, because the occlusion correction partly compensates. In the footprint model
alone, $`d_\mathrm{eff}/d`$ for an assumed 0.8 or 1.2 cm changes by −11% / +9% at 1 m range and −14% / +12% at 5 m relative to
1 cm. The validation matrix (Section 7.9, block D) repeated this over 4 seeds, for 0.5, 1 and 2 cm stems and
assumed diameters of 0.7–1.3 times the true one (Section 9.3 gives the required accuracy), and generated diameter
distributions: with ss, log-normal diameters of 10 ± 3 mm read 1.004 / 1.017 (0–0.1 / 0.1–0.2 m), and an equal
mixture of 6 and 14 mm stems 1.008 / 1.022, against 1.001 / 1.010 for stems of one diameter, all estimated with the
10 mm mean (an earlier model calculation had put the mixture's effect below 1%). So a spread of diameters about the
right mean adds up to about 1%, but the mean must be right.

### 7.6 Sensitivity to the minimum gap, the return convention, the voxel size and the beam

**Table 6.** Scenes with a minimum gap $`g`$ between stems, estimated assuming $`g = 0`$ or the true $`g`$ (4 scenes,
`density`, lowest 0.1 m).

| Stems/m² | $`g = 0.5`$, assume 0 | assume 0.5 | $`g = 1`$, assume 0 | assume 1 |
|---|---|---|---|---|
| 100 | 1.089 | 1.069 | 1.101 | 1.059 |
| 300 | 1.107 | 1.059 | 1.148 | 1.048 |
| 1000 | 1.176 | 1.025 | 1.344 | 1.004 |

**Return convention.** Scenes generated with `outline` read higher than `surface` scenes, with or without the
correction, and the difference grows with height (Section 8.2): 1–2% at 0–0.1 m and 2.5–4.5% at 0.1–0.2 m (4 scenes;
e.g. `density` 1.063 → 1.078 at 300 and 1.039 → 1.057 at 1000 stems/m²), and 5–10% at 0.2–0.5 m (one scene).
Estimating with the other convention changed the estimates by at most 0.2%.

**Table 7.** Voxel size (random stems, 4 scenes, default pool radius 0.05 m, so windows of $`W = 0.125`$, 0.15 and
0.30 m; `density` and, in brackets, `pooled_density`).

| Stems/m² | Voxel | Observed fraction 0–0.1 m | 0–0.1 m | 0.1–0.2 m |
|---|---|---|---|---|
| 300 | 2.5 cm | 0.99 | 1.071 (1.029) | 1.092 (0.999) |
| 300 | 5 cm | 1.00 | 1.063 (0.977) | 1.076 (0.987) |
| 300 | 10 cm | 1.00 | 1.022 (0.990) | 1.027 (1.001) |
| 1000 | 2.5 cm | 0.57 | 1.002 (0.976) | 1.071 (1.042) |
| 1000 | 5 cm | 0.93 | 1.039 (1.008) | 1.076 (0.999) |
| 1000 | 10 cm | 1.00 | 1.018 (0.989) | 1.026 (1.004) |

Larger voxels widen the window and reduce the high bias of `density` (Section 8.1), and are observed more fully;
smaller voxels resolve finer structure but more of them are unobserved at high density (the 1.002 at 2.5 cm and
1000 stems/m² is the cancellation of Section 8.3).

**Beam.** With the footprint correction disabled, 300 stems/m² read 2.30 (2.49 also without the occlusion
correction) and 100 stems/m² 2.47 (2.56). With a near-perfect beam (1 mm, 0.1 mrad), `density` read 1.031 (300) and
1.030 (1000), `pooled_density` 0.991 and 0.992.

### 7.7 Tiling and sloped ground

**Tiling.** A 300 stems/m² scene estimated in 3 × 3 tiles of 2 m with `scripts/raydensity_tiles.sh` gave cells
identical to a single 6 m tile (`density`, `pooled_density`, `hits`, `weighted_path`, `observed_fraction`,
`patchiness_ratio`; one scene). The test suite (`RayDensityTiles`) also checks agreement.

**Sloped ground** (300 stems/m², each scene sheared by a ground plane of slope 10% in x and −5% in y, 4 scenes):

| | 0–0.1 m | 0.1–0.2 m | 0.2–0.3 m |
|---|---|---|---|
| `raydensity --ground` / flat, `density` (largest cell difference) | 0.999 (2.4%) | 0.995 (2.6%) | 0.997 (3.5%) |
| `raydensity --ground` / flat, `pooled_density` | 1.000 | 0.990 | 0.984 |
| Absolute heights / flat, `density` | 0.763 | 0.716 | 0.952 |
| Flattened with `raytranslate subtract` / flat, `density` | 0.429 | 0.583 | 0.712 |

Only this slope has been tested.

### 7.8 Trunks and their separability

Scenario: 300 stems/m² with trunks at 0.1 per m² (mean diameter 0.2 m), so a cell holding a trunk has 0.18–0.25 m⁻¹
of trunk frontal area, about 6–8% of its pneumatophore $`N d`$. 1 × 1 m cells, 4 scenes (2–5 trunks in each plot),
`density` / pneumatophore truth (s.d. 0.2–3% in trunk cells, under 1% elsewhere):

| | Cells with a trunk, 0–0.1 / 0.1–0.2 / 0.2–0.3 m | Cells without |
|---|---|---|
| A: trunks not separated | 1.15 / 1.17 / 1.23 (1.06 / 1.07 / 1.06 of pneumatophores + trunk) | 1.07 / 1.08 / 1.12 |
| B: trunk returns removed (perfect labels) | 1.10 / 1.12 / 1.14 | 1.07 / 1.08 / 1.12 |
| C: trunks detected above the pneumatophores and removed | 1.10 / 1.11 / 1.14 | 1.07 / 1.08 / 1.12 |

**Impact.** Unseparated, a trunk cell's estimate is close to the combined pneumatophore and trunk frontal area, but
this is a coincidence: alone (scenes with almost no pneumatophores), trunks read 0.70 ± 0.12 of their frontal area
with the stem settings, 0.82 ± 0.15 with the trunk's own diameter, and 1.76 ± 0.31 without the footprint and
occlusion corrections, varying widely between scenes. A 20 cm trunk spans several voxels, so violates the
stem-scale assumptions (thin elements, footprint widening of d, the occlusion correction for d). Trunks should be
separated, and their frontal area added from their fitted diameters. A trunk also shadows the volume behind it.

**Separability.** In simulation, trunks were cleanly separable in the bottom 0.3 m. In each scene, every trunk within
4 m of the plot centre (5–8 per scene) was found by clustering returns in a band above the tallest pneumatophores
(about 0.6–1.3 m) and fitting a horizontal circle to each cluster: centres within 0.2 cm and diameters 0–0.4 cm too
large. Removing returns within each fitted circle plus 1.5 cm, at all heights, removed 99.8–100% of trunk returns and
0.03–0.06% of pneumatophore returns in the bottom 0.3 m. The pneumatophore estimate in trunk cells then exceeded
trunk-free cells by 2–3%, consistent with the trunk's 1–8% share of the cell's ground area: the estimate is per unit
of free (observable) volume, so it should be multiplied by (1 − trunk cross-section / cell area) for a per-cell
value. These results depend on noise-free ranges: real range noise and SLAM pose error (Section 9.1) blur the trunk
outlines, and the removal margin must exceed them.

**Caveats for real stands.** The simulation used perfectly vertical, circular trunks of uniform cross-section and
nothing else above the pneumatophores. Separation will be harder where trunks lean, branch low, flare at the base,
or are multi-stemmed; where saplings or seedlings (1–5 cm) fall between the pneumatophore and trunk scales; and where
low branches, prop roots, debris or epiphytes return in the detection band. Checking the fitted circle at several
heights would detect lean.

### 7.9 Wider validation matrix (stage 2)

The scenes above share one scan geometry, perfect ranges and poses, one stem diameter and one height distribution.
The validation matrix (`scripts/pneumatophore_validation_matrix.py`) varies these one block at a time around a
reference scene, and asks whether `density` meets a ±2% criterion with and without `--self_shadow_correction`.

**Design.** Reference: the stage-2 scenes of Section 7.1 (4 m plot, 56 250 rays per metre of walk at 1 m lines,
1 cm non-overlapping stems of 0.25 ± 0.08 m at 300 stems/m², the default beam with `surface` returns), estimated with
`--estimator pooled --no_flat_top`, 5 cm voxels and the default pool radius, without ("default") and with ("ss")
`--self_shadow_correction`.

| Block | Factor | Levels | Stems/m² | Configurations |
|---|---|---|---|---|
| A | Scan geometry | walk lines 0.5, 1, 2, 3 m apart; scanner at 0.8, 1.2, 2 m; zenith bands 0–101.5°, 30–90°, 60–95°; an 8 m plot; a single walk line across an 8 m plot | 300, 1000 (100 for the single line) | 20 |
| B | Range noise | $`\sigma`$ = 0, 0.5, 1, 2, 3 cm, with `surface` and `outline` returns (estimator matched) | 300, 1000 | 20 |
| C | Pose noise | $`\sigma`$ = 1, 2, 5 cm with knots every 0.5 or 2 m, and none | 300 | 7 |
| D | Diameter | true 0.5, 1, 2 cm, assumed 0.7–1.3 times it (pass tested at 1.0 only); a 6/14 mm mixture and log-normal 10 ± 3 mm, assumed 10 mm | 300, 1000 | 8 |
| E | Stem heights | normal 0.25 ± 0.08 and 0.4 ± 0.15 m; uniform 0.05–0.5 m (reported only); one height 0.6 m; 2 m (no tops in view) | 300, 1000 | 10 |
| F | Low observation (reported only) | stems of 0.25 ± 0.08 m or one height 0.6 m, lines 1 or 2 m apart | 600, 1000, 1500 | 12 |
| G | Clustering (reported only) | cable roots with lateral spread 1.5, 3, 6 cm, pool radius 0.05 and 0.1 m; patches of radius 0.2, 0.3, 0.5 m, cover 0.3 or 0.6, stems between patches 0 or 0.2 of the density within | 150, 500; 300 | 12 + 12 |
| H | Return convention | generator × estimator `surface` / `outline` (pass tested where they match); heights normal or 0.6 m | 300, 1000 | 8 scenes |
| I | Spacing | Poisson (overlapping) stems; minimum gap 0, 0.5, 1 diameters (estimator given the true gap) | 300, 1000 | 8 |
| J | Window | (voxel, pool radius) = (0.025, 0.05), (0.05, 0.05), (0.05, 0.1), (0.05, 0.15), (0.1, 0.1) m, over the 25 distinct 300 stems/m² scenes of A–C | 300 | 125 |
| X | Diagnostics of the failures (not part of the pass count) | stems of one height combined with outline returns, noise, 2 cm stems and the mixture; wider windows for 2 cm stems; tops on and off a layer boundary; no footprint | 300, 1000 | 13 |

**Metric and decision rule.** d/t90 is Σ `density` / Σ truth over the cells of a 0.1 m layer whose
`observed_fraction` is at least 0.9, as the mean ± s.e. over seeds 101–104. A configuration passes if
|d/t90 − 1| ≤ 0.02 in both 0–0.1 and 0.1–0.2 m; a layer with no cell observed to 0.9 is untested, and a
configuration with neither layer tested is untestable. The 37 configurations whose ss estimate failed with 4 seeds
were re-run with seeds 105–108 and judged on all 8. Blocks F and G, the uniform 0.05–0.5 m heights, the assumed
diameters other than the true one and the crossed conventions were reported, not judged: stage 1 had already shown
biases there that self-shadowing does not cause. The harness also reports `pooled_density` (p/t90), cell RMSE, the
observed fraction, and `density` / truth in bins of observed fraction.

**Table 11.** Pass counts and the range of d/t90 over the configurations judged.

| Block | Judged | Default: pass / fail | Default 0–0.1 / 0.1–0.2 m | ss: pass / fail | ss 0–0.1 m | ss 0.1–0.2 m |
|---|---|---|---|---|---|---|
| A | 20 | 0 / 20 | 1.037–1.073 / 1.061–1.091 | 20 / 0 | 0.983–1.006 | 0.997–1.020 |
| B | 20 | 0 / 20 | 1.029–1.084 / 1.043–1.123 | 9 / 11 | 0.973–1.016 | 0.984–1.061 |
| C | 7 | 0 / 7 | 1.025–1.068 / 0.951–1.078 | 5 / 2 | 0.961–1.001 | 0.890–1.010 |
| D | 8 | 0 / 8 | 1.040–1.078 / 1.038–1.114 | 6 / 2 | 0.984–1.008 | 0.983–1.027 |
| E | 5 (3 untestable) | 0 / 5 | 1.042–1.076 / 1.055–1.080 | 5 / 0 | 0.985–1.011 | 0.998–1.020 |
| H | 6 (2 untestable) | 0 / 6 | 1.040–1.084 / 1.068–1.123 | 4 / 2 | 0.984–1.019 | 1.004–1.061 |
| I | 8 | 0 / 8 | 0.925–1.067 / 0.970–1.078 | 5 / 3 | 0.874–1.001 | 0.915–1.018 |
| J | 125 | 0 / 125 | 0.993–1.089 / 0.916–1.129 | 90 / 35 | 0.951–1.016 | 0.885–1.045 |
| All | 199 | 0 / 199 | | 144 / 55 | | |

The default fails everywhere, mostly by its self-shadowing bias. With ss, every scan geometry passed, including lines
3 m apart, a 2 m scanner height, shallow rays only and a single walk line; so did range noise up to 2 cm with
`surface` returns, pose noise up to 2 cm, 0.5 and 1 cm stems, log-normal diameters, every height distribution
judged, a minimum gap with the true gap given (but for a near-miss at 1000 stems/m²), and every window but one for the
scenes that passed at the default window. Where it was testable, low observation (block F) read 0.985–1.008 at
0–0.1 m and 0.970–1.021 at 0.1–0.2 m with ss; no cell of the lowest 0.1 m was 0.9 observed at 600–1500 stems/m² with stems of one height, at 1000–1500 stems/m² with 2 m lines, or at 1500 stems/m².

**Table 12.** The configurations that failed with ss (d/t90, 0–0.1 / 0.1–0.2 m, 8 seeds) and their causes.

| Configuration | ss | Diagnosed cause |
|---|---|---|
| `outline` returns, 300 stems/m², range noise 0–3 cm (B, H) | 1.009–1.016 / 1.029–1.044 | stem tops; +1.2–1.9% where no stems end |
| `outline` returns, 1000 stems/m², range noise 0–3 cm (B, H) | 0.991–0.999 / 1.030–1.061 | the same |
| Range noise 3 cm, `surface`, 1000 stems/m² (B) | 0.973 / 0.984 | a near-miss at 1000 stems/m² |
| Pose noise 5 cm, knots every 0.5 / 2 m (C) | 0.971 / 0.926; 0.961 / 0.890 | trajectory error where stems end |
| 2 cm stems, 300 stems/m² (D) | 0.995 / 1.027 | the tops of thick stems |
| Mixture of 6 and 14 mm stems, assumed 10 mm (D) | 1.008 / 1.022 | the mixture itself (+0.7–1.3%) |
| Poisson (overlapping) stems, 300 / 1000 stems/m² (I) | 0.968 / 0.982; 0.874 / 0.915 | the occlusion correction assumes stems that don't overlap |
| Minimum gap 0.5, 1000 stems/m² (I) | 0.976 / 1.007 | a near-miss at 1000 stems/m² |
| Window (J): the outline and 5 cm pose scenes at every window (34), and shallow rays (60–95°) at 10 cm voxels (1) | 0.951–1.016 / 0.885–1.045 | as for their scenes; 0.996 / 0.976 for the shallow rays |

The diagnoses (block X and the layers above 0.2 m):

* **Stem-top discs.** Stem tops are not in the measurement model (Section 8.2), and the generator makes them solid
  discs of the widened outline. In one reference scene, the share of stem returns that came from top discs was 0.3%,
  2.5%, 6.7%, 12% and 15% in the 0.1 m layers from the ground up to 0.5 m, and ss read 1.001, 1.010, 1.052, 1.117 and
  1.203 in the same layers (block E, 4 seeds). With stems of one height (0.6 m, no tops below 0.5 m) ss read
  1.000–1.007 in every layer from 0 to 0.5 m; without a beam footprint the excess remained (1.001, 1.012, 1.048,
  1.101, 1.172). So the excess above 0.1 m is mostly the stem tops, and the 1% at 0.1–0.2 m with the default heights
  is their first appearance. The tops of 2 cm stems explain that failure: with one height they read 0.989 / 0.982,
  while wider windows, which reduce the self-shadowing, left 0.1–0.2 m at 1.032 (`--pool_radius 0.1`) and 1.035
  (0.15).
* **Outline returns.** With returns ranged to the widened outline, ss read 1.2–1.9% high at 0–0.4 m even for stems
  of one height (1.012–1.019, 300 stems/m², with and without 3 cm range noise), and 3–4% more at 0.1–0.2 m where
  stems of the reference heights end. The
  estimator's `--return_range` switch cannot remove this: estimating the same scenes with the other convention
  changed every value by 0.07–0.14% in the layers judged (up to 0.4% in the poorly observed 1000 stems/m² scenes of
  one height), because the response tables of the two conventions differ by under 0.2% at these solid fractions
  (Section 5). Which convention a real scanner follows has to be measured (Section 9.2, step 1).
* **Pose error.** 5 cm of trajectory error lowered both outputs alike (`pooled_density` 0.971 / 0.929 and
  0.958 / 0.890), so it is not a self-shadowing effect; with stems of one height it disappeared below 0.2 m
  (0.997 / 1.000) and left 0.981–0.994 at 0.2–0.5 m. It acts where stems end, smearing the tops and the hits near
  them vertically, and 2 cm or less was harmless.
* **Overlapping stems.** The occlusion correction assumes stems that don't overlap (Section 5). Multiplying the Poisson
  results by the response at their density, $`R = 1.046`$ at 300 and 1.157 at 1000 stems/m², gives 1.013 and 1.011 at
  0–0.1 m, so most of the deficit is that model mismatch, not the self-shadow correction. Real stems don't overlap.
* **Diameter mixture.** An equal mixture of 6 and 14 mm stems, estimated with the 10 mm mean, read 0.7–1.3% higher
  than stems of one diameter, also with stems of one height (1.020 / 1.017), so it is a property of the mixture:
  the footprint widens thin and thick stems by different amounts, so the mean effective width is not that of the
  mean diameter. Log-normal diameters with a 30% coefficient of variation added only 0.3–0.7%.
* **Near-misses at 1000 stems/m².** At 1000 stems/m² ss read 0.973–0.994 at 0–0.1 m and up to 1.020 at 0.1–0.2 m
  in blocks A, B, E and I, so a small further error tips a configuration over the line (3 cm range noise; a minimum
  gap of 0.5 with the true gap given). In the reference scene at this density the lowest layer is 0.92 observed and
  82% of its cells reach 0.9, so the cells judged are partly a selection (Section 8.3).
* **Layer-dependent error of the pooled value.** ss equals the pooled value times the ratio `density` /
  `pooled_density`, and `pooled_density` itself has an error that changes with the layer and the density:
  0.981 / 0.992 (0–0.1 / 0.1–0.2 m) at 300 stems/m² and 1.012 / 1.007 at 1000 stems/m² for the reference heights,
  and 1.027, 1.006, 0.987, 0.979 and 0.976 in the 0.1 m layers from the ground for stems of one height (300 stems/m²),
  where `density` with ss was flat (1.000–1.007). So `density` / `pooled_density` is not a clean measure of the
  self-shadowing in the lowest layer, and the 2% that the first-order correction leaves at 300 stems/m² (Section 8.1)
  is cancelled there by the pooled value's −1.9%. In single-voxel windows, `--occlusion_radius` removed the −1.9% at
  300 stems/m², but not the +1% at 1000 (Section 8.1.2); its cause at the default window is not established.
* **Grid-boundary artefact.** The layer just below the tops of stems of one height read 0.918 with ss (0.5–0.6 m,
  stems of 0.6 m), whether the grid ended at 0.6 m or extended to 0.7 m, and 1.047 when the tops lay inside the layer
  (stems of 0.55 m). So it is not an edge effect of the grid's top, but of tops that lie exactly on a voxel and layer
  boundary; the likely cause is that their top-disc returns fall in the voxel above, leaving the layer with neither
  top returns nor returns carried down from above (Section 8.2). Real stems vary in height, so this is an artefact
  of the uniform-height scenes; their layer holding the tops was left out of every result above.

**Observed fraction.** Pooled over the uniform-stem configurations of blocks A and F, cells at least 0.9 observed
read 0.995 / 1.010 with ss (0–0.1 / 0.1–0.2 m; 1683 / 2187 cells), but below 0.9 they read low in every bin:
0.971 / 0.965 at 0.8–0.9, 0.960 / 0.954 at 0.7–0.8, 0.968 / 0.944 at 0.5–0.7, 0.894 / 0.897 at 0.3–0.5 and
0.780 / 0.831 below 0.3 (default: 1.060 / 1.076 at 0.9–1, and 1.033 / 1.027 to 0.824 / 0.873 below). This confirms
the 0.9 threshold of Section 8.3 over the wider range of scans.

### 7.10 Beam footprint mismatch (block L)

Every scene so far was estimated with the footprint it was generated with. Block L of the validation matrix
generates scenes with a "true" beam and estimates them assuming another, to find how accurately $`w_0`$, $`\gamma`$ and $`\tau`$ must
be known and how to calibrate them (Section 9.2).

**Design.** Three true beams: default (10 mm, 3 mrad, $`\tau = 0.05`$), narrow (5 mm, 0.5 mrad, 0.05) and wide (20 mm,
5 mrad, 0.05), at 300 and 1000 stems/m² (reference scene of Section 7.9 otherwise, 4 seeds). Each scene was
estimated with the base (default), ss and sw estimators, assuming:

* the exit diameter or the divergence at 0.5, 0.75, 1.25, 1.5 and 2 times the true value (e0.5 … e2, g0.5 … g2);
* thresholds of 0.02, 0.03, 0.08, 0.1 and 0.15 (t0.02 … t0.15);
* joint errors that give the right beam diameter or $`d_\mathrm{eff}`$ at 2 m but not at other ranges: twice the divergence with
  a smaller exit (j_steep), no divergence with a larger exit (j_flat), and the exit at 0.5, 1.5 and 2 times with
  the threshold that restores $`d_\mathrm{eff}`$ at 2 m (jt_e0.5, jt_e1.5, jt_e2);
* no footprint correction (nofp);
* for the default beam, an assumed diameter of 0.8–1.2 times the true 1 cm crossed with the thresholds 0.02–0.15.

A further scene has one walk line along the middle of an 8 m plot (default beam, 300 stems/m²), and is binned by the
cells' distance from the line, to see the range dependence. The metric is d/t90 (Section 7.9) and, for the
sensitivities, the ratio r90 of each estimate to the same estimator with the true footprint, paired by seed.

**Table 15.** d/t90 at 0–0.1 m by assumed footprint, and the ratio $`Q`$ of the assumed to the true path-weighted mean
$`d_\mathrm{eff}`$ over the same cells (default beam, 300 stems/m²). The joint rows of the narrow and wide beams are their own
2 m matches (e.g. j_steep is 4 mm, 1 mrad for the narrow beam and 10 mm, 10 mrad for the wide one). s.e. ≤ 0.005.

| Assumed | Default beam: $`Q`$ | base | ss | sw | ss, 1000/m² | ss, narrow beam | ss, wide beam |
|---|---|---|---|---|---|---|---|
| true | 1.000 | 1.067 | 1.001 | 0.999 | 0.985 | 0.998 | 1.003 |
| e0.5 | 0.814 | 1.302 | 1.236 | 1.276 | 1.207 | 1.162 | 1.352 |
| e0.75 | 0.907 | 1.173 | 1.107 | 1.123 | 1.085 | 1.074 | 1.152 |
| e1.25 | 1.093 | 0.979 | 0.913 | 0.896 | 0.902 | 0.932 | 0.890 |
| e1.5 | 1.185 | 0.905 | 0.839 | 0.810 | 0.832 | 0.874 | 0.800 |
| e2 | 1.369 | 0.787 | 0.721 | 0.674 | 0.719 | 0.776 | 0.669 |
| g0.5 | 0.912 | 1.167 | 1.100 | 1.115 | 1.060 | 1.023 | 1.107 |
| g0.75 | 0.956 | 1.115 | 1.048 | 1.054 | 1.021 | 1.011 | 1.053 |
| g1.25 | 1.044 | 1.024 | 0.957 | 0.948 | 0.951 | 0.986 | 0.959 |
| g1.5 | 1.088 | 0.984 | 0.918 | 0.902 | 0.920 | 0.975 | 0.918 |
| g2 | 1.176 | 0.913 | 0.847 | 0.819 | 0.863 | 0.952 | 0.846 |
| t0.02 | 1.136 | 0.943 | 0.877 | 0.854 | 0.874 | 0.922 | 0.847 |
| t0.03 | 1.079 | 0.992 | 0.926 | 0.911 | 0.918 | 0.953 | 0.907 |
| t0.08 | 0.920 | 1.157 | 1.091 | 1.104 | 1.064 | 1.049 | 1.124 |
| t0.1 | 0.879 | 1.209 | 1.143 | 1.166 | 1.110 | 1.077 | 1.199 |
| t0.15 | 0.797 | 1.329 | 1.262 | 1.306 | 1.214 | 1.138 | 1.380 |
| j_steep (4 mm, 6 mrad) | 0.953 | 1.119 | 1.052 | 1.060 | 1.070 | 1.006 | 1.079 |
| j_flat (16 mm, 0 mrad) | 1.047 | 1.020 | 0.954 | 0.943 | 0.913 | 0.991 | 0.938 |
| jt_e0.5 (5 mm, $`\tau`$ 0.0084) | 0.979 | 1.090 | 1.024 | 1.026 | 1.022 | 1.004 | 1.043 |
| jt_e1.5 (15 mm, $`\tau`$ 0.104) | 1.012 | 1.055 | 0.989 | 0.984 | 0.966 | 0.996 | 0.980 |
| jt_e2 (20 mm, $`\tau`$ 0.150) | 1.020 | 1.046 | 0.980 | 0.974 | 0.952 | 0.995 | 0.962 |
| no footprint | 0.452 | 2.302 | 2.237 | 2.457 | 2.034 | 1.491 | 3.275 |

**Sensitivities.** As $`\mathrm{d}\ln(\mathrm{density}) / \mathrm{d}\ln(\mathrm{parameter})`$, with ss, the default beam at 300 stems/m² and 0–0.1 m
(local values from the ×0.75 / ×1.25 levels and the thresholds 0.03 / 0.08; least-squares fits over all levels
agree within 0.02):

| Parameter | Default beam | Narrow beam | Wide beam |
|---|---|---|---|
| Exit diameter $`w_0`$ | −0.38 | −0.28 | −0.51 |
| Divergence $`\gamma`$ | −0.18 | −0.05 | −0.18 |
| Threshold $`\tau`$ | +0.17 | +0.10 | +0.22 |
| Stem diameter $`d`$ (assumed) | +0.46 | | |

sw is about 1.2 times as sensitive ($`w_0`$ −0.44, $`\gamma`$ −0.21, $`\tau`$ +0.20) and the base estimator about 0.94 times. The
response is not symmetric: relative to the true-footprint estimate, e0.75 / e1.25 gave +10.6% / −8.8%, and t0.03 /
t0.08 −7.5% / +9.0%. The ranges of each parameter alone that keep the density within ±2% and ±5% are in
Section 9.2 (Table 18).

**One curve: the mean effective width.** Whatever combination of $`w_0`$, $`\gamma`$ and $`\tau`$ is wrong, the density error follows
the error $`Q`$ in the path-weighted mean $`d_\mathrm{eff}`$ (equivalently $`\bar\rho`$) at the ranges actually scanned.
Fitting $`\ln r_{90} = b_1 \ln Q + b_2 (\ln Q)^2`$ to all the single-parameter errors of each true beam, density, estimator and layer gave $`b_1`$ of
−0.96 to −1.05 for ss, −1.10 to −1.31 for sw and −0.87 to −0.98 for base, with RMS residuals of 0.01–0.15%. The
joint errors, not used in the fit, fall on the same curves within 0.74% (ss, default beam: within 0.23%; e.g.
j_flat at 1000 stems/m², $`Q`$ 1.082, r90 0.927 against 0.927 predicted). So, to a good approximation,

```math
\frac{\Delta a}{a} \approx -1.0\,\frac{\Delta \bar d_\mathrm{eff}}{\bar d_\mathrm{eff}}\ \ (\text{ss}), \qquad -1.2\times\ (\text{sw}), \qquad -0.95\times\ (\text{base}), \qquad (15)
```

and a beam calibration needs to fix only the mean $`d_\mathrm{eff} / d`$ over the survey's ranges, not each parameter. The
exception is the stem diameter: an assumed diameter also enters the occlusion and self-shadow terms, so the
diameter × threshold grid lies off the curve by up to 2.4–4.2% with ss (1.1–3.3% base, 7.6–9.1% sw).

**Diameter × threshold.** A wrong threshold can absorb a wrong diameter. With ss at 300 stems/m² (0–0.1 m):

| Assumed $`d`$ / true $`d`$ | $`\tau`$ 0.02 | $`\tau`$ 0.03 | $`\tau`$ 0.05 | $`\tau`$ 0.08 | $`\tau`$ 0.1 | $`\tau`$ 0.15 | $`\tau^*`$ (d/t90 = 1) |
|---|---|---|---|---|---|---|---|
| 0.8 | 0.776 | 0.823 | 0.897 | 0.986 | 1.039 | 1.163 | 0.085 |
| 0.9 | 0.829 | 0.877 | 0.952 | 1.042 | 1.095 | 1.217 | 0.064 |
| 1.0 | 0.877 | 0.926 | 1.001 | 1.091 | 1.143 | 1.262 | 0.050 |
| 1.1 | 0.919 | 0.968 | 1.044 | 1.133 | 1.184 | 1.301 | 0.037 |
| 1.2 | 0.957 | 1.007 | 1.082 | 1.169 | 1.220 | 1.335 | 0.028 |

A counted-plot calibration would pick $`\tau^*`$ for the diameter assumed. Within the same stand $`\tau^*`$ transfers: with ss, the
0.1–0.2 m layer then reads 1.009–1.010 for every diameter error, the model's assumed / true $`\bar\rho`$ is the same at 1, 2 and
4 m to within 0.4%, and for 2 cm stems it differs from that for 1 cm stems by at most 0.7%. But at 1000 stems/m² the
same $`\tau^*`$ gives 0.966–1.004 at 0–0.1 m (a 4% spread over the diameter errors), and with the base estimator $`\tau^*`$ also
absorbs its +7% self-shadowing bias ($`\tau^* = 0.032`$ at the true diameter), which then transfers less well (the $`\bar\rho`$ ratio
for 2 cm stems differs by up to 2.5%). $`\tau^*`$ is therefore an empirical fit, not the physical threshold, and must not be
reused for another scanner, stem size or density.

**Range.** On the single-line plot (ss, 0–0.1 m), r90 by the cells' distance from the walk line:

| Assumed | 0.5 m | 1.5 m | 2.5 m | 3.5 m |
|---|---|---|---|---|
| e0.75 / e1.25 | 1.106 / 0.912 | 1.099 / 0.917 | 1.090 / 0.923 | 1.082 / 0.929 |
| g0.75 / g1.25 | 1.046 / 0.958 | 1.059 / 0.947 | 1.075 / 0.934 | 1.091 / 0.923 |
| t0.03 / t0.08 | 0.925 / 1.089 | 0.922 / 1.093 | 0.918 / 1.098 | 0.915 / 1.104 |
| j_steep / j_flat | 1.061 / 0.946 | 0.994 / 1.007 | 0.925 / 1.089 | 0.867 / 1.181 |
| jt_e0.5 / jt_e2 | 1.027 / 0.975 | 0.997 / 1.003 | 0.964 / 1.039 | 0.933 / 1.079 |
| Observed fraction (true footprint) | 1.00 | 0.83 | 0.30 | 0.05 |

Exit and threshold errors are nearly flat with range (1–2.5% change over 0.5–3.5 m for these levels, 4–5% for e0.5
and t0.15), while divergence errors grow with it. Pairs matched only at 2 m are right near 1.5–2 m and wrong
elsewhere, by −13% to +18% at 3.5 m. In the survey geometry of the matrix (lines 1 m apart), such a pair still
gives errors of 1.2–5.2% at 300 stems/m² and 1.9–8.6% at 1000 (ss, default beam; up to 7.6% and 10.7% for the wide
beam). A calibration at a single range cannot detect this (Section 9.2).

## 8. Known biases and limitations

### 8.1 Self-shadowing bias of `density`, and an optional correction

`density` reads high by 4–8% for random stems at the default settings, by the same amount at 30 and 600 stems/m²
(Table 2), and whatever the number of rays (Section 7.4). `pooled_density`, the path-weighted mean, does not have
this bias for uniform stems.

**Mechanism.** A stem removes the free path behind it within its own pooling window, but keeps its own hits.
Windows that hold a stem therefore have less path than their true share and a higher hits / path, and the
equal-weight mean of the window ratios is high. The evidence (stage-1 study, Section 7.1):

* *Own-shadow test.* A 2D Monte Carlo of the window estimator (`scripts/pneumatophore_self_shadow_mc2d.cpp`:
  infinitely tall stems in a periodic domain, lines recorded only after a pre-path long enough that every stem
  upstream of a recorded point lies on the line's past, as for tall stems in 3D) can continue each line through the
  stem it returned from, within that stem's own window only. This removes the bias: at the default window
  ($`W = 0.15`$ m) and $`d_\mathrm{eff}`$ of 1–2.3 cm, the window mean relative to the path-weighted value fell from 1.03–1.09 to
  0.99–1.02 (at 600 stems/m², from 1.10 to 1.04).
* *Crediting the whole path* behind every hit to the window edge (the crossing estimator's idea at window scale)
  overcorrects: 0.48–0.96 of the truth in the same model.
* *Not ray count.* 0.5, 2 and 8 M rays gave the same `density` (Section 7.4), and a quarter of the rays changed it by
  under 2% at the default window (30 and 300 stems/m², stage-1 scenes).
* *Not stem tops or the ground.* Stems of one height (0.6 m) and stems taller than the scanner (2 m, no tops in
  view) read 1.060–1.075 at 0.1–0.4 m where fully observed (100 and 300 stems/m²), as the standard scenes do
  (1.067–1.087 at 0–0.2 m).
* *Smaller in very regular stands.* In the 2D model, stems at least one diameter apart gave +0.9% against +2.9% for
  random stems (1000 stems/m², $`d_\mathrm{eff} = 1`$ cm). A reduction predicted from the hard-core spacing of stems with no gap
  was not detected.

**First-order bias.** A hit's own stem shadows the rest of that ray's path to the edge of the window, on average the
mean exit distance $`s`$ of the window (the mean distance from a uniformly random point of the window to its edge along a
uniformly random horizontal direction). To first order (derived for tall stems and horizontal ray directions, not
fitted), the relative bias of the window mean is

```math
\beta = \frac{s\,d_\mathrm{eff}}{W^2} = c\,\frac{d_\mathrm{eff}}{W}, \qquad c = \frac{4\ln(1+\sqrt2) - \tfrac43(\sqrt2-1)}{2\pi} = 0.4732 \qquad (11)
```

for a square window of width $`W`$ ($`s = c W`$). Windows clipped at the grid's edge are rectangles: $`W^2`$ is replaced by their
area and $`s`$ by the rectangle's mean exit distance (`ray::meanExitDistance` in `raylib/raystemmodel.cpp`). To first
order the bias does not depend on the stem density or the number of rays.

**Table 8.** First-order bias against the measured bias, measured as `density` / `pooled_density` − 1 without the
correction (for uniform stems `pooled_density` has no self-shadowing bias, so the ratio leaves out the occlusion
correction's errors), random stems, lowest 0.1 m, well observed cells, stage-1 scenes.

| Case | $`W`$ (m) | $`d_\mathrm{eff}`$ (cm) | $`d_\mathrm{eff}/W`$ | Eq. (11) | Measured |
|---|---|---|---|---|---|
| `--pool_radius 0.15` | 0.35 | 2.2–2.3 | 0.064 | +3.0% | +2.7 to +2.9% |
| No footprint ($`d_\mathrm{eff} = d`$) | 0.15 | 1.0 | 0.067 | +3.2% | +3.5 to +3.7% (300, 1000 stems/m²) |
| `--pool_radius 0.1` | 0.25 | 2.2–2.3 | 0.089 | +4.2% | +4.4 to +4.7% |
| Default | 0.15 | 2.2–2.3 | 0.149 | +7.0% | +8.8% (300 stems/m²; 2D model also +8.8%) |
| Wide footprint (20 mm, 5 mrad) | 0.15 | 3.3 | 0.22 | +10.6% | +12 to +14% |
| Single voxel (`--pool_radius 0`) | 0.05 | 2.2–2.3 | 0.45 | +21% | +28 to +29% |

Eq. (11) is within 10% of the measured bias up to $`d_\mathrm{eff}/W \approx 0.1`$. Beyond that the bias is about 1.15–1.35 times
larger: the second-order term, derived below. (Over 30–600 stems/m² the default case measured +8.7 to +9.5% in both
layers, except +6.1% at 600 stems/m² and 0–0.1 m.) Consequences:

* More scanning will not reduce the bias, but it is predictable: it scales with $`d_\mathrm{eff}/W`$, so it is larger for thicker
  stems, longer ranges ($`d_\mathrm{eff}`$ of 3 cm at 5 m gives about +10% at $`W = 0.15`$ m) and smaller windows, and the +4–8% holds
  only for this configuration.
* The `patchiness_ratio` of about 1.09 for random stems is mostly this bias.
* A larger window (`--pool_radius 0.1`, or larger voxels) reduces it, but reads low where stems are clustered or
  patchy at a smaller scale than the window (Table 3). The choice of pool radius therefore trades this bias against
  the clustering bias; it is not a remedy for either.

**The correction (`--self_shadow_correction`, off by default).** With this option each voxel's neighbourhood density
of Eq. (8) is divided by its first-order factor before the occlusion correction of Eq. (9):

```math
m_v' = \frac{m_v}{1 + s_v\,\bar\rho_v\,d / A_v} \qquad (12)
```

where $`A_v`$ is the area of the neighbourhood ($`W^2`$ away from the grid's edge), $`s_v`$ its mean exit distance and
$`\bar\rho_v d`$ its mean effective stem width. `pooled_density` is left as it was. The option needs `--estimator pooled` and
`--stem_diameter`, and `rayrender` and `raydensity` refuse it otherwise. It is deterministic, added about 9% to
the run time, and the run prints the mean
factor (1.071 for a 300 stems/m² scene at the default settings). `raydensity` records it in the store's `.zattrs`
(`"self_shadow_correction": true`, written only when set).

**Table 9.** `density` / truth without and with the correction (stage-1 scenes, 4 seeds; cases whose layer is at
least 95% observed, except where stated).

| Case | Without, 0–0.1 m | With, 0–0.1 m | Without, 0.1–0.2 m | With, 0.1–0.2 m |
|---|---|---|---|---|
| Random, 30–600 stems/m² | 1.067–1.075 | 1.000–1.008 | 1.078–1.087 | 1.010–1.021 |
| No footprint, 300–1500 stems/m² | 1.019–1.031 | 0.992–1.001 | 1.039–1.044 | 1.010–1.016 |
| Wide footprint, 100–300 stems/m² | 1.100–1.110 | 1.003–1.006 | 1.121–1.126 | 1.020–1.022 |
| Stem diameter 0.5 and 2 cm | 1.057–1.100 | 1.000–1.003 | 1.057–1.128 | 1.008–1.028 |
| Voxels 2.5–10 cm, pool radius 0.05–0.15 m ($`W`$ 0.125–0.35 m) | 1.023–1.089 | 0.993–1.004 | 1.031–1.102 | 0.997–1.023 |
| One height (0.6 m) or taller than the scanner (2 m), 100–300 stems/m², 0.1–0.4 m | | | 1.060–1.075 | 0.994–1.007 |
| 1 cm range noise or 2 cm pose noise, 300 stems/m² | 1.061–1.066 | 0.995–1.000 | 1.058–1.076 | 0.991–1.008 |
| Overlapping (Poisson) stems, 300 stems/m² | 1.032 | 0.968 | 1.047 | 0.981 |
| Patches of 0.2–0.5 m radius, 30% cover, 300 stems/m² | 0.986–1.022 | 0.922–0.957 | 1.021–1.047 | 0.955–0.980 |
| Patches of 0.3 m radius, 30% cover, 1000 stems/m² | 0.834 | 0.786 | 0.931 | 0.877 |
| Cable roots, 150 stems/m² | 0.867 | 0.811 | 0.898 | 0.840 |
| Cable roots, 150 stems/m², `--pool_radius 0.1` | 0.757 | 0.726 | 0.789 | 0.757 |
| Cable roots, 500 stems/m² | 0.470 | 0.441 | 0.587 | 0.552 |

For uniform stems the correction leaves `density` within 1% of the truth in the lowest 0.1 m and 1–3% high at
0.1–0.2 m, the residual being mostly the second-order term (Section 8.1.1). Overlapping stems read low with it because the occlusion
correction assumes non-overlapping stems (`pooled_density` read 0.945 and 0.957 there). The correction does not
change the sensitivity to the assumed diameter: ±20% gave −10.4% / +8.1% with it, against −10.8% / +8.6% without
(300 stems/m², 0–0.1 m).

*Caveat at 0–0.1 m.* The near-unity values in the lowest layer are partly a cancellation. The correction falls about
2% short (`density` / `pooled_density` was still 1.02 with it at 300 stems/m²), and in that layer `pooled_density`
itself reads about 2% low (0.981 at 300 stems/m²), for a reason not yet found (it varies with the layer and the
density: Section 7.9). At 0.1–0.2 m, where `pooled_density`
is within 1%, the corrected `density` reads 1–2% high.

**Why it is off by default.** In clumped stands the correction makes `density` worse: patches read 0.92–0.96 and cable
roots 0.81 with it (Table 9), and over the 24 clustered configurations of the validation matrix (block G, default
window) patches read 0.922–0.984 and cable roots at 150 stems/m² 0.759–0.921, depending on their lateral spread
(Section 8.6). The cause is heterogeneity within the window. A window that straddles a patch edge or a root line
under-weights its dense part (its ratio is path-weighted, and rays penetrate the dense part less), so it reads low.
Without the correction, the self-shadowing bias cancelled this. **The current accuracy of the default `density` in
patches (Table 3) and the size of the cable-root bias (−13%) therefore rest on that cancellation**, not on the
estimator being right for those arrangements, and may not hold for other patch sizes, densities, beams or windows.
For clumped stands the remedy found is a smaller window (Section 8.1.2), used per cell by the hybrid
(Section 8.1.3), not the correction at the default window.

**Remedies considered and rejected (at the default window).**

* Crediting the whole path behind each hit: overcorrects (0.48–0.96, above).
* Removing each stem's own shadow by stem identity: right in the model, but real scans don't identify the stems.
* Weighting windows by their expected unshadowed path: about −20% in patches (derived, not run).
* The geometric mean of `density` and `pooled_density`: 0.56–0.59 in patches and along cable roots.
* Switching the correction by `patchiness_ratio`: right only by the same cancellation, and discontinuous.

#### 8.1.1 The second-order factor

**Derived.** Take one isolated stem in a window of width $`W`$, crossed by many rays from all horizontal directions. Its
shadow removes a fraction $`\delta g`$ of the window's free path, where $`\delta = d_\mathrm{eff}/W`$ and $`g`$ is the mean exit distance from the
stem's position to the window's edge (averaged over directions, in units of $`W`$). The window's ratio is therefore
$`1/(1-\delta g)`$ times its share, and the mean over windows, with the stem uniformly placed, is

```math
1 + \beta = E\!\left[\frac{1}{1 - \delta g}\right] = 1 + E[g]\,\delta + E[g^2]\,\delta^2 + E[g^3]\,\delta^3 + \dots
= 1 + 0.4732\,\delta + 0.2274\,\delta^2 + 0.1108\,\delta^3 + \dots \qquad (13)
```

The moments are close to powers of the first ($`c^2 = 0.2239`$, $`c^3 = 0.1060`$), so the series is close to
$`1 + x + x^2 + x^3 + \dots = 1/(1-x)`$, with $`x = c\delta`$ the first-order term of Eq. (11): within 0.4% for $`x \le 0.3`$. The
self-shadow factor is therefore $`1/(1-x)`$ to second order, not $`1 + x`$. At $`\delta = 0.45`$ (single 5 cm voxels, $`d_\mathrm{eff}`$
2.25 cm) the expectation is 1.272 and $`1/(1-x)`$ is 1.271, against 1.213 to first order and 1.275 from the 2D Monte
Carlo; at the default window ($`\delta = 0.15`$) $`1/(1-x)`$ and $`1 + x`$ differ by 0.5%.

**Fitted.** Three further terms were fitted to the 2D Monte Carlo of the window estimator
(`scripts/pneumatophore_self_shadow_window_mc.cpp`; 146 runs over $`\delta`$ = 0.1–1, $`\Lambda = N d_\mathrm{eff} w`$ = 0.03–2 and
$`d/d_\mathrm{eff}`$ = 0.3, 0.44 and 1, two seeds each). The factor used is

```math
f = \frac{1}{\max\!\left(0.2,\ 1 - x - k_f\,\dfrac{\pi w_r^2}{16A} + 0.76\,x^3\right)} + 0.35\,x\,\tau, \qquad
x = \frac{s\,d_\mathrm{eff}}{A}, \quad \tau = N d_\mathrm{eff}\, s \qquad (14)
```

with $`A`$ the window's area and $`s`$ its mean exit distance (`ray::selfShadowFactor`, used with
`--self_shadow_second_order`):

* *Footprint spread*, $`k_f \pi w_r^2/(16A)`$: the reported returns spread along the ray over a width $`w_r`$, the physical
  stem for `surface` returns ($`w_r = d`$, $`k_f = 0.55`$) and the widened outline for `outline` returns ($`w_r = d_\mathrm{eff}`$,
  $`k_f = 1.2`$, fitted to only two runs).
* *The cube*, $`0.76\,x^3`$, fitted together with $`k_f`$ to the runs with $`\Lambda \le 0.3`$. For stems widened by the footprint
  ($`d/d_\mathrm{eff} \le 0.44`$) the fit is within 1.5% up to $`\delta = 1`$; for unwidened stems ($`d = d_\mathrm{eff}`$) it is 2.3% low at $`\delta = 0.45`$ and
  up to 10% off at $`\delta`$ = 0.8–1.
* *Density term*, $`0.35\,x\,\tau`$, with $`\tau`$ the optical depth of the exit distance and $`N d_\mathrm{eff}`$ taken from the occlusion
  neighbourhood. With it, the factor is within 1.2% for widened stems at $`\delta`$ = 0.3 and 0.45 up to $`\Lambda = 1`$. It is valid to
  $`N d_\mathrm{eff} W \approx 1`$, about 900 stems/m² for single 5 cm voxels and a 2.25 cm effective width: at $`\Lambda = 2`$ the Monte Carlo
  factor was 1.68–1.88 ($`\delta = 0.45`$) against 1.33 fitted, and beyond $`\Lambda \approx 1`$ it grows without bound. This is a hard limit
  of single-voxel windows in dense stands, not only of the fit. The $`\max(0.2, \dots)`$ caps the first part of the factor at 5.

The derivation is exact for one isolated tall stem; the three terms are empirical, and the density term is specific
to the conditions simulated.

#### 8.1.2 Small windows (sw, experimental)

Single-voxel windows (`--pool_radius 0`) resolve patches and root lines instead of averaging across them, but with
the first-order correction alone (stage 1, "p0 ss") they read 1.024 for uniform stems with the standard rays and
0.923 with a quarter of them (300 stems/m², 0–0.1 m). Four effects add up in single voxels. Their sizes come from
the small-window study's ablations, each restoring one of them in the configuration below:

**Table 13.** Effects in single-voxel windows (300 stems/m², 0–0.1 m) at ×0.25 (×1) of the standard ray density, and
the option that removes each.

| Effect | Effect on `density` | Option |
|---|---|---|
| Neighbour top-up: voxels with fewer than 10 rays are blended with their neighbours by path, which under-weights the dense ones | −11% (−0.7%); −7 to −14% at ×0.25 over 300–1500 stems/m² | `--no_neighbour_priors` |
| Occlusion response inverted per voxel: the voxel's solid fraction is noisy and $`a\,R(a)`$ is convex | −10% (−6%) | `--occlusion_radius 0.05` (invert $`R`$ over the 3 × 3 voxels and apply its factor) |
| First-order factor in place of $`1/(1-x)`$: the factor is 4.6% too small | +5.3% (+5.3%) | `--self_shadow_second_order` |
| Finite-ray ratio bias: a ratio of sums over a few rays is biased high | +12.5% (+3.8%) | `--ratio_bias_correction` (delta method over the rays' segments) |

The configuration that combines the four, "sw", is

```console
--pool_radius 0 --self_shadow_correction --self_shadow_second_order --occlusion_radius 0.05 --ratio_bias_correction --no_neighbour_priors
```

(with `--estimator pooled`, `--stem_diameter` and the footprint as usual). It does not change the run time, and the
results are deterministic. In this configuration `--occlusion_radius` also removed the −1.9% of the pooled value in
the lowest layer at 300 stems/m², but not the +1% at 1000 stems/m² (Section 7.9).

**Table 14.** `density` / truth over all cells, 0–0.1 m, at ×0.25 / ×0.5 / ×1 of the standard ray density (4 seeds,
4 m plots), and the observed fraction of the layer with sw (×0.25 / ×1). Cable roots: 3 cm lateral spread; patches:
30% cover.

| Stems | Default | ss | sw | Observed, sw |
|---|---|---|---|---|
| Random, 30/m² | 1.088 / 1.067 / 1.073 | 1.015 / 0.995 / 1.000 | 1.006 / 0.984 / 0.991 | 1.00 / 1.00 |
| Random, 100/m² | 1.069 / 1.074 / 1.075 | 0.999 / 1.003 / 1.004 | 1.000 / 0.993 / 0.997 | 1.00 / 1.00 |
| Random, 300/m² | 1.055 / 1.066 / 1.067 | 0.990 / 0.999 / 1.001 | 0.973 / 0.996 / 0.999 | 0.91 / 1.00 |
| Random, 600/m² | 1.049 / 1.055 / 1.070 | 0.988 / 0.993 / 1.008 | 0.910 / 0.973 / 1.005 | 0.51 / 0.96 |
| Random, 1000/m² | 0.976 / 1.027 / 1.041 | 0.923 / 0.972 / 0.985 | 0.750 / 0.858 / 0.932 | 0.14 / 0.68 |
| Random, 1500/m² | 0.892 / 0.974 / 1.011 | 0.849 / 0.926 / 0.961 | 0.707 / 0.757 / 0.820 | 0.02 / 0.26 |
| Patches, 0.2 m radius | 0.982 / 0.984 / 0.986 | 0.918 / 0.920 / 0.922 | 0.978 / 0.988 / 0.989 | 0.99 / 1.00 |
| Patches, 0.3 m radius | 0.996 / 1.006 / 1.007 | 0.932 / 0.942 / 0.942 | 0.967 / 0.997 / 0.989 | 0.99 / 1.00 |
| Patches, 0.5 m radius | 1.019 / 1.024 / 1.022 | 0.954 / 0.959 / 0.957 | 0.967 / 1.000 / 0.998 | 0.98 / 1.00 |
| Cable roots, 150/m² | 0.847 / 0.867 / 0.867 | 0.792 / 0.811 / 0.811 | 0.921 / 0.970 / 0.978 | 0.98 / 1.00 |
| Cable roots, 500/m² | 0.386 / 0.432 / 0.470 | 0.362 / 0.405 / 0.441 | 0.483 / 0.645 / 0.760 | 0.82 / 0.92 |

At 0.1–0.2 m, sw read 0.995–1.012 for random stems from 30 to 1000 stems/m² with ×0.5 and ×1 rays, 0.994–1.018 in
patches and 0.964–0.990 along cable roots at 150 stems/m².

**Strengths.** sw is the best single-window estimator found for clustered stems (the hybrid of Section 8.1.3 uses it): cable roots at 150 stems/m² read 0.92–0.99
(0.978 / 0.987 at 0–0.1 / 0.1–0.2 m with the standard rays), against 0.85–0.90 with the default and 0.79–0.84 with
ss, and patches read 0.967–1.018. For uniform stands it read 0.98–1.02 over the cells at least 0.9 observed in every
case of the study (random stems at 30–1500 stems/m², 0.5 and 2 cm stems, no footprint and a wide one, stems of one
height), and with the standard rays it is within 1% at 0–0.2 m for random stems at 30–600 stems/m². It is a little
less sensitive to the assumed diameter: ±20% gave −8.8% / +6.1% (300 stems/m², 0–0.1 m), against −10.4% / +8.1%
with ss.

**Failures.** sw is worse than ss for dense stems under sparse scans. Without the neighbour top-up, voxels with 2 rays
or fewer are unobserved; these are mostly the dense voxels, so the voxels that remain are a selection that reads low
(Section 8.3). At 1000 stems/m² the lowest layer was 0.14, 0.38 and 0.68 observed with sw at ×0.25, ×0.5 and ×1 rays
(0.62, 0.80 and 0.92 with the default), and read 0.75, 0.86 and 0.93; at 1500 stems/m², 0.71–0.82. The same happens
with a wide footprint at 1000 stems/m² (0.74 with the standard rays), 2 cm stems at 300 stems/m² with a quarter of
the rays (0.85) and stems of one height (0.6 m) at 300 stems/m² with a quarter of the rays (0.82 at 0–0.1 m); at
1000 stems/m² with stems of one height every layer up to 0.4 m read 0.33–0.91. Cable roots at 500 stems/m² still
read 0.48–0.76. Beyond about 900 stems/m² at 5 cm voxels the second-order factor's density term is also outside its
range. So sw needs each voxel's neighbourhood to be well observed. Overlapping stems read as with ss (0.969), and sw
has not been validated with `outline` returns, or with 2.5 cm voxels, which are ray-starved: of the small windows
tried, only 5 cm voxels with pool radius 0 are usable.

**Alternatives rejected for small windows.**

* Falling back to a larger window when a voxel has few rays: the selection is correlated with the voxel's content.
* The same fallback triggered by the neighbourhood's ray count: patches and cable roots lose 5–20%.
* Treating voxels with only 1–2 rays as observed.
* An expected-path estimator: its cross-voxel self-shadow is +9 to +17%, and cable roots read 0.91.

The hybrid of Section 8.1.3 keeps sw's strengths and avoids these failures by choosing per cell.

#### 8.1.3 The hybrid (recommended)

ss is accurate in uniform stands, also dense and poorly penetrated ones, but reads low in clumped ones; sw resolves
clumps but reads low where rays are sparse. The hybrid (`--self_shadow_correction --hybrid`) uses each where it
works.

**Why a per-voxel switch fails.** The obvious rule, sw in voxels with enough rays and ss elsewhere, is biased,
because a voxel's ray count carries a content signal: stems shadow their neighbours, so dense places get fewer
rays. From voxel dumps of the block K scenes (per 5 cm layer, ray counts of the sw estimate):

* at 300 stems/m² with a quarter of the standard rays (lowest 5 cm), voxels with 2 rays or fewer held 1.37 times the
  layer's mean density, and voxels with 10–20 rays 0.63 times;
* voxels holding a stem got 0.81–0.90 of the rays of empty voxels at 100–300 stems/m² (10–19% fewer), and 0.64–0.83
  at 1000 stems/m² (17–36% fewer);
* within each ray-count bin, sw followed the truth (where it observed the voxel), while ss was pulled towards the
  neighbourhood mean: in the example above, ss read 1.18 against a truth of 1.37 in the bin of 2 rays or fewer, and 0.80
  against 0.63 in the 10–20-ray bin.

ss's errors therefore cancel over a cell, but not within a ray-count bin. A per-voxel weight $`w(n)`$ on sw keeps ss's
error $`(1-w)(\mathrm{ss} - \mathrm{truth})`$ where $`n`$ is small, which is negative, and drops the compensating positive error where $`n`$
is large. Ramping the weight from 0 at 2 rays to 1 at 10, 20 or 40 rays gave 0.994, 0.978 and 0.961 at 300
stems/m² with the standard rays (the hybrid 0.999), and along cable roots at 150 stems/m² 0.962, 0.944 and 0.905
with the standard rays (sw 0.978) and 0.874, 0.814 and 0.797 with a quarter of them (sw 0.921, the hybrid 0.938).
The stage-2 fallback triggered by the neighbourhood's ray count lost 5–20% in the same way (Section 8.1.2).

**The per-cell rule.** The two estimates are combined in two steps, both from ray counts only, never from hits:

1. Per voxel: the sw density where sw observes the voxel (more than 2 rays, without the neighbour top-up), the ss
   density otherwise.
2. Per cell: these voxel values if sw observes at least 0.9 of the cell's voxels, otherwise the ss values
   throughout.

The sw part is the configuration of Section 8.1.2 (pool radius 0, the second-order factor, the occlusion response
from the `--pool_radius` neighbourhood, the finite-ray correction, no neighbour top-up), and the ss part is the
default window with the first-order correction. The 0.9 is `DensityOptions::hybrid_min_observed` in the library;
it is not a command-line option.

**The selection effect.** Because ray counts carry content, the cells that pass the 0.9 test are a selection. In
the uniform scenes where both kinds of cell occurred (300 stems/m² with a quarter of the rays, 600 and 1000 stems/m²;
0–0.2 m), the cells taking sw were 7–11% less dense than those taking ss. Each estimator is nearly unbiased on its
own subset (sw-cells: hybrid 1.006–1.017 of their truth; ss-cells: ss 0.979–0.993), so the total is unbiased: the
hybrid's plot totals were within 0.5% of ss's in these scenes. A per-cell choice moves whole cells between two
estimators that are each right for the cell, rather than mixing an estimator's errors by content.

**Outputs.** `density` is the hybrid; `pooled_density`, `hits`, `weighted_path` and `observed_fraction` are those of
ss. `hybrid_weight` (Section 6) is 0 for ss cells and at least 0.9 for sw cells. The rendered image uses the per-voxel
values of step 1 in every column (it has no cells).

**Cost.** One pass: the voxel grid is copied after the ray walk and both estimates are computed from it, so the run
needs about twice the memory of one estimate (from the code) and 5% more time than ss (0.50 s against 0.48 s per
render, mean over the study's 504 runs; the default took 0.46 s). The tiled output equals a single run with the
existing halo (the sw part's occlusion neighbourhood is the pool radius, so the default halo covers it), the result
is deterministic, and `--hybrid` is refused without `--self_shadow_correction` or with any of the four single-voxel
options (tests `RayDensityHybrid` and `RayDensityObservedMask`).

**Table 16.** `density` / truth over all cells, 0–0.1 m, at ×0.25 / ×0.5 / ×1 of the standard ray density (4 seeds,
4 m plots; the scenes of Table 14).

| Stems | Default | ss | sw | Hybrid |
|---|---|---|---|---|
| Random, 30/m² | 1.088 / 1.067 / 1.073 | 1.015 / 0.995 / 1.000 | 1.006 / 0.984 / 0.991 | 1.006 / 0.984 / 0.991 |
| Random, 100/m² | 1.069 / 1.074 / 1.075 | 0.999 / 1.003 / 1.004 | 1.000 / 0.993 / 0.997 | 1.002 / 0.993 / 0.997 |
| Random, 300/m² | 1.055 / 1.066 / 1.067 | 0.990 / 0.999 / 1.001 | 0.973 / 0.996 / 0.999 | 0.992 / 0.999 / 0.999 |
| Random, 600/m² | 1.049 / 1.055 / 1.070 | 0.988 / 0.993 / 1.008 | 0.910 / 0.973 / 1.005 | 0.988 / 0.994 / 1.013 |
| Random, 1000/m² | 0.976 / 1.027 / 1.041 | 0.923 / 0.972 / 0.985 | 0.750 / 0.858 / 0.932 | 0.923 / 0.972 / 0.985 |
| Random, 1500/m² | 0.892 / 0.974 / 1.011 | 0.849 / 0.926 / 0.961 | 0.707 / 0.757 / 0.820 | 0.849 / 0.926 / 0.961 |
| Patches, 0.2 m radius | 0.982 / 0.984 / 0.986 | 0.918 / 0.920 / 0.922 | 0.978 / 0.988 / 0.989 | 0.994 / 0.990 / 0.990 |
| Patches, 0.3 m radius | 0.996 / 1.006 / 1.007 | 0.932 / 0.942 / 0.942 | 0.967 / 0.997 / 0.989 | 0.989 / 1.000 / 0.989 |
| Patches, 0.5 m radius | 1.019 / 1.024 / 1.022 | 0.954 / 0.959 / 0.957 | 0.967 / 1.000 / 0.998 | 1.001 / 1.006 / 0.998 |
| Cable roots, 150/m² | 0.847 / 0.867 / 0.867 | 0.792 / 0.811 / 0.811 | 0.921 / 0.970 / 0.978 | 0.938 / 0.976 / 0.980 |
| Cable roots, 500/m² | 0.386 / 0.432 / 0.470 | 0.362 / 0.405 / 0.441 | 0.483 / 0.645 / 0.760 | 0.433 / 0.556 / 0.707 |

At 0.1–0.2 m the hybrid read 0.976–1.016 for random stems (0.976 at 1500 stems/m² with a quarter of the rays, as ss),
0.998–1.024 in patches, 0.973–0.992 along cable roots at 150 stems/m² and 0.615–0.901 at 500. Over the cells at
least 0.9 observed, cable roots at 500 stems/m² read 0.467 / 0.596 / 0.749 at 0–0.1 m.

**Table 17.** The hybrid's errors, the basis of the recommendation (Section 10).

| Case | Hybrid | ss | Source |
|---|---|---|---|
| Uniform stands where ss passes ±2% | within ±2%; worst 1.019 ± 0.004 (0.5 cm stems, 300/m², ×0.25 rays, 8 seeds; ss 1.011) | within ±2% | sets H, H2, H3; block S |
| Patches (0.2–0.5 m radius, 30% cover, 300/m²) | 0.989–1.006 (0.998–1.024 at 0.1–0.2 m) | 0.918–0.959 | set H |
| Cable roots, 150/m² | 0.938–0.980 (0.979 with lines 2 m apart) | 0.792–0.811 | set H; block S |
| Cable roots, 500/m² | 0.433–0.707 (0.467–0.749 over cells ≥ 0.9 observed) | 0.362–0.441 | set H |
| Dense, poorly penetrated stands (1000–1500/m², stems of one height) | the same as ss, including its shortfall (e.g. 0.923 at 1000/m², ×0.25 rays) | | sets H, H2, H3 |
| Range noise 1 cm (300/m²) | 0.990 / 0.990 | 1.000 / 1.008 | block S |
| Range noise 2 cm | 0.970 / 0.966 | 0.999 / 1.001 | block S |
| Pose noise 2 cm (knots every 0.5 m) | 0.944 / 0.932 | 0.998 / 0.999 | block S |

(Noise rows: d/t90 at 0–0.1 / 0.1–0.2 m, 4 seeds.) The other block S scans (lines 2 m apart, scanner at 2 m, at 300
and 1000 stems/m²) read 0.987–1.013 with the hybrid where the cells were at least 0.9 observed (1000 stems/m² with
lines 2 m apart had none at 0–0.1 m). The noise sensitivity is sw's: sw read the same as the hybrid within 0.001 in
the noise scenes, where sw observed at least 0.99 of every cell, so every cell took sw. Its cause has not been
isolated (probably noise moving hits and path between single voxels, which a 0.15 m window averages out), and
nothing in the hybrid's rule detects noise, so where the range or pose error may exceed about 1–2 cm use ss instead.

### 8.2 Stem tops and the return convention

The model has no stem tops. Steep rays end on stem tops, which counts as interception by the side of the stem, and
stems end in every layer above 0.1 m for the height distribution simulated. The generator's treatment of the tops
is itself a choice: the top is a solid disc of the widened outline, and a side hit is tested at the height where the
beam enters the widened outline but, under `surface`, reported lower.

The return convention interacts with the tops. Deep inside a stand of uniform height, the longer paths of `surface`
returns are balanced by as many extra returns carried down from the more strongly illuminated layer above, so $`R`$ and
the estimates are convention-independent. Where stems end, there are no returns to carry down from above, so a
`surface` scene reads lower than an `outline` scene, by up to $`N d_\mathrm{eff} R`$ × 7 mm in the layer holding the tops (5% at
300, 15–19% at 1000 stems/m²). This was confirmed with stems of uniform height 0.6 m (one scene each, 1 cm stems):
below the top layer the conventions agreed within about 1% at 300 stems/m² (`pooled_density` 1.023/1.031,
1.009/1.016, 0.993/1.004, 0.977/0.980, 0.979/0.984 for 0–0.5 m) and within about 2% at 1000, while in the layer
holding the tops `outline` read 7.5% (300) and 15% (1000) higher. With the default heights, the difference grows
with height:

| `pooled_density`, one scene (seed 101) | 0–0.1 m | 0.1–0.2 m | 0.2–0.3 m | 0.3–0.4 m | 0.4–0.5 m |
|---|---|---|---|---|---|
| 300/m², `surface` | 0.968 | 0.982 | 1.002 | 1.079 | 1.180 |
| 300/m², `outline` | 0.977 | 1.007 | 1.054 | 1.162 | 1.300 |
| 1000/m², `surface` | 1.009 | 0.999 | 0.952 | 0.999 | 1.172 |
| 1000/m², `outline` | 1.021 | 1.040 | 1.014 | 1.079 | 1.285 |

So the layers near the stem tops are convention-dependent: every result above about 0.1 m depends on the return
convention and on the treatment of stem tops, by 5–15%, and is not validated. In `density` at 300 stems/m², the
upper layers read about 2%, 5%, 11% and 18% above the lowest layer for 0.1–0.2, 0.2–0.3, 0.3–0.4 and 0.4–0.5 m
(`surface`, one scene). Near-horizontal rays alone are no remedy (Section 7.4).

The validation matrix separated the tops from the rest (Section 7.9). Top-disc returns were 0.3%, 2.5%, 6.7%, 12%
and 15% of the stem returns in the 0.1 m layers up to 0.5 m (one reference scene), and with ss `density` read 1.001,
1.010, 1.052, 1.117 and 1.203 there; stems of one height read 1.000–1.007 in every layer below their tops. So with
the self-shadowing corrected, the remaining excess above 0.1 m is mostly the stem tops (5–20% at 0.2–0.5 m). The
tops also account for most of the failures with `outline` returns, 2 cm stems and 5 cm pose error. A model of the
tops (an end-cap term in $`G`$, or excluding returns at stem tops) has not been tried; real stem tops, which taper and
are rounded, need checking against the dowel arrays of Section 9.2.

### 8.3 Poorly observed cells read low

Where rays penetrate poorly, the estimates read low. A stage-1 analysis at voxel resolution
(`scripts/pneumatophore_low_observation_study.py`: dense stems, stems of one height or taller than the scanner, and
2 m line spacing; 8 scenes on 4 m plots, binned over the 1 × 1 × 0.1 m cells below the stem tops) separates two
components:

1. *Selection.* The observed voxels are the less occluded places. Below an observed fraction of 0.9, the true
   density of the windows they sample was 1.4–4.1% below the cell's (0.959–0.986), so even an exact estimate of the
   observed part reads low.
2. *Poor penetration.* Windows that are observed but poorly penetrated also read low: with the self-shadow
   correction, their estimates were 0.870–0.972 of the truth of the windows sampled, below an observed fraction of
   0.9.

**Table 10.** `density` / truth binned by the cell's `observed_fraction` (truth-weighted means over the cells).

| Observed fraction | Cells | Without correction | With `--self_shadow_correction` |
|---|---|---|---|
| 0.9–1.0 | 64 | 1.062–1.063 | 1.005–1.008 |
| 0.5–0.9 | 203 | 0.997–1.012 | 0.945–0.957 |
| below 0.5 | 85 | 0.878–0.954 | 0.834–0.906 |

Without the correction, cells observed 0.5–0.9 read right only because the self-shadowing bias cancels this low
bias. The same cancellation gives the "+1%" at 1500 stems/m² in Table 2 (observed fraction 0.68; in the stage-1
scenes 1.011 without the correction and 0.961 with it). With stems of one height (0.6 m) at 1000 stems/m², the
observed fraction was 0.28, 0.41, 0.58, 0.79 and 0.97 in the 0.1 m layers from the ground, and `density` read
0.893, 0.960, 0.988, 1.022 and 1.044 without the correction and 0.848, 0.911, 0.937, 0.969 and 0.989 with it
(stage-1 scenes, 4 seeds).

Neither remedy tried helped: filling unobserved voxels from their observed neighbours changed the binned values by
−0.6 to +1.2 points, and scaling the self-shadow factor by the window's observed voxel fraction recovered only
0.6–2.1 points below 0.9. Cells below an observed fraction of 0.9 are therefore masked rather than corrected. The
validation matrix confirmed the threshold over a wider range of scans: with ss, cells at least 0.9 observed read
0.995 / 1.010 (0–0.1 / 0.1–0.2 m), and every bin below 0.9 read low, from 0.97 at 0.8–0.9 to 0.78–0.83 below 0.3
(Section 7.9). sw leaves more cells below 0.9, because it has no neighbour top-up (Section 8.1.2); the hybrid's
`observed_fraction` is that of ss.

**The mask.** `--min_observed_fraction 0.9` (Section 6; off by default in both `rayrender` and `raydensity`) sets
the densities of the cells below 0.9 to NaN, with or without the correction, and keeps their `observed_fraction`,
`hits`, `weighted_path` and `hybrid_weight`. 0.9 is the recommended value. Caveat: the mask applies to each cell or
layer as a whole. Profile layers (`--profile`) span the whole plot, and with absolute heights over uneven ground a
low layer includes voxels below the ground, which no ray reaches, so its observed fraction falls below 0.9 and the
mask blanks it although its voxels above the ground are well observed. The same holds for cells crossed by the
ground surface. Over uneven ground, use heights above the ground (`raydensity --ground`, Section 6) when masking;
`rayrender` has no `--ground`, so there mask by hand only the layers and cells that lie wholly above the ground.

The stage-1 scenes used one scan geometry (scanner at 1.2 m, lines 1 m apart, about 56 000 rays/m²) and stems of
0.25 m mean height, which is a favourable case. In stage 1, lines 2 m apart from a scanner at 2 m kept the lowest
layer fully observed at 300 stems/m² (1.064 without the correction, 0.992 with it), but at 1500 stems/m² lines 2 m
apart left it 0.38 observed (0.975 and 0.927). In the validation matrix the lowest layer at 1000 stems/m² was 0.92
observed with lines 1 m apart, and 0.66 and 0.59 with lines 2 and 3 m apart (0.82, 0 and 0.36 of its cells at least
0.9 observed); with stems of one height (0.6 m) at 600–1500 stems/m², no cell of it was. Field scans with wider line
spacing, longer ranges, or taller or denser stems will be more penetration-limited.

### 8.4 Inputs: diameter and threshold

The estimate is about as sensitive to the assumed mean diameter (Table 5: −11% / +9% for 0.8 / 1.2 cm at 300
stems/m², −8% / +5% at 1000) as to the detection threshold (Table 4: −11 to −12% / +22 to +27% for 0.02 / 0.15 when it
is 0.05). The validation matrix (block D, ss, 0–0.1 m, 4 seeds) shows how the diameter sensitivity depends on the
stem size relative to the footprint: an assumed diameter 10% low or high changed `density` by −6.6% / +6.3% for
0.5 cm stems, −4.9% / +4.3% for 1 cm stems (−3.6% / +2.7% at 1000 stems/m²) and −1.6% / +1.1% for 2 cm stems, all at
300 stems/m² unless stated; ±20% gave −13.6% / +12.2%, −10.4% / +8.1% and −3.9% / +1.7%. Thin stems are seen mostly
through the footprint, so the hit rate is nearly proportional to $`N`$, and $`N d`$ inherits the relative error in $`d`$.
A threshold calibrated on a counted plot also absorbs any error in the diameter (and, without the correction, the
bias of Section 8.1). Block L (Section 7.10) measured how far such a calibration carries: with ss, a threshold fitted
to compensate a diameter error of ±10–20% held within 1% across layers, ranges of 1–4 m and 2 cm stems (the last
two from the model), but spread by 4% at 1000 stems/m², and with the base estimator it carried the +7% self-shadowing
bias with it. It is an empirical fit, not the physical threshold, so it must not be reused for another scanner,
stem size or density (Section 9.2). The beam's exit diameter and divergence act on the same effective width:
Section 7.10 gives their sensitivities and Section 9.3 turns all of these into the accuracy a field campaign needs.

### 8.5 `pooled_density` is path-weighted horizontally and vertically

`pooled_density` weights each voxel by its path, and rays penetrate dense parts less, so it reads low wherever the
density varies within the cell: about half the truth for patches (Table 3), and also in the vertical, where the
upper voxels of a layer in which the density falls with height carry more path and less density. This is why it
read 0.95 at 0.2–0.3 m at 1000 stems/m² (Table 1): per 5 cm layer it read 1.043 (0.20–0.25 m, weighted path 23 400)
and 1.056 (0.25–0.30 m, 54 500), with no dip. It is within 3% only where the stems are uniform within the cell
horizontally and vertically, and it also carries the stem-top bias: at 300 stems/m² it read 1.075 at 0.3–0.4 m
(Section 7.4) and, per 5 cm layer, rose steadily from 0.966 at the ground to 1.224 at 0.45–0.50 m (one scene). Even
for stems of one height its error depends on the layer (1.027 at 0–0.1 m to 0.976 at 0.4–0.5 m, 300 stems/m²;
Section 7.9), so `density` / `pooled_density` is only a rough self-shadow signature.

### 8.6 Clustering at scales below the pooling window

Stems clustered more tightly than the window read low, as any free-path estimate does when pooled over uneven
cover. This is the normal arrangement of pneumatophores, along the cable roots from which they rise, so it matters.
At the default window, cable roots with a 3 cm lateral spread read 0.87 (and 0.77 with a 0.25 m window), with a
`patchiness_ratio` of 1.6–1.8 in the lowest 0.2 m (Table 3). The validation matrix (block G, 4 seeds, 0–0.1 /
0.1–0.2 m) varied the spread:

| Cable roots | Spread 1.5 cm | 3 cm | 6 cm |
|---|---|---|---|
| 150 stems/m², default | 0.812 / 0.830 | 0.867 / 0.898 | 0.984 / 1.013 |
| 150 stems/m², ss | 0.759 / 0.776 | 0.811 / 0.840 | 0.921 / 0.948 |
| 150 stems/m², default, `--pool_radius 0.1` | 0.725 / 0.741 | 0.757 / 0.789 | 0.888 / 0.930 |
| 500 stems/m², default | 0.387 / 0.441 | 0.475 / 0.587 | 0.751 / 0.884 |
| 500 stems/m², ss | 0.363 / 0.413 | 0.446 / 0.552 | 0.707 / 0.833 |

So the default's −13% at a 3 cm spread is not a general figure: the bias ranges from −19% to near zero over 1.5–6 cm
at 150 stems/m², and is −25 to −61% at 500 stems/m², and part of it is cancelled by the self-shadowing bias (with
ss it is 2–6 points lower). Patches (0.2–0.5 m radius, 30% or 60% cover, stems between patches 0 or 0.2 of the
density within) read 0.986–1.051 with the default and 0.922–0.984 with ss.

Single-voxel windows change this picture (Section 8.1.2): with sw, cable roots at 150 stems/m² (3 cm spread) read
0.978 / 0.987 and patches 0.989–1.005 with the standard rays, because the window no longer straddles root lines and
patch edges. At 500 stems/m² sw still read 0.76 / 0.91, and it needs the lower layers to be well observed. Other
spreads than 3 cm have not been run with sw or the hybrid; the hybrid read 0.938–0.980 at 150 stems/m² (3 cm spread,
×0.25–×1 rays) and 0.979 with lines 2 m apart (Section 8.1.3).

### 8.7 Summary of residual biases (synthetic, 5 cm voxels, 1 × 1 × 0.1 m cells)

| Bias | Size | Applies to | Section |
|---|---|---|---|
| Self-shadowing | +4–9% (first order $`0.4732\,d_\mathrm{eff}/W`$; $`1/(1-x)`$ to second order) | `density`, default | 8.1, 7.9 |
| Residual with ss | 0.983–1.006 at 0–0.1 m and 0.997–1.020 at 0.1–0.2 m over the scan geometries of the matrix (partly a cancellation at 0–0.1 m) | `density`, ss | 8.1, 7.9 |
| Heterogeneity within the window (patch edges, root lines) | hidden by the self-shadowing bias in the default; with ss, patches 0.92–0.98, cable roots 0.76–0.92 (150/m²) | `density` | 8.1, 8.6 |
| Sparse rays in dense stands | sw 0.75–0.93 at 1000/m², 0–0.1 m (×0.25–×1 rays); the hybrid takes ss there | `density`, sw | 8.1.2 |
| Hybrid | uniform within ±2% where ss passes (worst 1.019); patches 0.989–1.006; cable roots 0.94–0.98 (150/m²), 0.43–0.71 (500/m²); as ss in dense, poorly penetrated stands | `density`, hybrid | 8.1.3 |
| Range or pose noise, with the hybrid | 1 cm range noise 0.990; 2 cm range noise 0.966–0.970; 2 cm pose noise 0.93–0.94 (ss passes) | `density`, hybrid, sw | 8.1.3 |
| Beam footprint error | about −1.0 × the error in the mean $`d_\mathrm{eff}`$ (ss; −1.2× sw, −0.95× base); $`w_0`$ ±25%: −9 / +11%; $`\gamma`$ ±25%: −4 / +5%; datasheet misreadings: $`w_0`$ +18 / −21%, $`\gamma`$ +10 / −15%; matched only at 2 m: −13 to +18% at 3.5 m | all | 7.10 |
| Stem tops | 5–20% at 0.2–0.5 m with ss, growing with height | both | 8.2, 7.9 |
| Return convention (`outline` data) | +1.2–1.9% where no stems end, +3–4% more at 0.1–0.2 m; not removed by `--return_range` | both | 7.9, 8.2 |
| Poor penetration (observed fraction below 0.9) | default about 1.00–1.03 at 0.5–0.9 by cancellation, 0.82–0.95 below 0.5; ss 0.94–0.97 and 0.78–0.91 | both | 8.3, 7.9 |
| Pose error 5 cm | −3 to −11% (ss), where stems end | both | 7.9 |
| Range noise 3 cm | ss 0.994 / 0.993 (300/m²), 0.973 / 0.984 (1000/m²), at 0–0.1 / 0.1–0.2 m | both | 7.9 |
| Diameter error ±10% (±20%) | 1 cm stems −4.9% / +4.3% (−10.4% / +8.1%), ss, 300/m² | both | 8.4 |
| Diameter mixture (6/14 mm) | +0.7 to +1.3% | both | 7.5, 7.9 |
| Threshold error (0.4×, 3×) | −11 to −12% / +22 to +27% | both | 8.4 |
| Minimum gap ignored ($`g = 1`$) | +4% (100/m²) to +34% (1000/m²) | both | 5 |
| Overlapping stems | ss 0.968 (300/m²), 0.874 (1000/m²): the occlusion correction assumes no overlap | both | 7.9 |
| Patchiness within the cell | about −50% | `pooled_density` | 8.5 |
| Vertical path weighting | e.g. −5% at 0.2–0.3 m (1000/m²) | `pooled_density` | 8.5 |
| Clustering below the window (cable roots) | default −19% to −2% at 150/m² (spread 1.5–6 cm), −25 to −61% at 500/m²; ss 2–6 points lower; sw −2% (3 cm, 150/m²) | `density` | 8.6 |

## 9. What was not simulated, and a calibration and validation protocol

### 9.1 Not simulated

The synthetic validation uses the same beam, stem and return models as the corrections, so it tests the estimator's
consistency, not the physical models. The validation matrix (Section 7.9) added range noise up to 3 cm, smooth
trajectory error up to 5 cm, diameter error and two diameter distributions, five height distributions, scan
geometries from 0.5 to 3 m line spacing and 0.8 to 2 m scanner height, and overlapping or spaced stems; block L
(Section 7.10) added a footprint and threshold that differ from the ones assumed, but only within the Gaussian beam
model. It did not include:

* A beam that is not Gaussian, or whose threshold varies with range, other scanners' error models, the hybrid with
  a mismatched footprint (block L ran base, ss and sw; with Eq. (15) the hybrid should lie between ss and sw), and
  the effect of range and pose error on trunk removal.
* Mixed stem and ground returns at the stem base, and imperfect ground classification: the ground is flat (except
  for the slope test) and perfectly labelled, yet the main layer is 0–0.1 m.
* Stem taper, lean, curvature and branching.
* Reflectance and incidence dependence of the threshold, wet or muddy stems, and the scanner's actual return
  detection (multiple returns, return mode).
* Real stem tops (the generator's tops are a modelling choice, Section 8.2), which cause most of the remaining error
  above 0.1 m.
* Water, debris, seedlings, crabs and other non-stem returns in the stem layer.

None of the method has been checked against field data.

### 9.2 Calibration and validation protocol

The aim is to calibrate the beam model from the physics (not from the density result) and then validate the whole
chain on independent data of known $`N d`$. Section 9.3 gives the accuracy each measurement needs.

1. **Beam footprint, threshold and return convention (rod test).** Scan isolated vertical rods or dowels of known
   diameter, with a surface like the stems' (painted or wetted to match their reflectance), and measure the
   effective width $`d_\mathrm{eff}(r)`$ at each range as twice the half-width of the band of beam offsets that return from the
   rod (the hit fraction of the rays passing within a few centimetres of its axis). Fit $`w_0`$, $`\gamma`$ and $`\tau`$ of Eq. (4)–(6)
   to all the $`d_\mathrm{eff}(r)`$ points; the design below follows from a Fisher analysis of this fit and a synthetic rod test
   (`scripts/pneumatophore_rod_calibration.py`).

   * *Degeneracies.* With one rod diameter, the beam scale and the threshold cannot be separated: a wider beam with
     a higher threshold gives the same $`d_\mathrm{eff}`$ at every range (correlations of the fitted parameters ≥ 0.998, even
     with 9 ranges from 1 to 5 m). Scanning at several ranges separates only the exit diameter from the divergence.
     So a single-rod, single-range calibration of the three parameters is not possible, and with $`w_0`$ and $`\gamma`$ taken
     from the datasheet it only moves $`\tau`$.
   * *Several rod diameters* break the degeneracy, because the footprint widens thin rods relatively more: use rods
     of about the survey's mean stem diameter, half of it and twice it, plus thin rods (3–5 mm), which constrain the
     threshold best. Fitted standard errors ($`\ln w_0`$ / $`\ln\gamma`$ / $`\ln\tau`$, at 0.5 mm per $`d_\mathrm{eff}`$ point): 0.058 / 0.079 / 0.135
     for 0.3, 0.5 and 1 cm rods at 1, 3 and 5 m; 0.15 / 0.18 / 0.50 for 0.5, 1 and 2 cm rods at 1, 2 and 4 m.
   * *Ranges* of about 1, 2.5 and 4–5 m, spanning the survey's, with 1–5 × 10⁴ beams per rod and range (12 500–54 000
     in the synthetic test).
   * *Joint fit and χ² check.* Fit all three parameters jointly to all rods and ranges, and use χ² per degree of
     freedom to flag a datasheet value that is wrong: with $`w_0`$ fixed 25% high, a threshold-only fit gave χ²/dof 4.7,
     and with $`\gamma`$ fixed 50% high 8.2 (1.3 with the right values).
   * *What one rod gives.* The density needs only the mean $`d_\mathrm{eff}/d`$ at the survey's ranges (Eq. 15), and a single
     1 cm rod at 1–5 m pins that to 0.9% (s.e.), although each parameter on its own is unconstrained. In the
     synthetic test the three rods (0.5, 1, 2 cm) recovered 10.01 mm, 3.08 mrad and $`\tau = 0.053`$ (true 10 mm, 3 mrad,
     0.05), while the 1 cm rod alone gave 13.6 mm, 5.4 mrad and $`\tau = 0.137`$, wrong parameters but the right mean
     $`d_\mathrm{eff}/d`$ (0.999 of the true value over 1–3 m, 1.000 at 4 m).
   * *The short-range trap.* A rod scanned only at 1.5–2.5 m hid a divergence error: with $`\gamma`$ 50% high and $`\tau`$ fitted,
     χ²/dof was 0.09 and the mean $`d_\mathrm{eff}`$ over 1–3 m was right to 0.3%, but at 4 m it was 4.3% high, so the density
     there would read about 4% low. Never calibrate at a single range or a narrow band of ranges.
   * *Datasheet pitfalls.* Eq. (4) needs the 1/e² exit diameter and the full divergence angle. Reading an exit
     diameter quoted as FWHM as if it were 1/e² (the 1/e² diameter is 1.70 times the FWHM) gives +18% with ss
     (default beam, 300 stems/m²), and the reverse conversion −21%; taking a quoted half-angle divergence as the full
     angle gives +10%, and the reverse −15% (interpolated from Table 15).

   Propagate the fitted mean $`d_\mathrm{eff}/d`$ through Eq. (15), then run the known-density array test of step 2.

   *Open check.* In the synthetic rod test the measured $`d_\mathrm{eff}`$ lay 0.6–1.3% below the model (the script's port of
   `effectiveStemWidth`) in the 1.27 and 1.73 m range bins, for every rod diameter (−0.9 to −1.9 s.e.), and scattered
   on both sides of it elsewhere (within 2.1 s.e.). This is probably the way
   range was defined when extracting $`d_\mathrm{eff}`$ from the rays, not a library error, but it has to be resolved before the
   rod test is used quantitatively, as it is half the ±2% budget.

   **Table 18.** Accuracy needed for each beam parameter alone, as the range of assumed / true values that keeps the
   density within ±2% and ±5% (0–0.1 m, 300 stems/m²; thresholds absolute, true 0.05).

   | Beam, estimator | $`w_0`$, ±2% | $`w_0`$, ±5% | $`\gamma`$, ±2% | $`\gamma`$, ±5% | $`\tau`$, ±2% | $`\tau`$, ±5% |
   |---|---|---|---|---|---|---|
   | Default, ss | 0.945–1.050 | 0.870–1.133 | 0.884–1.107 | 0.735–1.288 | 0.044–0.056 | 0.036–0.065 |
   | Default, sw | 0.953–1.042 | 0.887–1.111 | 0.900–1.090 | 0.771–1.246 | 0.045–0.055 | 0.038–0.063 |
   | Default, ss, 1000/m² | 0.943–1.053 | 0.864–1.139 | 0.853–1.138 | 0.652–1.367 | 0.043–0.056 | 0.035–0.067 |
   | Narrow, ss | 0.925–1.068 | 0.826–1.181 | 0.587–1.412 | beyond 0.5–2 | 0.040–0.060 | 0.028–0.079 |
   | Wide, ss | 0.960–1.038 | 0.904–1.100 | 0.888–1.104 | 0.746–1.280 | 0.045–0.054 | 0.039–0.061 |

   Combined, the mean $`d_\mathrm{eff}`$ must be within ±2% at the survey's ranges for ±2% in the density: about ±0.45 mm on the
   23 mm $`d_\mathrm{eff}`$ of 1 cm stems at 2 m with the default beam. ss already reads 1–2% high at 0.1–0.2 m (Section 8.1), so
   the footprint should take only about ±1% of that: roughly $`w_0`$ to ±2.5%, $`\gamma`$ to ±5% and $`\tau`$ to ±6%, each of which alone
   moves the density by about 1%.

   *Return convention.* With `surface` ranging, grazing returns lie close to the plane through the rod's axis normal
   to the beam; with `outline` ranging they lie up to $`(d_\mathrm{eff} - d)/2`$ nearer the scanner. The difference (about 6–10 mm)
   is below single-shot range noise, so it needs many returns averaged; it matters, because `outline` data read 1–2%
   high (3–4% more where stems end) whichever convention the estimator is told (Section 7.9). The spread of returns
   also measures the range noise, and a pair of rods a few diameters apart shows how close stems shadow each other.
   Scan rods of different heights too, to see how real tops return.
2. **Known-$`N d`$ array (end-to-end test).** Build arrays of dowels of known diameter and positions (random,
   non-overlapping; at least two densities, e.g. a few hundred and about 1000 per m²; uniform and varied heights;
   and rows imitating cable roots) over a few m² of flat ground, and scan them with the survey's hand-held SLAM
   workflow (same walking speed, line spacing and ground classification). Compare cell and layer estimates with
   the known $`N d`$, with the hybrid and ss. This tests range noise, pose error, ground classification, stem tops and
   the return convention together, with the beam parameters fixed from step 1.
3. **Field quadrats.** In the survey area, count stems and measure their diameters at several heights (e.g. at the
   middle of each 0.1 m layer, to capture taper) in 1 × 1 m quadrats covering the range of densities and
   arrangements (random, patchy, cable roots). Record the minimum spacing between stems (for `--min_gap`), the stem
   heights, and for cable roots the lateral spread of the stems about the root lines and their density along them.
   Use the measured mean diameter per layer as `--stem_diameter`, and compare the estimates with the quadrat $`N d`$,
   per layer. Keep the quadrats used for checking separate from any used to tune parameters: calibrating $`\tau`$ on
   quadrats absorbs the diameter error (and, without the correction, the bias of Section 8.1), and the fitted value
   is not the physical threshold: it spread by 4% between densities in simulation, and must not be reused for
   another scanner or stem size (Section 7.10).
4. **Internal checks on survey data.** Compare estimates from near and far rays (or from different walking lines)
   over the same cells, which tests the footprint model; check the distribution of `observed_fraction`
   (Section 8.3), the `patchiness_ratio` (Section 6) and the number of voxels at the solid fraction limit (printed).
   Check the trajectory's consistency between passes (e.g. a trunk's position from different lines), as pose error
   above about 2 cm matters.

### 9.3 What a field campaign must measure, and how well

Ranked by the size of the error a plausible mistake causes (ss, 1 cm stems, default beam, 300 stems/m², 0–0.1 m,
unless stated; the accuracy for ±2% is interpolated from the synthetic sensitivities). The beam rows (1, 3, 4, 8,
10, 12) come from block L (Section 7.10), and the rest from the earlier stages.

| Rank | Quantity | Effect in simulation | Needed for ±2% | Measured by (Section 9.2) |
|---|---|---|---|---|
| 1 | Footprint correction itself | none: ×1.45–3.3 (narrow to wide beam) | always give `--footprint` | step 1 |
| 2 | Arrangement: lateral spread of stems about cable roots, their density along them, patch sizes | ss: patches −2 to −8%, cable roots −8 to −24% (150/m², 1.5–6 cm spread) and −29 to −64% (500/m²); hybrid at 3 cm: −2 to −6% (150/m²), −29 to −57% (500/m²) (Sections 8.6, 8.1.3) | decides the configuration; with ss, the spread to about 0.5 cm | step 3 |
| 3 | Footprint datasheet definitions | FWHM exit diameter read as 1/e²: +18% (reverse −21%); half-angle divergence read as full: +10% (reverse −15%) | the 1/e² diameter and full angle, or the rod test | step 1 |
| 4 | Detection threshold $`\tau`$ | 0.02–0.15 when it is 0.05: −12 / +26% (−11 to −12 / +22 to +27% over 100–1000/m², Table 4) | 0.044–0.056 alone (±12%); ±6% within the beam's ±1% share | step 1 |
| 5 | Minimum gap between stems | assuming no gap when stems are 1 diameter apart: +10% (300/m²), +34% (1000/m²) (Table 6) | about 0.5, 0.2 and 0.07 diameters at 100, 300 and 1000 stems/m² | step 3 |
| 6 | Overlapping (Poisson) stems | −13% at 1000/m², −3% at 300/m² (Section 7.9) | — (a model check: real stems don't overlap) | — |
| 7 | Stem heights and tops | +5 to +20% at 0.2–0.5 m from the tops (Section 8.2) | not reached for layers holding tops; record the height distribution | steps 1–3 |
| 8 | Exit diameter $`w_0`$ | ±25%: −9 / +11% | ±5% alone; ±2.5% within the beam's share | step 1 |
| 9 | Trajectory (pose) error | ss: ≤ 2 cm 0.994–1.007, 5 cm −3 to −11% (Section 7.9); hybrid: 2 cm −5.6 / −6.8% (Section 8.1.3) | below about 2 cm (ss), about 1 cm (hybrid) | step 4; the SLAM system |
| 10 | Beam parameters matched at one range only | −13 to +18% at 3.5 m from the line; 1.2–5.2% over a survey with lines 1 m apart (1.9–8.6% at 1000/m²) | calibrate at 1–5 m | step 1 |
| 11 | Mean stem diameter per layer | ±20%: −10 / +8%; ±10%: −6.6 / +6.3% (0.5 cm stems), −4.9 / +4.3% (1 cm), −1.6 / +1.1% (2 cm) (Section 8.4) | about ±3% (0.5 cm), ±4% (1 cm, 0.4 mm), ±15% (2 cm) | step 3 |
| 12 | Divergence $`\gamma`$ | ±25%: −4 / +5% | ±11% alone; ±5% within the beam's share | step 1 |
| 13 | Return convention | `outline` data: +1.2–1.9%, +3–4% more where stems end (Section 7.9) | which of the two | step 1 |
| 14 | Range noise | ss: ≤ 2 cm changes it by at most 1.3%, 3 cm by −0.7 / −1.7% (300/m²) and −1.1 / −3.3% (1000/m²) (Section 7.9); hybrid: 1 cm −1%, 2 cm −3 to −3.4% (Section 8.1.3) | below about 2 cm (ss), about 1 cm (hybrid) | step 1 |
| 15 | Diameter distribution | log-normal CV 30%: +0.3–0.7%; 6/14 mm mixture: +0.7–1.3% (Section 7.5) | the mean matters; note strong bimodality | step 3 |
| 16 | Scan geometry (line spacing 0.5–3 m, scanner height 0.8–2 m, ray angles) | ss within ±2% at 0–0.2 m in every geometry of block A (Section 7.9) | — | planning |
| — | Observed fraction (scan design) | below 0.9, cells read 3–22% low (Section 8.3) | lines about 1 m apart at ≥ 1000 stems/m²; `--min_observed_fraction 0.9` | planning; step 4 |

So for uniform stands, where the arrangement does not apply, the beam footprint without a rod calibration is the
largest error source: in the lowest 0.2 m, a datasheet misread or an unknown threshold alone exceeds every other
plausible mistake except the minimum gap at high density. Its parameters act only through the mean $`d_\mathrm{eff}`$ (Eq. 15), so the rod test needs to
fix that mean to about ±1% over the survey's ranges (Table 18). The threshold and the diameter both enter through
$`d_\mathrm{eff}`$, so a threshold calibrated on a counted plot absorbs a diameter error (Section 8.4); calibrate the beam from
the rod test and measure the diameter directly. The minimum gap matters only at high densities, and the arrangement
decides which configuration to use rather than entering as a parameter.

## 10. Recommended settings and workflow

**Settings.** `--density_model vertical_cylinder --estimator pooled --no_flat_top` (the `raydensity` defaults), 5 cm
voxels (`--pixel_width 0.05` / `--voxel 0.05`), 1 × 1 × 0.1 m cells, `--stem_diameter` the measured mean diameter at
the heights analysed, `--footprint w0,divergence,threshold` from the rod test (Section 9.2), `--return_range
surface` unless the rod test shows otherwise, and `--min_gap` from the field. Use `density` as the estimate, in one
of the configurations below. Always use only cells with `observed_fraction` of at least 0.9:
`--min_observed_fraction 0.9` blanks the others (it is off by default; over uneven ground use heights above the
ground, Section 8.3). Treat a `patchiness_ratio` above about 1.3 in the lowest 0.2 m (of a run at the default
window) as patchiness or clustering within the cell (`pooled_density` is then unusable), and treat layers holding
stem tops as uncertain by 5–20% (Section 8.2).

**Configuration.**

* **The hybrid (`--self_shadow_correction --hybrid`): recommended**, for uniform and clustered stands alike, where
  the range and pose errors are about 1 cm or less. Each cell takes the single-voxel estimate (sw) where it observes
  at least 0.9 of the cell's voxels, and the first-order correction at the default window (ss) otherwise, chosen
  from ray counts only (Section 8.1.3). Its errors in simulation (Table 17; `density` / truth, 0–0.1 m):

  | Case | Hybrid |
  |---|---|
  | Uniform stands, where ss passes ±2% | within ±2%; worst 1.019 ± 0.004 (0.5 cm stems, ×0.25 rays) |
  | Patches | 0.99–1.01 |
  | Cable roots, 150 stems/m² | 0.94–0.98 |
  | Cable roots, 500 stems/m² | 0.43–0.71 (0.47–0.75 over the cells at least 0.9 observed) |
  | Dense, poorly penetrated stands | the same as ss, including its shortfall where few cells are well observed |
  | Range noise 1 cm / 2 cm | 0.990 / 0.966–0.970 (ss passes) |
  | Pose noise 2 cm | 0.93–0.94 (ss passes) |

  `hybrid_weight` shows which cells took sw. The hybrid has not been run with a mismatched footprint, with
  `outline` returns or at voxel sizes other than 5 cm.
* **ss (`--self_shadow_correction`, default pool radius): the robust fallback**, where range or pose error exceeds
  about 1–2 cm, or cannot be bounded. It met ±2% in every scan geometry of the validation matrix, with range and pose
  error up to 2 cm, and from 100 to 1000 stems/m² where the cells were at least 0.9 observed (Section 7.9). It reads
  low in clustered stands (patches 0.92–0.98, cable roots 0.76–0.92); a stand is uniform enough for it where the
  `patchiness_ratio` of a run at the default window without the correction is close to that of random stems
  (1.03–1.09 in the lowest 0.2 m).
* **The default (no correction)** only with awareness of the cancellation: it reads 4–9% high for uniform stems, and
  its near-unity values in patches, and its −13% along cable roots at a 3 cm spread (−2 to −19% for 1.5–6 cm), rest
  on the self-shadowing bias cancelling the low bias of uneven windows (Section 8.1), which need not hold elsewhere.
  Below an observed fraction of 0.9 the same cancellation hides a low bias; mask those cells in every configuration.

sw on its own (the single-voxel options of Section 8.1.2) is superseded by the hybrid, which equals it where it is
well observed and falls back to ss elsewhere. Don't choose between configurations cell by cell by their values,
which is right only by the cancellation of Section 8.1: the hybrid chooses by ray counts, not values. Report the
configuration used: `raydensity` records each option in the store's `.zattrs` (`"hybrid": true`,
`"min_observed_fraction": 0.9`).

For example, for a plot:

```console
rayrender cloud.ply top density --density_model vertical_cylinder --estimator pooled --no_flat_top \
  --stem_diameter 0.01 --footprint 0.01,0.003,0.05 --pixel_width 0.05 --self_shadow_correction --hybrid \
  --min_observed_fraction 0.9 --cells cells.csv --cell_size 1,1,0.1
```

and for a tile of a domain, `raydensity ... --stem_diameter 0.01 --footprint 0.01,0.003,0.05
--self_shadow_correction --hybrid --min_observed_fraction 0.9` (with the measured diameter and the calibrated
footprint in place of these values).

**Workflow.**

1. Import the cloud with its trajectory (`rayimport`), so each point has its sensor origin.
2. Extract a ground mesh and set ground returns to zero alpha, so they count as unobstructed rays:

   ```console
   raydecimate cloud.ply 3 cm                                  # cloud_decimated.ply, for a lighter mesh
   rayextract terrain cloud_decimated.ply                      # cloud_decimated_mesh.ply
   raysplit cloud.ply cloud_decimated_mesh.ply distance 0.005  # cloud_inside.ply: ends > 5 mm above the mesh
   raycolour cloud_outside.ply alpha 0                         # cloud_outside_coloured.ply: the rest, as misses
   raycombine all cloud_inside.ply cloud_outside_coloured.ply --output cloud_no_ground.ply
   ```

   `raysplit` puts every ray whose end is more than the distance above the mesh in `cloud_inside.ply`, and the
   rest (ends below the mesh, within the distance of it, or beyond its extent) in `cloud_outside.ply`. This was
   checked on one synthetic scene only (`raycreate pneumatophores 101 --stems_per_m2 300 --footprint
   0.01,0.003,0.05 --plot_width 4 --rays 900000 --ground_returns`, flat noise-free ground): it removed every ground
   return, and the layers above 0.1 m were unchanged, but stem returns near the ground are removed with them, so the
   lowest 0.1 m read 0.954 of the truth against 1.062 with perfectly labelled ground (8.8% of its stem returns
   lost). A larger distance loses more (0.912 at 1 cm, 0.818 at 2 cm); keeping the ground returns read 2.84.
   Decimating first matters for speed: with the mesh of the full cloud (308 000 triangles, against 155 000),
   `raysplit` had not finished on the 0.9 M rays after 4 minutes, against 32 s with the mesh of the 3 cm decimated
   cloud. On real ground (rough, sloping, wet, with range noise) the distance must exceed the noise, and the loss at
   the stem bases will be larger: check it on the known-$`N d`$ array of Section 9.2.
3. Detect trunks above the pneumatophore layer (a circle fit per cluster, or `rayextract trunks`), extrapolate them
   to the ground, and set returns within each trunk's outline plus a margin larger than the range noise to zero
   alpha.
4. Split into tiles and estimate (`scripts/raydensity_tiles.sh`, which runs `raysplit ... tiles` and `raydensity`
   per tile, with `--ground` for heights above the mesh), writing one Zarr store; see
   [`large_scale_density.md`](large_scale_density.md). For a plot, `rayrender ... density --cells` does the same.
5. Check `observed_fraction` (mask cells below 0.9, e.g. with `--min_observed_fraction 0.9`), `hybrid_weight`,
   `patchiness_ratio` and the saturation warnings, and check that the range and pose errors allow the hybrid (else
   use ss), as above.
6. Add each trunk's frontal area D·Δz/V to its cell(s), with D corrected for the footprint if needed, and multiply
   the pneumatophore density by the cell's free area fraction.
7. Validate against quadrats (Section 9.2).

## 11. Equations and code

| Eq. | Content | Code |
|---|---|---|
| (2) | $`G = \sin\theta`$ | `raylib/raydensitygrid.h`: `gFunction` |
| (3), (8) | Beer–Lambert hits / weighted path, per voxel and pooled | `raylib/raydensitygrid.h`: `DensityGrid::operator()`, `Voxel::density`, `DensityGrid::voxelDensity`; `raylib/raydensitygrid.cpp`: `calculateDensities`, `applyOverlapCorrection` |
| (4)–(6) | Gaussian footprint, effective width | `raylib/raystemmodel.cpp`: `stripInterceptFraction`, `effectiveStemWidth`, `StemWidthModel` |
| (7) | Crossing extension | `raylib/raydensitygrid.cpp`: `DensityGrid::calculateDensities` |
| (9) | Response R and its inversion | `raylib/raystemmodel.cpp`: `freePathResponse`, `StemResponseTable`, `invertStemResponse`; `raylib/raystemmodel.h`: `invertDensityResponse`; `raylib/raystemresponsetable.h` |
| (10) | Cell outputs | `raylib/raydensitygrid.h`: `DensityCell`; `raylib/raydensitygrid.cpp`: `aggregateDensityCells`; `raylib/rayrenderer.cpp`: `writeDensityCells`, `writeDensityProfile`; `raylib/raydensity.cpp`: `estimateDensityTile` |
| (11), (12) | Mean exit distance; self-shadow correction | `raylib/raystemmodel.cpp`: `meanExitDistance`; `raylib/raydensitygrid.cpp`: `applyOverlapCorrection`; `raylib/raydensitygrid.h`: `DensityOptions::self_shadow_correction` |
| (13), (14) | Second-order self-shadow factor; occlusion neighbourhood; finite-ray correction | `raylib/raystemmodel.cpp`: `selfShadowFactor`; `raylib/raydensitygrid.cpp`: `applyOverlapCorrection`, `enableRatioBiasCorrection`; `raylib/raydensitygrid.h`: `DensityOptions::self_shadow_second_order`, `occlusion_radius`, `ratio_bias_correction` |
| — | Neighbour top-up, observed voxels | `raylib/raydensitygrid.cpp`: `addNeighbourPriors`, `estimateDensityGrid` (`DensityOptions::pooled_neighbour_priors`) |
| — | Hybrid estimate (Section 8.1.3) | `raylib/raydensitygrid.cpp`: `estimateDensityGrid` (the single-voxel copy), `DensityGrid::combineHybrid`, `aggregateDensityCells` (the per-cell choice); `raylib/raydensitygrid.h`: `DensityOptions::hybrid`, `hybrid_min_observed`, `DensityGrid::voxelDensity(size_t)`, `DensityCell::hybridWeight` |
| — | Observed-fraction mask | `raylib/raydensitygrid.h`: `DensityCell::mask`; `raylib/rayrenderer.cpp`: `writeDensityCells`, `writeDensityProfile`; `raylib/raydensity.cpp`: `estimateDensityTile` |
| — | Synthetic scenes, diameter distributions | `raylib/raypneumatophoregen.cpp`: `PneumatophoreGen::generate`, `castRay` |

## 12. Reproducibility

The scripts are in `scripts/`. The Python ones need the raycloudtools tools on the PATH (or a build's tool
directories given with `--bin`) and numpy, and write into the current directory or `--out`; the C++ programs build
against raylib as described in their headers. Unless stated, scenes are `raycreate pneumatophores <seed> <scene
options> --footprint 0.01,0.003,0.05 --truth_cells 1,1,0.1`, and estimates `rayrender pneumatophores.ply top density
--density_model vertical_cylinder --estimator pooled --no_flat_top --pixel_width 0.05 --min_bound -3,-3,0
--max_bound 3,3,<0.1 × layers> --cells cells.csv --footprint 0.01,0.003,0.05 --stem_diameter 0.01 <variant
options>`, compared with `pneumatophores_truth_cells.csv` by summing the 36 cells of each 0.1 m layer.
`scripts/pneumatophore_report_study.py` does this for lists of scenes and variants, generating each scene once per
seed and estimating it with every variant; for example, Table 2:

```console
scripts/pneumatophore_report_study.py --out table2 --layers 2 --seeds 4 \
  --scene "random 30/m2=--stems_per_m2 30" --scene "random 100/m2=--stems_per_m2 100" \
  --scene "random 300/m2=--stems_per_m2 300" --scene "random 600/m2=--stems_per_m2 600" \
  --scene "random 1000/m2=--stems_per_m2 1000" --scene "random 1500/m2=--stems_per_m2 1500" \
  --variant "pooled corrected=" --variant "pooled uncorrected=--no_overlap_correction" \
  --variant "crossing=--estimator crossing" --variant "beer_lambert=--estimator beer_lambert --no_overlap_correction"
```

In the table, "study" is `scripts/pneumatophore_report_study.py` and "none" the variant with no extra options.

| Result | Seeds | Produced by |
|---|---|---|
| Table 1 | 101–108 | `scripts/pneumatophore_seed_study.py --seeds 8` (and `--extra "--pool_radius 0.1"`) |
| Table 2, Section 5 (correction effect) | 101–104 | study `--layers 2`, as above |
| Table 3 | 101–104 (0.3 m patches: Table 1 run) | study; scenes `--stems_per_m2 300 --patchy 0.5,0.3,0`, `--stems_per_m2 300 --patchy 0.2,0.3,0`, `--stems_per_m2 150 --clustered --cluster_spread 0.03`; variants none, `--pool_radius 0.1`, `--no_overlap_correction`. Profile layers: `scripts/pneumatophore_validation_plots.py` (seed 11) |
| Section 7.4, angles | 101–104 | study `--layers 4`; scenes 300 and 1000 stems/m²; variants none, `--min_sin_zenith 0.9` |
| Section 7.4, ray count | 101 | study `--seeds 1`; scenes `--stems_per_m2 300` with `--rays 500000`, the default and `--rays 8000000` |
| Table 4 | 101–104 | study `--layers 2`; scenes 100, 300, 1000 stems/m²; variants `--footprint 0.01,0.003,0.02 / 0.05 / 0.15` |
| Table 5 | 101 | study `--seeds 1`; variants `--stem_diameter 0.008 / 0.012` |
| Table 6 | 101–104 | study `--layers 2`; scenes 100, 300, 1000 stems/m² with `--min_gap 0.5` or `--min_gap 1`; variants `--min_gap 0` and the true gap |
| Return convention, Section 7.6 | 101–104 | study `--layers 2`; scenes and variants `--return_range surface / outline` |
| Table 7 | 101–104 | study; variants `--pixel_width 0.025 / 0.05 / 0.1` |
| Beam, Section 7.6 | 101–104 | study; variant `--footprint 0,0,0.05` (with and without `--no_overlap_correction`); scenes and variant `--footprint 0.001,0.0001,0.05` |
| Section 7.7, tiling and slope | 101–104 (tiling: 101) | `scripts/pneumatophore_slope_tiles_study.py --seeds 4` (it runs `scripts/raydensity_tiles.sh` for the tiling) |
| Section 7.8 | 31–34 | `scripts/pneumatophore_trunk_study.py --seed <s>` (and `--pooled`) in `seed_<s>`, aggregated with `scripts/pneumatophore_trunk_aggregate.py`; trunks alone: `scripts/pneumatophore_trunk_only_study.py` |
| Section 8.2 table, uniform 0.6 m stems | 101 | study `--seeds 1 --layers 5`; scenes 300 and 1000 stems/m² with `--return_range surface / outline`, and with `--height 0.6,0` for stems of one height (see below) |
| Section 8.5, 5 cm layers | 101 | the Section 8.2 scenes, estimated as above with `--max_bound 3,3,0.5 --cell_size 6,6,0.05` |
| Section 5, response accuracy | 6 Monte Carlo seeds | `scripts/pneumatophore_response_noise.cpp` |
| Section 8.1 (Tables 8 and 9), Sections 6, 8.3 and 9.1 (stage-1 values) | 101–104 | `scripts/pneumatophore_self_shadow_study.py --sets density,tall,beam,window,diameter,patches,small_window,rays,spacing,diameter_error,lowobs,noise --seeds 101,102,103,104` (4 m plots; `--sets all` now also runs the K sets); quarter rays: `--sets rays --rays_per_m 14062` |
| Sections 7.9, 8.2–8.7 and 9.3 (validation matrix: Tables 11 and 12, the block G table, the observed-fraction bins) | 101–104; failing ss configurations also 105–108 | `scripts/pneumatophore_validation_matrix.py --blocks all --jobs 4 --out matrix`, then `--blocks all --rerun_failures 105,106,107,108 --out matrix` and `--blocks X --out matrix` (diagnostics). It writes `summary.csv`, `summary.txt`, `bins.csv` and `line_distance.csv`; the pass counts apply the decision rule of Section 7.9 to `summary.csv` (variants `base` and `ss`; block D at the true diameter, block H at matched conventions) |
| Section 8.1.2 (Tables 13 and 14, sw) | 101–104 | `scripts/pneumatophore_self_shadow_study.py --sets K2 --seeds 101,102,103,104 --jobs 2 --out small_window` (Tables 13 and 14; the ablations are its variants `sw priors`, `sw coupled`, `sw first order` and `sw no ratio`), `--sets K3,K4` (diameters, footprints, assumed diameter ±20%) and `--sets K5` (stems of one height, overlapping stems); `--sets K` for the voxel × pool radius grid |
| Section 8.1.3 (Tables 16 and 17, hybrid) | 101–104; the 0.5 cm stems at ×0.25 rays 101–108 | `scripts/pneumatophore_self_shadow_study.py --sets H,H2,H3 --seeds 101,102,103,104 --jobs 2 --out hybrid` (H: the block K scenes at ×0.25, ×0.5 and ×1 rays; H2: stem diameters and footprints; H3: stems of one height), and `--sets H2 --cases "d 0.005 m, 300/m2, rays x0.25" --seeds 101,102,103,104,105,106,107,108 --out hybrid_d005` (the worst case) |
| Section 8.1.3, block S (noise and scan geometry rows) | 101–104 | `scripts/pneumatophore_validation_matrix.py --blocks S --jobs 2 --out matrix_s` (S runs only when named; variants base, ss, sw and hybrid) |
| Section 8.1.3, per-voxel rules and the selection effect | 101–104 | voxel dumps of the set H scenes from a build with a temporary dump of both estimates per voxel, which is not in the repository, so these numbers can't be reproduced from it |
| Section 7.10 (Table 15, sensitivities, Eq. 15, diameter × threshold, range), Table 18 and the beam rows of Section 9.3 | 101–104 | `scripts/pneumatophore_validation_matrix.py --blocks L --jobs 2 --out footprint` (L runs only when named), then `scripts/pneumatophore_footprint_mismatch.py footprint`, which writes `footprint_tables.txt`, `footprint_sens.txt`, `footprint_deff.txt`, `footprint_grid.txt`, `footprint_range.txt` and `footprint_rows.csv` into `footprint/` |
| Section 9.2, rod test | 201–202 (synthetic rods) | `scripts/pneumatophore_rod_calibration.py design --out rod` (Fisher analysis: `rod_design.txt`) and `scripts/pneumatophore_rod_calibration.py simulate --out rod` (4 M rays per scene, needs `raycreate`; `rod_simulation.txt`) |
| Section 8.1.1 (Eqs. 13, 14) | 1–2 per configuration | `scripts/pneumatophore_self_shadow_window_mc.cpp` (key–value arguments, e.g. `N 300 d 0.01 deff 0.0225 w 0.05 conv 0 M 6e6 seed 1 D 0.08 L 3`), over $`\delta = d_\mathrm{eff}/w`$ in 0.1, 0.2, 0.3, 0.45, 0.6, 0.8, 1, $`N d_\mathrm{eff} w`$ in 0.03, 0.3, 1, 2 and $`d/d_\mathrm{eff}`$ in 1, 0.44, 0.3 ($`\phi \le 0.25`$), seeds 1 and 2; the moments of Eq. (13) by numerical integration over the square |
| Section 8.1, 2D Monte Carlo | 3 per configuration (6 for one) | `scripts/pneumatophore_self_shadow_mc2d.cpp` (key–value arguments, e.g. `N 300 d 0.01 deff 0.0234 D 0.22 w 0.05 k 1`) |
| Table 10 and the remedies of Section 8.3 | 101 (and 102 at 1000 and 1500 stems/m²) | `scripts/pneumatophore_low_observation_study.py` |
| Section 10, ground removal | 101 | the commands of step 2 on `raycreate pneumatophores 101 --stems_per_m2 300 --footprint 0.01,0.003,0.05 --plot_width 4 --rays 900000 --truth_cells 1,1,0.1 --ground_returns`, estimated over `--min_bound -2,-2,0 --max_bound 2,2,0.3`, against the same scene without `--ground_returns` |

The uniform 0.6 m scenes quoted in Section 8.2 were made with a small program before `raycreate` had `--height`.
`--height 0.6,0` gives the same stem positions and truth below 0.5 m, but not a bit-identical scan (at 1000
stems/m² the number of stem hits differs by 32 of 1.58 M), and layer values within 0.002 of those quoted. Plot totals are deterministic for a seed, so
the numbers are reproducible exactly with the same build. Run with one seed, the scripts reproduced their earlier
outputs exactly (the Table 2 cells at 300 stems/m², the slope and tiling comparison, trunks alone, the response
accuracy and the 300 stems/m² rows of Table 10's analysis).
