# Frontal area density of mangrove pneumatophores from hand-held lidar: method, synthetic validation and known biases

Technical report on the vertical-cylinder density estimate of raycloudtools (`--density_model vertical_cylinder` in
`rayrender` and `raydensity`). It describes the measurement model, the corrections, the cell outputs, a synthetic
validation, the known biases and what has not been tested, and gives a protocol for calibration and validation
with real data. Tiling, the Zarr output and large domains are described in
[`large_scale_density.md`](large_scale_density.md).

## Abstract

Hydrodynamic models of flow through mangrove pneumatophores need the frontal area per unit volume, a = N·d (stems
per m² of ground times stem diameter, m⁻¹), in cells of about 1 × 1 × 0.1 m over forests of about 300 × 300 m. We
estimate a from hand-held SLAM lidar ray clouds by treating each ray as a probe of a medium of vertical cylinders:
the interception rate per metre of ray is a·sin θ, which gives a G-weighted Beer–Lambert maximum-likelihood
estimate. Two corrections make it usable for 1 cm stems: (i) a Gaussian beam footprint model, which widens thin
stems to an effective width d_eff(r) of 2–3 cm and would otherwise overestimate a by 2.1–2.6 times; and (ii) an
occlusion (space-occupancy) correction for non-overlapping stems, from a deterministic Monte Carlo response table
R(φ, ρ, g), applied to hits and path pooled over a horizontal neighbourhood of voxels. Voxel estimates are averaged
into cells (`density`, the main output) and also pooled (`pooled_density`), and their ratio flags patchiness.

On synthetic scans of randomly placed stems (30–1000 stems/m², 5 cm voxels), `density` in the lowest 0.1 m reads
4–8% high, a geometric self-shadowing bias of about +0.47·d_eff/W (W the pooling window width) that does not depend
on the number of rays or the stem density and is under investigation; `pooled_density` is within 3%, but reads
about half the truth where stems are patchy within the cell. Stems in rows along cable roots read 13% low. The
estimate is about as sensitive to the assumed stem diameter (−11% / +9% for ±20% at 300 stems/m²) as to the beam
detection threshold (−11 to −12% / +22 to +27% for 0.4 / 3 times the true threshold), and to the minimum gap
between stems at high density (up to +34%). Layers above 0.1 m depend on how stem tops and partial returns are
treated by 5–15%, and poorly penetrated cells (observed fraction below about 0.9) read low. All results are
synthetic, from scenes generated with the same idealised beam and stem models that the corrections assume; the
method has not been checked against field data, and Section 9 gives a procedure for doing so.

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

$$a(z) = N(z)\,d, \qquad \phi(z) = N(z)\,\frac{\pi d^2}{4} = \frac{\pi d}{4}\,a(z), \tag{1}$$

where N(z) is the number of stems per m² of ground that reach height z, d the stem diameter, a the frontal area per
unit volume (m⁻¹) and φ the solid (volume) fraction. The drag force per unit volume of a canopy model is
proportional to a, and the frontal area index is ∫a dz. For stems placed independently at random (a Poisson
process, so stems may overlap), the probability that a horizontal line of length L misses every stem is
exp(−a L). Real stems don't overlap; for a non-overlapping (hard-core) arrangement the interception rate per unit
of free path is larger by R(φ, 1) (Section 5), so the gap fraction is approximately exp(−a R(φ, 1) L): the exponent
is 3% larger at 300 and 10% larger at 1000 stems/m² of 1 cm stems.

The target is a in cells of the hydrodynamic grid, e.g. 1 × 1 × 0.1 m (300 × 300 × 6 cells for a 300 m domain and
the lowest 0.6 m), from hand-held SLAM lidar. A ray cloud stores each return (or the end of a ray with no return)
with the sensor position at that time, so both where the beam stopped and the free space it crossed are known. The
estimator works on a voxel grid (5 cm by default) and aggregates the voxels into cells.

## 2. Measurement model

**Rays.** Each ray runs from the scanner position to its end point. A ray with nonzero alpha ended on a return
(a hit); a ray with zero alpha either had no return (its end is at the maximum range) or its return was removed
(e.g. classified as ground). Rays are clipped to the grid (`raylib/raydensitygrid.cpp`,
`DensityGrid::calculateDensities`), and the voxels they cross are walked (`DensityGrid::operator()` in
`raylib/raydensitygrid.h`).

**Interception by vertical cylinders.** A ray at zenith angle θ (measured from vertical) moves sin θ horizontally
per metre of ray. It intercepts a vertical cylinder of diameter d when the stem's centre lies within d/2 of the
ray's horizontal track, so for N stems per m², placed at random, the interception rate per metre of ray is

$$\lambda(\theta) = G(\theta)\,a, \qquad G(\theta) = \sin\theta. \tag{2}$$

This replaces the default G = 0.5 of raycloudtools, which applies to leaves with a spherical (isotropic) angle
distribution and gives leaf area per volume. For vertical stems it would bias the estimate by a factor 2 sin θ,
about 2 for the near-horizontal rays of a hand-held scanner (test `RayDensityVerticalCylinders`: at θ = 90° the
estimate was 1.003 of the truth with G = sin θ, and 2.005 with G = 0.5). The end caps (stem tops) are not part of
the model (Section 8.2).

**Likelihood.** If hits along a ray form a Poisson process with rate G a, a ray i that crosses a voxel with path
length L_i inside it contributes the log-likelihood h_i ln(G_i a) − a G_i L_i, where h_i is 1 if it ends with a hit
in the voxel and 0 otherwise. Maximising over a gives the G-weighted Beer–Lambert estimate

$$\hat a = \frac{H}{\sum_i G_i L_i}, \tag{3}$$

with H the number of hits in the voxel (or neighbourhood, Section 4). With the beam footprint (Section 3), G_i is
replaced by G_i ρ_i.

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
* *Ray angle selection.* Rays with sin θ below `--min_sin_zenith` (default 0.1) or above `--max_sin_zenith`
  (default 1) are ignored. Near-vertical rays carry almost no information on vertical stems, and selecting rays on
  direction alone does not bias the estimate.

## 3. Beam footprint and effective width

A lidar pulse has a finite, range-dependent spot. We model it as a Gaussian beam of 1/e² diameter

$$w(r) = w_0 + \gamma\,r, \tag{4}$$

with exit diameter w₀ and full divergence γ (both at the 1/e² level, from the scanner's specification), so its
standard deviation is σ = w/4. A vertical cylinder is seen from any horizontal direction as a strip of width d. If
the beam axis passes at horizontal offset x from the stem's centre, the fraction of the beam's power that the
stem intercepts is

$$F(x) = \tfrac12\left[\operatorname{erf}\!\left(\frac{x + d/2}{\sqrt2\,\sigma}\right) - \operatorname{erf}\!\left(\frac{x - d/2}{\sqrt2\,\sigma}\right)\right]. \tag{5}$$

A return is triggered when F ≥ τ, the detection threshold (a fraction of the beam power). The stem is therefore
detected over an effective width

$$d_\mathrm{eff}(r) = 2x^*, \quad F(x^*) = \tau, \qquad \rho(r) = d_\mathrm{eff}(r)/d, \tag{6}$$

found by bisection (`raylib/raystemmodel.cpp`, `stripInterceptFraction`, `effectiveStemWidth`), and zero when
F(0) < τ (a stem too thin ever to be detected). For 1 cm stems, a 10 mm, 3 mrad beam and τ = 0.05, d_eff is 2.07 cm
at 1 m and 3.05 cm at 5 m; the path-weighted mean ρ over the voxels of the synthetic scans was 2.13–2.42. The
interception rate becomes N d_eff sin θ = a ρ G, so each ray segment is weighted by G ρ(r), with r the range of the
segment's midpoint from the ray's original start (`StemWidthModel`, tabulated every 1 cm up to 100 m and on a
logarithmic grid beyond), and the estimate stays a = N d for the true mean diameter. Without this correction the
synthetic 1 cm stems were overestimated 2.1–2.6 times.

**Inputs.** w₀ and γ come from the scanner's beam specification. The threshold τ is not a published parameter: it
depends on the detector, the return mode, the target's reflectance and the range, and has to be calibrated
(Section 9). The diameter d is an input: it should be the mean diameter at the heights analysed (stems taper). The
result is about as sensitive to d as to τ (Section 7.5), because with a 13–25 mm footprint the hit rate is roughly
proportional to N (d + 1–2 cm), so the data mainly measure the stem count, and N d inherits about half of any
relative error in d.

**Return range convention.** When the beam axis misses a stem but enough power hits it, where along the beam is the
return reported? `--return_range` selects the model, in the generator and in the correction:

* `surface` (default): at the range of the physical stem surface (the nearest illuminated part of the stem), or of
  the beam axis's closest approach to the stem when the axis misses it. A time-of-flight scanner ranges to the
  illuminated surface, so this is the physically reasonable default.
* `outline`: where the beam axis enters the footprint-widened outline (diameter d_eff).

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

**`beer_lambert`** (the raycloudtools default): Eq. (3) per voxel, multiplied by the small-sample factor (n − 1)/n
of the original estimator, for the voxel's n rays (`DensityGrid::Voxel::density`). It read 1.24–1.31 of the truth
in the lowest 0.1 m and 1.31–1.39 at 0.1–0.2 m (cell means). The likely cause is the self-shadowing bias of
Section 8.1 at a window of one voxel (W = 5 cm, about +20%).

**`crossing`**: hit rays count their path to where they would have left the voxel, plus the mean chord of the
widened stem that remains beyond the return, divided by sin θ:

$$\Delta L = \frac{1}{\sin\theta}\times\begin{cases}\pi d_\mathrm{eff}/4 & \text{outline}\\ \pi (d_\mathrm{eff}^2 + d^2)/(8 d_\mathrm{eff}) & \text{surface}\end{cases} \tag{7}$$

(`DensityGrid::calculateDensities`). Without a footprint, density = hits / crossing length does not assume random
placement within the voxel. But it has no occlusion correction, and with a footprint the entries into a widened
stem that lie inside another widened stem are never observed. It read 1.05 at 30, 0.97 at 300, 0.82 at 1000 and
0.74 at 1500 stems/m², lower than with a near-perfect beam (1 mm, 0.1 mrad: 1.00 at 300 and 0.93 at 1000). It is not
recommended with a footprint.

**`pooled`** (recommended): the Beer–Lambert hits and weighted path are summed over a horizontal neighbourhood
𝒩(v) of each voxel v before taking the ratio,

$$m_v = \frac{\sum_{u\in\mathcal N(v)} H_u}{\sum_{u\in\mathcal N(v)} \sum_i G_i \rho_i L_i}, \qquad \bar\rho_v = \frac{\sum_{u}\sum_i \rho_i\,(G_i \rho_i L_i)}{\sum_u \sum_i G_i \rho_i L_i}, \tag{8}$$

and m_v is corrected for occlusion (Section 5) using the neighbourhood's mean widening ρ̄_v. The voxel takes the
corrected value a_v and keeps its own path length (so `pooled_density`, Section 6, weights it by its own path).
Voxels that no ray reached are not assigned a value. 𝒩(v) is the (2k+1) × (2k+1) voxels centred on v in its
horizontal layer, with k = round(`--pool_radius` / voxel width), so it spans ±(k + ½) voxels and its width is
W = (2k + 1) w: 0.15 m (±0.075 m) for the default 0.05 m at 5 cm voxels, 0.25 m for 0.1 m, and 0.30 m for 0.05 m at
10 cm voxels (the run prints the window). It changes in steps with the voxel size. The window is computed with
summed-area tables (`DensityGrid::applyOverlapCorrection`).

Why pool: the Beer–Lambert estimate is unbiased for randomly placed stems however few there are in a voxel, provided
hits and path are counted over a region, not averaged as ratios of single voxels; pooling moves the ratio from one
voxel's few stems to a neighbourhood's, and the occlusion correction handles non-overlapping stems. Pooling trades
resolution for bias: it smooths over the window, so patch edges spread by about the window, stems clustered more
tightly than the window read low (Section 8.6), and the remaining self-shadowing bias falls as the window grows
(Section 8.1). Pooled corrected values were 1.01–1.08 of the truth (cell means; 0.98–1.01 path-weighted).

With `beer_lambert` and the occlusion correction on (a stem diameter given), the neighbourhood ratio is used only to
compute the correction factor, which is applied to the voxel's own estimate.

## 5. Occlusion and space-occupancy correction

The pooled estimate measures hits per unit of free path. For randomly placed, overlapping discs (a Boolean model)
the free path to the next hit is exponential, and the hit rate per unit free path is exactly N d_eff, however much
the widened discs overlap. Real stems don't overlap. In a non-overlapping (hard-core) arrangement the
footprint-widened stems leave less free space than overlapping ones would, so the measured rate exceeds N d_eff.
With a perfect beam (d_eff = d) it is at least a/(1 − φ), because the stems' own cross-sections are excluded from
the free path, and a little more because rays enter the stem layer at random points, so sample the gaps between
stems in proportion to their length, and hard-core gaps are more regular than exponential.

**Response.** We define the response as the measured rate relative to N d_eff,

$$R(\phi, \rho, g) = \frac{\text{hits per unit free path}}{N\,d_\mathrm{eff}}, \qquad m = a\,R\!\left(\tfrac{\pi d}{4} a,\ \bar\rho,\ g\right), \tag{9}$$

where g is the minimum gap between stem surfaces in diameters (stem centres at least d(1 + g) apart; `--min_gap`).
R is computed by a 2D Monte Carlo simulation (`raylib/raystemmodel.cpp`, `freePathResponse`): in units of d,
1500 non-overlapping discs are placed by random sequential addition in a square sized for φ; 10 000 lines start
from uniformly random points of the free space in the central 80% of the square, in random directions, and run to
their first reported return or the edge of the square; R = Σ hits / (Σ path · N ρ), pooled over 64 independent
realisations. The free space follows the return convention: under `outline`, outside the widened outlines; under
`surface`, outside the physical stems, except inside a widened outline beyond its reported return (a ray there
would already have returned). A response is saturated (infinite) when under 1% of candidate start points are free.

**Table.** For g = 0, R is tabulated for both conventions over 13 solid fractions (0, 0.01, 0.02, 0.04, 0.07, 0.1,
0.13, 0.16, 0.2, 0.25, 0.3, 0.35, 0.4) and 17 width ratios (1 to 5 in steps of 0.25), embedded in the library
(`raylib/raystemresponsetable.h`, generated by `scripts/generate_stem_response_table.cpp`) and interpolated
bilinearly (`StemResponseTable::response`). For g > 0 the table is over the solid fraction of the exclusion discs,
φ(1 + g)², and its nodes are simulated when first needed. Examples: R(0.024, 2.3, 0) = 1.046 (300 stems/m² of 1 cm
stems), R(0.079, 2.3, 0) = 1.157 (1000 stems/m²) and R(0.2, 2.3, 0) = 1.51. With a perfect beam R(φ, 1, 0) is 1.03
and 1.10 at 300 and 1000 stems/m², so about a third of the excess above 1 comes from the footprint widening.

**Accuracy of the table.** Rerunning the simulation with 6 other seeds varied R by 0.2–0.5% (s.d.) at φ ≤ 0.2, and
by 0.6% (surface) and 0.9% (outline) at φ = 0.25. Interpolation differed from direct simulation by at most 0.2% for
φ ≤ 0.2, and 0.4% at φ = 0.3. The two conventions agree within 0.2% at φ ≤ 0.1 and 0.7% at φ = 0.2, and diverge
above: 1.6% at φ = 0.25 (ρ = 3) and 3.8% at φ = 0.3 (ρ = 2.1), surface lower.

**Inversion.** For each voxel, Eq. (9) is solved for a given the neighbourhood's m and ρ̄ by bracketing upwards
from m (doubling, so R < 1 would also be handled) and bisecting (`invertDensityResponse`, `invertStemResponse` in
`raylib/raystemmodel.h`/`.cpp`).

**Limits.** The exclusion-disc solid fraction is capped at 0.4 (below the 0.547 random sequential addition limit
for discs), so N d ≤ 1.6/(π d (1 + g)²): 51 m⁻¹ (5100 stems/m²) for 1 cm stems with no gap. Voxels beyond this are
set to the limit, and their number is printed. In the synthetic scans the correction lowered the estimates: the
uncorrected cell means were 1–6% higher at 30–300 stems/m², 11% at 600 and 18–30% at 1000–1500.

**Minimum gap.** With stems at least one diameter apart, assuming no gap overestimated the cell mean by 34% at 1000
stems/m², 10% at 300 and 4% at 100 (half a diameter apart: 15%, 5% and 2%), relative to assuming the right gap
(Table 6). g should be set from field observations of how closely stems pack.

**Determinism.** The simulation uses fixed seeds and its own uniform random numbers from `std::mt19937`, so R, the
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
voxel that no ray reached is not filled. A voxel that still has 2 or fewer rays after the top-up is unobserved. In
practice the unobserved voxels are those no ray reached (below the ground, beyond the scan, or completely hidden),
or whose neighbours are also nearly empty.

**Outputs.** Over the cell's n_obs observed voxels, with corrected voxel densities a_v and weighted path lengths
P_v = Σ_i G_i ρ_i L_i:

$$\texttt{density} = \frac{1}{n_\mathrm{obs}}\sum_v a_v, \qquad \texttt{pooled\_density} = \frac{\sum_v a_v P_v}{\sum_v P_v}, \qquad \texttt{patchiness\_ratio} = \frac{\texttt{density}}{\texttt{pooled\_density}}. \tag{10}$$

* `density` (main output): the mean of the observed voxel densities, each voxel weighted equally (by volume), so
  the volume average that a hydrodynamic cell needs. It follows patchy stems, but carries the self-shadowing bias
  of Section 8.1 (+4–8% for random stems in the lowest 0.1 m at the default settings).
* `pooled_density`: the cell's corrected hits over its weighted path, i.e. the path-weighted mean of the voxel
  densities. Rays penetrate dense parts less, so they get less weight, horizontally and vertically: it is within 3%
  only where the stems are uniform within the cell in both directions, about half the truth for patches within the
  cell, and low in layers where the density falls with height (Section 8.5).
* `patchiness_ratio`: in the lowest 0.2 m, 1.03–1.09 for random stems at the default pool (mostly the bias of
  `density`), 1.6–2.4 for patches and cable roots. Higher up, the vertical profile within the cell also raises it
  (up to 1.17 for random stems at 0.2–0.3 m).
* `hits` (Σ a_v P_v, the corrected hits), `weighted_path` (Σ P_v) and `observed_fraction` (n_obs over the cell's
  voxels). The profile (`--profile`) gives the same per horizontal layer of voxels, with `observed_voxels`,
  `total_voxels` and `mean_hits_per_voxel`.

A cell without observed voxels is NaN (with `observed_fraction` 0 in the Zarr store). Both densities use only the
observed voxels, so they describe the observed part of the cell: where much of the cell is unobserved, they read
low (Section 8.3).

**Coarser cells.** Cells can be aggregated further from the stored values: `pooled_density` from the summed `hits`
and `weighted_path`, and `density` as the mean weighted by the number of observed voxels (`observed_fraction` times
the voxels per cell). When the cell dimensions are multiples of the voxel width, every cell holds the same number
of voxels (also with `--ground`, as the grid spans the full relief), so weighting by `observed_fraction` is exact;
otherwise the voxel count per cell varies with the ground height, and isn't stored.

## 7. Synthetic validation

### 7.1 Generator and scenes

`raycreate pneumatophores <seed>` (`raylib/raypneumatophoregen.cpp`) generates vertical cylinders of one diameter
(1 cm) on flat ground at z = 0, with heights drawn from a normal distribution of mean 0.25 m and s.d. 0.08 m
(truncated to be positive; 3%, 24%, 47% and 24% of stems end in successive 0.1 m layers), and scans them with a
hand-held lidar walked along a lawnmower path:

* Arrangement: random, non-overlapping (random sequential addition with the exact number of stems; `--min_gap`
  for a minimum gap, `--overlapping` for Poisson positions); `--patchy radius,cover,matrix_ratio` (stems only in
  random circular patches); or `--clustered` (stems along cable roots radiating from trees, 0.04 trees/m², 12 roots
  of 6 m per tree, stems offset from the root lines with s.d. `--cluster_spread`).
* Scan: a 6 × 6 m plot (stems generated 8 m beyond it), scanner at 1.2 m height walking lines 1 m apart, 2 000 000
  rays (about 56 000 per m² of plot) uniform over zenith angles 0–101.5°, 30 m range for rays with no return; ground
  returns are given zero alpha (perfectly classified).
* Beam: Gaussian, 10 mm exit diameter, 3 mrad divergence, τ = 0.05, `surface` return range; stem tops are solid
  discs of the widened outline. No range noise, pose error or intensity.
* Truth: N d in each 1 × 1 × 0.1 m cell and each 1 cm layer, from the stems' positions and heights.

Scenes depend only on the seed and options, on any platform.

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

| Stems/m² | Observed fraction | `density` | `pooled_density` | Uncorrected `density` | `crossing` | `beer_lambert` (uncorrected) | Mean ρ |
|---|---|---|---|---|---|---|---|
| 30 | 1.00 | 1.080 (1.089) | 0.987 (0.999) | 1.093 | 1.055 | 1.299 | 2.42 |
| 100 | 1.00 | 1.075 (1.087) | 0.984 (0.996) | 1.101 | 1.032 | 1.299 | 2.35 |
| 300 | 1.00 | 1.063 (1.076) | 0.977 (0.987) | 1.123 | 0.974 | 1.313 | 2.25 |
| 600 | 1.00 | 1.061 (1.079) | 0.995 (0.989) | 1.175 | 0.898 | 1.296 | 2.19 |
| 1000 | 0.93 | 1.039 (1.076) | 1.008 (0.999) | 1.230 | 0.815 | 1.240 | 2.15 |
| 1500 | 0.68 | 1.011 (1.060) | 0.981 (1.014) | 1.317 | 0.739 | 1.246 | 2.13 |

`density` reads 6–8% high wherever the cells are fully observed, independently of the density (Section 8.1). The
lower values at 1000 and 1500 stems/m² are not better accuracy: there 7% and 32% of the lowest layer's voxels are
unobserved, and the low bias of poorly observed cells (Section 8.3) cancels part of the high bias.

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
clustered more tightly (3 cm) than the window (Section 8.6). This rests on one configuration (one density, one
spread, 4 scenes). `pooled_density` is about half the truth in every non-random arrangement, and lower still over a
whole 6 m profile layer (0.30 for 0.5 m patches and 0.34 for cable roots, one scene each).

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

**Table 4.** Scenes generated with τ = 0.05, estimated assuming τ (4 scenes, `density`, lowest 0.1 m).

| Stems/m² | Assume τ = 0.02 | τ = 0.05 (true) | τ = 0.15 |
|---|---|---|---|
| 100 | 0.943 (−12%) | 1.075 | 1.358 (+26%) |
| 300 | 0.939 (−12%) | 1.063 | 1.325 (+25%) |
| 1000 | 0.929 (−11%) | 1.039 | 1.265 (+22%) |

Relative changes are against the true-τ estimate; at 0.1–0.2 m and for `pooled_density` they are similar (up to
+27%).

**Table 5.** Scenes of 1 cm stems estimated assuming another diameter (one scene, seed 101, `density`, lowest
0.1 m).

| Stems/m² | Assume d = 0.8 cm | d = 1 cm (true) | d = 1.2 cm |
|---|---|---|---|
| 300 | 0.941 (−11%) | 1.056 | 1.147 (+9%) |
| 1000 | 0.954 (−8%) | 1.038 | 1.094 (+5%) |

The effect is smaller at high density, because the occlusion correction partly compensates. In the footprint model
alone, d_eff/d for an assumed 0.8 or 1.2 cm changes by −11% / +9% at 1 m range and −14% / +12% at 5 m relative to
1 cm, and a mixture of 6 and 14 mm stems estimated with the 10 mm mean differs by under 1% (model calculations, not
scene runs). So a spread of diameters about the right mean matters little, but the mean must be right.

### 7.6 Sensitivity to the minimum gap, the return convention, the voxel size and the beam

**Table 6.** Scenes with a minimum gap g between stems, estimated assuming g = 0 or the true g (4 scenes,
`density`, lowest 0.1 m).

| Stems/m² | g = 0.5, assume 0 | assume 0.5 | g = 1, assume 0 | assume 1 |
|---|---|---|---|---|
| 100 | 1.089 | 1.069 | 1.101 | 1.059 |
| 300 | 1.107 | 1.059 | 1.148 | 1.048 |
| 1000 | 1.176 | 1.025 | 1.344 | 1.004 |

**Return convention.** Scenes generated with `outline` read higher than `surface` scenes, with or without the
correction, and the difference grows with height (Section 8.2): 1–2% at 0–0.1 m and 2.5–4.5% at 0.1–0.2 m (4 scenes;
e.g. `density` 1.063 → 1.078 at 300 and 1.039 → 1.057 at 1000 stems/m²), and 5–10% at 0.2–0.5 m (one scene).
Estimating with the other convention changed the estimates by at most 0.2%.

**Table 7.** Voxel size (random stems, 4 scenes, default pool radius 0.05 m, so windows of W = 0.125, 0.15 and
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
of trunk frontal area, about 6–8% of its pneumatophore N d. 1 × 1 m cells, 4 scenes (2–5 trunks in each plot),
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

## 8. Known biases and limitations

### 8.1 Self-shadowing bias of `density` (open)

`density` reads high by 4–8% for random stems at the default settings, by the same amount at 30 and 600 stems/m²
(Table 2), and whatever the number of rays (Section 7.4). This is not the noise of ratios of few hits. The likely
cause is geometric self-shadowing: a stem removes the free path behind it within its own pooling window, so windows
that hold stems have less path and a higher hits/path than their true share, and the equal-weight mean of the
window ratios is high. A derivation (by an independent reviewer; the mean distance from a random point of a square
to its edge along a random direction is about 0.47 W) gives a relative bias of about

$$b \approx 0.47\,\frac{d_\mathrm{eff}}{W}, \tag{11}$$

independent of the stem density and of the scan density. Predicted (using the run's mean ρ̄, d_eff = ρ̄ · 1 cm)
against measured (`density` − 1, lowest 0.1 m):

| Case | W (m) | d_eff (cm) | Predicted | Measured | Observed fraction |
|---|---|---|---|---|---|
| 30 stems/m², 5 cm voxels | 0.15 | 2.42 | +7.6% | +8.0% | 1.00 |
| 100 stems/m² | 0.15 | 2.35 | +7.4% | +7.5% | 1.00 |
| 300 stems/m² | 0.15 | 2.25 | +7.1% | +6.3% | 1.00 |
| 600 stems/m² | 0.15 | 2.19 | +6.9% | +6.1% | 1.00 |
| 300 stems/m², `--pool_radius 0.1` | 0.25 | 2.25 | +4.2% | +3.3% | 1.00 |
| 100 stems/m², `--pool_radius 0.1` | 0.25 | 2.35 | +4.4% | +4.2% | 1.00 |
| 300 stems/m², 2.5 cm voxels | 0.125 | 2.25 | +8.5% | +7.1% | 0.99 |
| 300 stems/m², 10 cm voxels | 0.30 | 2.25 | +3.5% | +2.2% | 1.00 |
| 300 stems/m², near-perfect beam | 0.15 | 1.10 | +3.4% | +3.1% | 1.00 |
| 1000 stems/m², near-perfect beam | 0.15 | 1.10 | +3.4% | +3.0% | 1.00 |
| 1000 stems/m², 5 cm voxels | 0.15 | 2.15 | +6.7% | +3.9% | 0.93 |
| 1500 stems/m², 5 cm voxels | 0.15 | 2.13 | +6.7% | +1.1% | 0.68 |

Where the cells are fully observed, Eq. (11) predicts the measurements to within about 1.5 percentage points, and
over-predicts slightly at the narrowest and widest windows; where they are not, the low bias of Section 8.3 cancels part of it.
Consequences:

* More scanning will not reduce the bias, but it is predictable: it scales with d_eff/W, so it is larger for thicker
  stems, longer ranges (d_eff 3 cm at 5 m gives about +10% at W = 0.15 m) and smaller windows, and the +4–8% holds
  only for this configuration.
* The `patchiness_ratio` of about 1.09 for random stems is mostly this bias.
* A larger window (`--pool_radius 0.1`, or larger voxels) reduces it, but reads low where stems are clustered or
  patchy at a smaller scale than the window (Table 3). The choice of pool radius therefore trades this bias against
  the clustering bias; it is not a remedy for either.

This is a known, deterministic bias that is under active investigation; a correction is being developed separately.
Until then, treat `density` as reading high by about 0.47·d_eff/W, and note that Eq. (11) has been checked only for
the cases in the table above (one stem diameter, one scan geometry).

### 8.2 Stem tops and the return convention

The model has no stem tops. Steep rays end on stem tops, which counts as interception by the side of the stem, and
stems end in every layer above 0.1 m for the height distribution simulated. The generator's treatment of the tops
is itself a choice: the top is a solid disc of the widened outline, and a side hit is tested at the height where the
beam enters the widened outline but, under `surface`, reported lower.

The return convention interacts with the tops. Deep inside a stand of uniform height, the longer paths of `surface`
returns are balanced by as many extra returns carried down from the more strongly illuminated layer above, so R and
the estimates are convention-independent. Where stems end, there are no returns to carry down from above, so a
`surface` scene reads lower than an `outline` scene, by up to N d_eff R × 7 mm in the layer holding the tops (5% at
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

### 8.3 Poorly observed cells read low

Where rays penetrate poorly, the estimates read low, in both columns, probably because the observed voxels are the
less occluded ones (the mechanism has not been isolated). At 1500 stems/m², 32% of the lowest layer's voxels were
unobserved (mean `observed_fraction` 0.68), and `density` read 1.011, about 5% below the fully observed value; the
"+1%" is the cancellation of this with Section 8.1. With stems of uniform height 0.6 m at 1000 stems/m² (one scene),
`observed_fraction` was 0.28, 0.41, 0.58, 0.80 and 0.97 in the 0.1 m layers from the ground, and `density` read
0.911, 0.948, 0.990, 1.030 and 1.038 (`pooled_density` 0.873, 0.907, 0.955, 1.000, 0.990): 9–12% low in the two
lowest layers relative to the well observed ones. Treat cells with `observed_fraction` below about 0.9 as reading low
(about −5% at 0.7 and −10% or more at 0.3–0.4 in simulation). All validation scenes used one scan geometry (scanner
at 1.2 m, lines 1 m apart, about 56 000 rays/m²) and stems of 0.25 m mean height, which is a favourable case; field
scans with wider line spacing, longer ranges, or taller or denser stems will be more penetration-limited.

### 8.4 Inputs: diameter and threshold

The estimate is about as sensitive to the assumed mean diameter (Table 5: −11% / +9% for 0.8 / 1.2 cm at 300
stems/m², −8% / +5% at 1000) as to the detection threshold (Table 4: −11 to −12% / +22 to +27% for 0.02 / 0.15 when it
is 0.05). A threshold calibrated on a counted plot also absorbs any error in the diameter (and the bias of
Section 8.1), and such a calibration will not transfer to other ranges, stem sizes or densities (Section 9.2).

### 8.5 `pooled_density` is path-weighted horizontally and vertically

`pooled_density` weights each voxel by its path, and rays penetrate dense parts less, so it reads low wherever the
density varies within the cell: about half the truth for patches (Table 3), and also in the vertical, where the
upper voxels of a layer in which the density falls with height carry more path and less density. This is why it
read 0.95 at 0.2–0.3 m at 1000 stems/m² (Table 1): per 5 cm layer it read 1.043 (0.20–0.25 m, weighted path 23 400)
and 1.056 (0.25–0.30 m, 54 500), with no dip. It is within 3% only where the stems are uniform within the cell
horizontally and vertically, and it also carries the stem-top bias: at 300 stems/m² it read 1.075 at 0.3–0.4 m
(Section 7.4) and, per 5 cm layer, rose steadily from 0.966 at the ground to 1.224 at 0.45–0.50 m (one scene).

### 8.6 Clustering at scales below the pooling window

Stems clustered more tightly than the window read low, as any free-path estimate does when pooled over uneven
cover: cable roots with a 3 cm lateral spread read 0.87 (default window, 0.15 m) and 0.77 (0.25 m window), with a
`patchiness_ratio` of 1.6–1.8 in the lowest 0.2 m (Table 3). This is the normal arrangement of pneumatophores, so
it matters; only one configuration has been simulated.

### 8.7 Summary of residual biases (synthetic, default settings, 5 cm voxels, 1 × 1 × 0.1 m cells)

| Bias | Size | Applies to | Section |
|---|---|---|---|
| Self-shadowing (open) | +4–8% (≈ +0.47 d_eff/W) | `density`, all arrangements | 8.1 |
| Stem tops and return convention | 5–15% above 0.1 m, growing with height | both | 8.2 |
| Poor penetration | about −5% at observed fraction 0.7, −10% or more at 0.3–0.4 | both | 8.3 |
| Diameter error ±20% | −11% / +9% (300/m²), −8% / +5% (1000/m²) | both | 8.4 |
| Threshold error (0.4×, 3×) | −11 to −12% / +22 to +27% | both | 8.4 |
| Minimum gap ignored (g = 1) | +4% (100/m²) to +34% (1000/m²) | both | 5 |
| Patchiness within the cell | about −50% | `pooled_density` | 8.5 |
| Vertical path weighting | e.g. −5% at 0.2–0.3 m (1000/m²) | `pooled_density` | 8.5 |
| Clustering below the window (cable roots) | −13% (default), −23% (`--pool_radius 0.1`) | `density` | 8.6 |

## 9. What was not simulated, and a calibration and validation protocol

### 9.1 Not simulated

The synthetic validation uses the same beam, stem and return models as the corrections, so it tests the estimator's
consistency, not the physical models. In particular, it did not include:

* Range noise and SLAM pose error (typically 1–3 cm, larger than the 1 cm stems, the 6–7 mm convention difference
  and the 1.5 cm trunk removal margin). The generator's rays are exact.
* Mixed stem and ground returns at the stem base, and imperfect ground classification: the ground is flat (except
  for the slope test) and perfectly labelled, yet the main layer is 0–0.1 m.
* Stem taper, lean, curvature, branching, and a distribution of diameters (only one diameter was generated; the
  diameter mixture result of Section 7.5 is a model calculation).
* Reflectance and incidence dependence of the threshold, wet or muddy stems, and the scanner's actual return
  detection (multiple returns, return mode).
* More than one scan geometry, height distribution, or scanner (Section 8.3).
* Real stem tops (the generator's tops are a modelling choice, Section 8.2).
* Water, debris, seedlings, crabs and other non-stem returns in the stem layer.

None of the method has been checked against field data.

### 9.2 Calibration and validation protocol

The aim is to calibrate the beam model from the physics (not from the density result) and then validate the whole
chain on independent data of known N d.

1. **Beam footprint, threshold and return convention (static rod test).** Scan an isolated vertical rod or dowel of
   known diameter, similar to the stems in size and surface (e.g. 1 cm, painted or wetted to match their
   reflectance), from a static scanner at several ranges spanning those of the survey (e.g. 1–5 m), and a pair of
   rods a few diameters apart. The lateral spread of the returns across the rod, against range, gives d_eff(r),
   and fitting Eq. (5)–(6) with the specified w₀ and γ gives τ (or fit all three). For the return convention: with
   `surface` ranging, grazing returns lie close to the plane through the rod's axis normal to the beam; with
   `outline` ranging they lie up to (d_eff − d)/2 nearer the scanner. The difference (about 6–10 mm) is below
   single-shot range noise, so it needs many returns averaged. The spread of returns also measures the range noise.
2. **Known-N d array (end-to-end test).** Build arrays of dowels of known diameter and positions (random,
   non-overlapping; at least two densities, e.g. a few hundred and about 1000 per m²; uniform and varied heights;
   and rows imitating cable roots) over a few m² of flat ground, and scan them with the survey's hand-held SLAM
   workflow (same walking speed, line spacing and ground classification). Compare cell and layer estimates with
   the known N d. This tests range noise, pose error, ground classification, stem tops and the return convention
   together, with the beam parameters fixed from step 1.
3. **Field quadrats.** In the survey area, count stems and measure their diameters at several heights (e.g. at the
   middle of each 0.1 m layer, to capture taper) in 1 × 1 m quadrats covering the range of densities and
   arrangements (random, patchy, cable roots), and record the minimum spacing between stems (for `--min_gap`).
   Use the measured mean diameter per layer as `--stem_diameter`, and compare the estimates with the quadrat N d,
   per layer. Keep the quadrats used for checking separate from any used to tune parameters: calibrating τ on
   quadrats absorbs the diameter error and the bias of Section 8.1, which will not transfer to other ranges,
   densities or stem sizes.
4. **Internal checks on survey data.** Compare estimates from near and far rays (or from different walking lines)
   over the same cells, which tests the footprint model; check the distribution of `observed_fraction`
   (Section 8.3), the `patchiness_ratio` (Section 6) and the number of voxels at the solid fraction limit (printed).

## 10. Recommended settings and workflow

**Settings.** `--density_model vertical_cylinder --estimator pooled --no_flat_top` (the `raydensity` defaults), 5 cm
voxels (`--pixel_width 0.05` / `--voxel 0.05`), 1 × 1 × 0.1 m cells, the default `--pool_radius 0.05`,
`--stem_diameter` the measured mean diameter at the heights analysed, `--footprint w0,divergence,threshold` with the
threshold calibrated (Section 9.2), `--return_range surface`, and `--min_gap` from the field. Use `density` as the
estimate, keeping in mind its high bias (Section 8.1); use cells with `observed_fraction` of at least about 0.9;
treat a `patchiness_ratio` above about 1.3 in the lowest 0.2 m as patchiness within the cell (`pooled_density` is
then unusable, and `density` may read low if stems are clustered below the window); treat layers near the stem tops
as convention-dependent.

**Workflow.**

1. Import the cloud with its trajectory (`rayimport`), so each point has its sensor origin.
2. Extract a ground mesh (`rayextract terrain`) and set ground returns to zero alpha, so they count as unobstructed
   rays (for example by splitting returns near the mesh with `raysplit`, setting their alpha to zero with
   `raycolour ... alpha`, and recombining with `raycombine`; this step is not covered by the simulation).
3. Detect trunks above the pneumatophore layer (a circle fit per cluster, or `rayextract trunks`), extrapolate them
   to the ground, and set returns within each trunk's outline plus a margin larger than the range noise to zero
   alpha.
4. Split into tiles and estimate (`scripts/raydensity_tiles.sh`, which runs `raysplit ... tiles` and `raydensity`
   per tile, with `--ground` for heights above the mesh), writing one Zarr store; see
   [`large_scale_density.md`](large_scale_density.md). For a plot, `rayrender ... density --cells` does the same.
5. Check `observed_fraction`, `patchiness_ratio` and the saturation warnings.
6. Add each trunk's frontal area D·Δz/V to its cell(s), with D corrected for the footprint if needed, and multiply
   the pneumatophore density by the cell's free area fraction.
7. Validate against quadrats (Section 9.2).

## 11. Equations and code

| Eq. | Content | Code |
|---|---|---|
| (2) | G = sin θ | `raylib/raydensitygrid.h`: `gFunction` |
| (3), (8) | Beer–Lambert hits / weighted path, per voxel and pooled | `raylib/raydensitygrid.h`: `DensityGrid::operator()`, `Voxel::density`, `DensityGrid::voxelDensity`; `raylib/raydensitygrid.cpp`: `calculateDensities`, `applyOverlapCorrection` |
| (4)–(6) | Gaussian footprint, effective width | `raylib/raystemmodel.cpp`: `stripInterceptFraction`, `effectiveStemWidth`, `StemWidthModel` |
| (7) | Crossing extension | `raylib/raydensitygrid.cpp`: `DensityGrid::calculateDensities` |
| (9) | Response R and its inversion | `raylib/raystemmodel.cpp`: `freePathResponse`, `StemResponseTable`, `invertStemResponse`; `raylib/raystemmodel.h`: `invertDensityResponse`; `raylib/raystemresponsetable.h` |
| (10) | Cell outputs | `raylib/raydensitygrid.h`: `DensityCell`; `raylib/raydensitygrid.cpp`: `aggregateDensityCells`; `raylib/rayrenderer.cpp`: `writeDensityCells`, `writeDensityProfile`; `raylib/raydensity.cpp`: `estimateDensityTile` |
| — | Neighbour top-up, observed voxels | `raylib/raydensitygrid.cpp`: `addNeighbourPriors`, `estimateDensityGrid` |
| — | Synthetic scenes | `raylib/raypneumatophoregen.cpp`: `PneumatophoreGen::generate`, `castRay` |

## 12. Reproducibility

All scenes are `raycreate pneumatophores <seed> <scene options> --footprint 0.01,0.003,0.05 --truth_cells 1,1,0.1`,
and all estimates `rayrender pneumatophores.ply top density --density_model vertical_cylinder --estimator pooled
--no_flat_top --pixel_width 0.05 --min_bound -3,-3,0 --max_bound 3,3,<0.1 × layers> --cells cells.csv
--footprint 0.01,0.003,0.05 --stem_diameter 0.01 <variant options>`, compared with `pneumatophores_truth_cells.csv`
by summing the 36 cells of each 0.1 m layer, as in `scripts/pneumatophore_seed_study.py`.

| Result | Seeds | Produced by |
|---|---|---|
| Table 1 | 101–108 | `scripts/pneumatophore_seed_study.py --seeds 8` (and `--extra "--pool_radius 0.1"`) |
| Table 2, Section 5 (correction effect) | 101–104 | scenes `--stems_per_m2 30/100/300/600/1000/1500`; variants: none, `--no_overlap_correction`, `--estimator crossing`, `--estimator beer_lambert --no_overlap_correction` |
| Table 3 | 101–104 (0.3 m patches: Table 1 run) | scenes `--stems_per_m2 300 --patchy 0.5,0.3,0`, `--patchy 0.2,0.3,0`, `--stems_per_m2 150 --clustered --cluster_spread 0.03`; variants: none, `--pool_radius 0.1`, `--no_overlap_correction`. Profile layers: `scripts/pneumatophore_validation_plots.py` (seed 11) |
| Section 7.4, angles | 101–104 | scenes 300, 1000; variants: none, `--min_sin_zenith 0.9`, 4 layers |
| Section 7.4, ray count | 101 | scene 300 with `--rays 500000`, default, `--rays 8000000` |
| Table 4 | 101–104 | variants `--footprint 0.01,0.003,0.02 / 0.05 / 0.15` |
| Table 5 | 101 | variants `--stem_diameter 0.008 / 0.012` |
| Table 6 | 101–104 | scenes with `--min_gap 0.5`, `--min_gap 1`; variants `--min_gap 0` and the true gap |
| Return convention, Section 7.6 | 101–104 | scenes and variants `--return_range surface / outline` |
| Table 7 | 101–104 | variants `--pixel_width 0.025 / 0.05 / 0.1` |
| Beam, Section 7.6 | 101–104 | variant `--footprint 0,0,0.05` (with and without `--no_overlap_correction`); scenes and estimate with `--footprint 0.001,0.0001,0.05` |
| Section 7.7, tiling and slope | 101–104 (tiling: 101) | `raydensity` on each scene and on it sheared by z += 0.1 x − 0.05 y, with `--ground` (a plane mesh), absolute, or flattened with `raytranslate subtract`; `scripts/raydensity_tiles.sh` with 2 m tiles |
| Section 7.8 | 31–34 | `scripts/pneumatophore_trunk_study.py --seed <s>`; trunks alone: `raycreate ... --stems_per_m2 0.01 --trunks 0.1,0.2` |
| Section 8.2 table, Section 8.5 (5 cm layers) | 101 | scenes 300 and 1000, `--return_range surface / outline`, layers to 0.5 m; `--cell_size 6,6,0.05` |
| Uniform 0.6 m stems (Sections 8.2, 8.3) | 101 | a small program setting `PneumatophoreParams::mean_height = 0.6`, `height_sd = 0` (raycreate has no height option) |
| Section 5, response accuracy | 6 Monte Carlo seeds | `ray::freePathResponse` with other seeds, and against `StemResponseTable::response` |

The multi-variant runs above used a driver equivalent to `scripts/pneumatophore_seed_study.py` (generate once per
seed, estimate with each variant), which is not part of the repository. Plot totals are deterministic for a seed,
so the numbers are reproducible exactly with the same build.
