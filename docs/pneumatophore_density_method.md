# Pneumatophore frontal area density from hand-held lidar: method, corrections and trunks

A brief technical note on the methodological basis of the density estimate (`--density_model vertical_cylinder` in
`rayrender` and `raydensity`), its corrections, how they were validated, and the effect of sparse trunks. All
numbers are from synthetic hand-held scans (`raycreate pneumatophores`, 1 cm stems, a 10 mm, 3 mrad beam and a 5%
detection threshold unless stated), estimated with the pooled estimator at 5 cm voxels, and are reproducible with
the scripts in `scripts/`. Numbers quoted as mean ± s.d. are over independent scenes (seeds). None of this has been
checked against field data, and the scenes are generated with the same idealised stem and beam models that the
corrections assume (§3, §5), so the validation tests the estimator's consistency, not the physical models.

## 1. Target quantity

For vertical stems of diameter d at N stems per m² of ground, the hydrodynamically relevant quantity is the frontal
area per unit volume, a = N d (m⁻¹), as a function of height. Related quantities follow from it: the frontal area
index ∫a dz, the solid (volume) fraction φ = N π d²/4 = π d a/4, and, for randomly placed stems, the probability of
an unobstructed horizontal path of length L, exp(−a L).

## 2. Measurement model

Each lidar ray is a probe of the medium along its path. For a ray at zenith angle θ, the horizontal distance
travelled per metre of ray is sin θ, and a vertical stem is struck if its centre lies within d/2 of the ray's
horizontal track, so the interception rate per metre is

&nbsp;&nbsp;&nbsp;&nbsp;λ(θ) = G(θ) a, with G(θ) = sin θ.

This replaces raycloudtools' default G = 0.5, which applies to isotropically oriented leaves; for stems it would
bias estimates by a factor 2 sin θ (about 2× for the near-horizontal rays of a hand-held scanner). Since G depends
on each ray's direction, each ray's path length in a voxel is weighted by its own G, and the maximum-likelihood
(Beer–Lambert) estimate for a voxel or neighbourhood is

&nbsp;&nbsp;&nbsp;&nbsp;â = H / Σᵢ Gᵢ Lᵢ,

with H the number of rays ending (hitting) within it, and Lᵢ the path length of ray i up to its end or exit. Rays
with no return, and ground returns (given zero alpha), contribute path but no hit, so they act as censored
observations. Rays whose end lies outside the region analysed are treated as misses within it. Rays steeper than
`--min_sin_zenith` (default 0.1) are ignored.

## 3. Beam footprint and return range

A pulse has a finite, range-dependent spot (1/e² diameter w(r) = w₀ + γ r). Modelling the beam as Gaussian and the
stem as a strip of width d, a return is triggered when the intercepted fraction of the beam's power exceeds a
threshold τ. Thin stems are therefore detected over an effective width d_eff(r) > d (for 1 cm stems: 2.07 cm at 1 m
and 3.05 cm at 5 m; ρ = d_eff/d averaged 2.1–2.4 over the voxels of the scans), and the interception rate is
a·(d_eff/d)·sin θ. Each ray segment is weighted by G·d_eff(r)/d, using the range from the ray's original start, so
the estimate remains a = N d for the true mean diameter d. Without this correction, 1 cm stems were overestimated by
2.1–2.6×.

w₀ and γ come from the scanner's beam specification. τ is not a published parameter (it depends on the detector,
return mode, reflectance and range) and must be calibrated: scenes generated with τ = 0.05 and estimated assuming
0.02 or 0.15 read 11–12% lower or 22–27% higher. Calibration options: counted reference quadrats; the apparent width
of isolated stems against range in the point cloud; or requiring consistent estimates from near and far rays over
the same cells.

**Return range.** For a partly intercepted stem (the beam axis misses the stem but enough power hits it), the range
reported by a scanner is an instrument property that has not been verified. `--return_range` selects the model, both
in the generator and in the correction: `surface` (default), the range of the physical stem surface (or of closest
approach), or `outline`, where the beam axis enters the footprint-widened outline. Simulated consistently, the
occlusion response (§5) is nearly the same for both (within 0.2% at φ ≤ 0.1, 0.7% at φ = 0.2), and estimating a
scene with the other convention changed the estimates by at most 0.2%. The generated data do differ, though: scenes
generated with `outline` read 1–2% higher in the lowest 0.1 m and 2.5–4.5% higher at 0.1–0.2 m than `surface`
scenes, with or without the correction, so the instrument's behaviour matters at that level and the option does not
remove it (the cause of the difference was not investigated). (An earlier concern, that the correction over-corrects
by 5–19% if ranges are reported at the stem surface, came from a response simulation whose lines started only
outside the widened outlines; with start points consistent with surface returns, the two conventions agree.)

## 4. Estimators, pooling and cell outputs

The voxel size (5 cm) is comparable to the stem spacing, so voxels hold zero, one or a few discrete stems, which
violates the assumption behind the standard per-voxel Beer–Lambert estimate (many small, randomly placed elements).
Plot totals in the lowest 0.1 m over 30–1500 random stems/m² (4 scenes each, cell means):

| Estimator | Estimate / truth |
|---|---|
| Per-voxel Beer–Lambert (`beer_lambert`), uncorrected | 1.24–1.31 (1.31–1.39 at 0.1–0.2 m): a voxel holding one stem intercepts many crossing rays, and averaging noisy ratios of few hits is biased high |
| Crossing (`crossing`: hit rays count path to the voxel exit, extended through the stem's remaining widened chord, π d_eff/4 from the outline or π(d_eff² + d²)/(8 d_eff) from the surface) | 1.05 at 30, 0.97 at 300, 0.82 at 1000 and 0.74 at 1500 stems/m². It has no occlusion correction, and with a footprint the entries into a widened stem that lie inside another widened stem are never observed, so it reads lower with the footprint than with a near-perfect beam (1 mm, 0.1 mrad), which gave 1.00 at 300 and 0.93 at 1000 stems/m². Not recommended with a footprint |
| **Pooled** (`pooled`, recommended): Beer–Lambert hits and path pooled over a horizontal neighbourhood, then corrected (§5) | 1.01–1.08 (cell means; 0.98–1.01 path-weighted, below) |

**Pool window.** The neighbourhood is the (2k+1) × (2k+1) voxels centred on each voxel in its horizontal layer, with
k = round(`--pool_radius` / voxel width), so it extends ±(k + ½) voxels: ±0.075 m for the default 0.05 m at 5 cm
voxels, ±0.125 m for 0.1 m. Its size changes in steps with the voxel size (the run prints it). Pooling trades
resolution for bias: it smooths over the window, so patch edges spread by about that much, and stems clustered more
tightly than the window (e.g. along cable roots) read low, as any free-path estimate does when pooled over uneven
cover.

**Cells and profiles.** `--cells` (and `raydensity`) aggregate voxels into cells (e.g. 1 × 1 × 0.1 m), and
`--profile` into whole layers, with three densities:

* `density` (main): the mean of the observed voxel densities, each voxel weighted equally (by volume). It follows
  patchy stems, but averaging noisy voxel ratios biases it high where neighbourhoods have few hits: +4–8% for random
  stems at 30–1000 stems/m² in the lowest 0.1 m with the default pool (+1% at 1500), +2–4% with `--pool_radius 0.1`.
* `pooled_density`: total hits over total weighted path of the cell's observed voxels (as for one voxel the size of
  the cell). It is within 2–3% for random stems in the lowest 0.2 m, but rays penetrate dense patches less, so
  patches get little weight: it read 0.47–0.55 of the truth in the lowest 0.1 m for 0.2–0.5 m patches in 1 m cells,
  and 0.3 over a 6 m profile layer (one scene).
* `patchiness_ratio` = `density` / `pooled_density`: 1.03–1.17 for random stems (1.09 at 100–300 stems/m² with the
  default pool, 1.05 with `--pool_radius 0.1`, so mostly the noise bias above), 1.6–2.4 in the lowest 0.2 m for
  patchy stems and cable roots. A ratio near 2 flags cells where `pooled_density` is unusable and `density` should
  be checked against the pool radius.

A voxel with fewer than 10 rays is topped up from its neighbours (a 3 × 3 × 3 prior), but for stems a voxel no ray
reached is not, so voxels below the ground or beyond the scanned volume stay unobserved. A voxel with 2 or fewer
rays is unobserved; a cell's densities are the means or sums over its observed voxels only (NaN if none), and
`observed_fraction` reports their share. The per-voxel (n−1)/n small-sample factor of the original estimator is not
applied to pooled values (it used the voxel's own ray count, not the neighbourhood's).

Voxels are assigned to cells by their centres. With `rayrender` the voxel and cell grids share the min bound, so
layers that are multiples of the voxel align exactly. With `raydensity --ground`, voxels are assigned by the height
of their centre above the ground beneath the column's centre, so layers are misregistered by up to half a voxel (2.5
cm) wherever the ground height isn't a multiple of the voxel; use voxels no larger than half the layer height.

## 5. Occlusion (space occupancy) correction

The pooled estimate measures hits per unit of free path. For randomly placed, overlapping discs (a Boolean model)
the hit rate per unit free path is exactly N d_eff, whatever the overlap. Real stems don't overlap, and in a
non-overlapping (hard-core) arrangement the footprint-widened stems leave less free space than overlapping ones
would, so the measured rate exceeds N d_eff. With a perfect beam (d_eff = d) it is at least a/(1 − φ), because
the stems' own cross-sections are excluded from the free path. It is a little more because rays enter the stem
layer at random points, so sample the gaps between stems in proportion to their length, and hard-core gaps are
more regular than exponential.

The response R(φ, ρ, g) = measured / true N d depends on the solid fraction φ, the widening ratio ρ = d_eff/d, the
minimum gap g between stem surfaces (in diameters; centres at least d(1 + g) apart), and nominally the return range
convention (§3). It is computed by a deterministic 2D Monte Carlo simulation: random sequential addition of 1500
non-overlapping discs and 10 000 lines from random start points, pooled over 64 independent realisations. For g = 0,
a table over 13 solid fractions (0–0.4) and 17 width ratios (1–5) for both conventions is embedded in the library
(`raylib/raystemresponsetable.h`, generated by `scripts/generate_stem_response_table.cpp`) and bilinearly
interpolated; other gaps are simulated when first needed. For example R(0.024, 2.3, 0) = 1.046 (300 stems/m²),
R(0.079, 2.3, 0) = 1.157 (1000 stems/m²) and R(0.2, 2.3, 0) = 1.51. Rerunning the simulation with 6 other seeds
varied R by 0.2–0.6% (s.d.) at φ ≤ 0.2 and up to 0.9% at φ = 0.25, so the table carries about that much noise;
interpolation differed from direct simulation by at most 0.2% for φ ≤ 0.2. For each voxel, the neighbourhood's
measured density m and path-weighted mean ρ are used to solve m = a·R(π d a/4, ρ, g) for a, bracketing upwards from
m (so R < 1 would also be handled).

Limits: the solid fraction of the exclusion discs is capped at 0.4 (below the 0.547 random packing limit), so N d ≤
1.6/(π d (1 + g)²), i.e. 51 m⁻¹ (5100 stems/m²) for 1 cm stems with no gap. Above this, voxels are clamped and
reported. The correction needs the mean stem diameter, and becomes sensitive to it at high densities. In the scans
it lowered estimates by 1–6% at 30–300 stems/m², 11% at 600 and 18–30% at 1000–1500 stems/m². About a third of the
response comes from the footprint widening: R(φ, 1, 0) is 1.03 at 300 and 1.10 at 1000 stems/m² (from the table).

The minimum gap matters at high density: with stems placed at least one diameter apart, assuming no gap
overestimated the cell mean by 34% at 1000 stems/m², 10% at 300 and 4% at 100 (half a diameter apart: 15%, 5% and
2%), relative to assuming the right gap. It should be set from field observations of how closely stems pack.

## 6. Validation (synthetic, pooled estimator, 5 cm voxels, 1 × 1 × 0.1 m cells)

Plot totals (estimate / truth over the 36 cells of each layer), as the mean ± s.d. over 8 independent scenes
(`scripts/pneumatophore_seed_study.py`; stem heights 0.25 ± 0.08 m):

| Scene | `density`, 0–0.1 m | 0.1–0.2 m | 0.2–0.3 m | 0–0.1 m, `--pool_radius 0.1` | `pooled_density`, 0–0.1 m | No correction, 0–0.1 m |
|---|---|---|---|---|---|---|
| Random 100 stems/m² | 1.073 ± 0.006 | 1.085 ± 0.006 | 1.123 ± 0.006 | 1.042 ± 0.007 | 0.983 ± 0.007 | 1.099 |
| Random 300 stems/m² | 1.065 ± 0.005 | 1.077 ± 0.003 | 1.113 ± 0.003 | 1.033 ± 0.006 | 0.979 ± 0.006 | 1.126 |
| Random 1000 stems/m² | 1.040 ± 0.005 | 1.075 ± 0.003 | 1.109 ± 0.002 | 1.021 ± 0.006 | 1.010 ± 0.005 | 1.231 |
| Patches of 0.3 m radius (40% cover), 300 stems/m² within | 1.015 ± 0.005 | 1.040 ± 0.005 | 1.093 ± 0.005 | 0.925 ± 0.008 | 0.547 ± 0.026 | 1.065 |

The standard deviations between scenes are under 1% (except `pooled_density` in patches), so the differences between
rows are systematic biases, not sampling noise. Individual 1 m cells have a relative RMSE of 5–9% in the lowest two
layers (`density`). The pool radius trades two biases: a larger window has more hits per neighbourhood, reducing the
upward bias of noisy ratios, but averages over patches smaller than the window, which reads low. `--pool_radius 0.1`
suits fairly uniform stands, and the default patchy ones.

Other cases (4 scenes each, `density` in the lowest 0.1 m and 0.1–0.2 m; s.d. under 1%):

| Case | Default pool | `--pool_radius 0.1` |
|---|---|---|
| Patches of 0.5 m radius, 30% cover, 300 stems/m² within | 1.03, 1.05 | 0.95, 0.99 |
| Patches of 0.2 m radius, 30% cover | 0.99, 1.01 | 0.88, 0.92 |
| Cable roots, 150 stems/m², 3 cm lateral spread | 0.87, 0.90 | 0.77, 0.80 |
| Sloped ground (10% in x, 5% in y; 300 stems/m²), `raydensity --ground` | within 0.5% of flat ground over the plot, cells within 3.5% | |
| Same, heights left absolute | 0.72–0.95 of flat ground | |
| Same, flattened with `raytranslate subtract` | 0.43–0.71 of flat ground (it shifts each whole ray by the ground at its end) | |
| Tiled (`scripts/raydensity_tiles.sh`, 3 × 3 tiles) versus a single run (one scene) | identical | |

Residual biases: (i) the noise bias of `density` above. (ii) Stem tops: steep rays end on stem tops, which the model
doesn't include, and this biases every layer in which stems end, not only the top one. With 3%, 24%, 47% and 24% of
stems ending in successive 0.1 m layers, `pooled_density` at 300 stems/m² read 0.98, 0.99, 1.00 and 1.08. Using only
shallow rays (`--min_sin_zenith 0.9`) lowered every layer at 300 stems/m² (by 1–10%), not only the upper ones, and
at 1000 stems/m² so few shallow rays reach the ground that the lowest layer read 0.35, so it is not a remedy. (iii)
At 1000 stems/m², `pooled_density` read 0.95 at 0.2–0.3 m (cause not established).

## 7. Sparse trunks in the bottom 0.3 m

Scenario: 300 stems/m² with trunks at 0.1 per m² (mean diameter 0.2 m), so a cell holding a trunk has 0.18–0.25 m⁻¹
of trunk frontal area, about 6–8% of its pneumatophore N d. Results for 1 × 1 m cells, 4 scenes (2–5 trunks in each
6 × 6 m plot), `density` / pneumatophore truth (s.d. 0.5–3% in trunk cells, under 1% elsewhere):

| | Cells with a trunk, 0–0.1 / 0.1–0.2 / 0.2–0.3 m | Cells without |
|---|---|---|
| A: trunks not separated | 1.15 / 1.17 / 1.23 (1.06 / 1.07 / 1.06 of pneumatophores + trunk) | 1.07 / 1.08 / 1.12 |
| B: trunk returns removed (perfect labels) | 1.10 / 1.12 / 1.14 | 1.07 / 1.08 / 1.12 |
| C: trunks detected above the pneumatophores and removed | 1.10 / 1.11 / 1.14 | 1.07 / 1.08 / 1.12 |

**Impact.** Unseparated, a trunk cell's estimate is close to the combined pneumatophore and trunk frontal area, but
this is a coincidence: alone (scenes with almost no pneumatophores), trunks read 0.70 ± 0.12 of their frontal area
with the stem settings, 0.82 ± 0.15 with the trunk's own diameter, and 1.76 ± 0.31 without footprint or occlusion
correction, varying widely between scenes. A 20 cm trunk spans several voxels, so violates the stem-scale
assumptions (thin elements, footprint widening of d, and the occlusion correction for d). Trunks should therefore be
separated, and their frontal area added from their fitted diameters. A trunk also shadows the volume behind it,
reducing coverage nearby.

**Separability.** In simulation, trunks were cleanly separable in the bottom 0.3 m. In each scene, every trunk
within 4 m of the plot centre (5–8 per scene) was found by clustering returns in a band above the tallest
pneumatophores (about 0.6–1.3 m) and fitting a horizontal circle to each cluster: centres within 0.2 cm and
diameters 0–0.4 cm too large (returns are reported at the physical surface, §3; with `outline` they would be widened
by the footprint). Removing returns within each fitted circle plus 1.5 cm, at all heights, removed 99.8–100% of
trunk returns and 0.03–0.06% of pneumatophore returns in the bottom 0.3 m. The pneumatophore estimate in trunk cells
then exceeded trunk-free cells by 2–3%, consistent with the trunk's 1–8% share of the cell's ground area: the
estimate is per unit of free (observable) volume, so it should be multiplied by (1 − trunk cross-sectional area /
cell area) for a per-cell value.

Recommended procedure: (1) detect trunks above the pneumatophore layer (a circle fit per cluster, or `rayextract
trunks`), and extrapolate them vertically to the ground; (2) set returns within each trunk's fitted outline (plus a
small margin) to zero alpha, so rays stopped by trunks count as unobstructed up to the trunk; (3) estimate
pneumatophores; (4) add each trunk's frontal area D·Δz/V to its cell(s), with D corrected for the footprint if
needed, and scale the pneumatophore density by the cell's free area fraction.

**Caveats for real stands.** The simulation used perfectly vertical, circular trunks of uniform cross-section and
nothing else above the pneumatophores. Separation will be harder where trunks lean, branch low, flare at the base,
or are multi-stemmed; where saplings or seedlings (1–5 cm) fall between the pneumatophore and trunk scales (those
taller than the pneumatophores would be detected and handled as trunks, shorter ones would be counted as
pneumatophores of the wrong diameter); and where low branches, prop-like roots, debris or epiphytes return in the
detection band. Vertical extrapolation assumes a trunk is visible above the pneumatophores at the same position as
at its base; checking the fitted circle at several heights would detect lean.

## 8. Main assumptions and limits

Stems are vertical cylinders of one (mean) diameter, randomly placed or with a known minimum gap; ground returns are
identified and given zero alpha; heights are measured above a ground mesh rather than by flattening the cloud; the
beam is Gaussian with a known divergence and a calibrated detection threshold; partial returns are reported as the
chosen `--return_range` convention (unverified for real scanners); stem tops are not modelled; and every voxel used
is observed by enough rays (cells report the observed fraction). The validation is synthetic and uses the same beam,
stem and return models as the corrections, so it cannot detect errors in those models: field validation against
counted quadrats is still needed, in particular of τ, the return range behaviour, and the high-bias of `density`.
