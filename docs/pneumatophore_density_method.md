# Pneumatophore frontal area density from hand-held lidar: method, corrections and trunks

A brief technical note on the methodological basis of the density estimate (`--density_model vertical_cylinder`
in `rayrender` and `raydensity`), its corrections, how they were validated, and the effect of sparse trunks.
All numbers are from synthetic hand-held scans (`raycreate pneumatophores`), with a 10 mm, 3 mrad beam and 5%
detection threshold unless stated, and are reproducible with the scripts in `scripts/`. They have not yet been
checked against field data.

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

with H the number of rays ending (hitting) within it, and Lᵢ the path length of ray i up to its end or exit.
Rays with no return, and ground returns (given zero alpha), contribute path but no hit, so they act as censored
observations. Rays whose end lies outside the region analysed are also treated as misses within it.

## 3. Beam footprint

A pulse has a finite, range-dependent spot (1/e² diameter w(r) = w₀ + γ r). Modelling the beam as Gaussian and the
stem as a strip of width d, a return is triggered when the intercepted fraction of the beam's power exceeds a
threshold τ. Thin stems are therefore detected over an effective width d_eff(r) > d (for 1 cm stems: about 2.1 cm
at 1 m and 3.0 cm at 5 m), and the interception rate is a·(d_eff/d)·sin θ. Each ray segment is weighted by
G·d_eff(r)/d, so the estimate remains a = N d for the true mean diameter d. Without this correction, 1 cm stems
were overestimated by about 2.5×.

w₀ and γ come from the scanner's beam specification. τ is not a published parameter (it depends on the detector,
return mode, reflectance and range) and must be calibrated: changing it from 0.05 to 0.02 or 0.15 changed estimates
by about −13% and +30%. Calibration options: counted reference quadrats; the apparent width of isolated stems
against range in the point cloud; or requiring consistent estimates from near and far rays over the same cells.

## 4. Estimators and their biases

The voxel size (5 cm) is comparable to the stem spacing, so voxels hold zero, one or a few discrete stems, which
violates the assumption behind the standard per-voxel Beer–Lambert estimate (many small, randomly placed elements).
Three behaviours were found:

| Estimator | Behaviour |
|---|---|
| Per-voxel Beer–Lambert, averaged | +20 to +60%: a voxel holding one stem intercepts a large fraction of crossing rays, and averaging noisy ratios of few hits is biased high |
| Crossing (hit rays count the full voxel chord, extended through the stem's mean interior chord π d_eff/4) | Exact for discrete stems that don't hide each other within a voxel; −15% at 10 cm voxels and 300 stems/m², −20% at 600 stems/m² with the footprint, from stems hiding each other |
| **Pooled** (recommended): Beer–Lambert path, with hits and path pooled over a ±5 cm horizontal neighbourhood, then corrected (§5) | Insensitive to voxel size; within about ±8% across 30–1500 stems/m² and patchy layouts |

Pooling trades resolution for bias: it smooths over ±5 cm, so patch edges spread by about that much, and stems
clustered along cable roots more tightly than the window are underestimated (−13% near the ground for roots with
3 cm lateral spread), as any free-path estimate is when pooled over uneven cover.

## 5. Occlusion (space occupancy) correction

The pooled estimate measures hits per unit of free path. For randomly placed, overlapping discs (a Boolean model)
the hit rate per unit free path is exactly N d_eff, whatever the overlap. Real stems don't overlap, and in a
non-overlapping (hard-core) arrangement the footprint-widened stems leave less free space than overlapping ones
would, so the measured rate exceeds N d_eff. With a perfect beam (d_eff = d) it is at least a/(1 − φ), because
the stems' own cross-sections are excluded from the free path. It is a little more because rays enter the stem
layer at random points, so sample the gaps between stems in proportion to their length, and hard-core gaps are
more regular than exponential.

The response R(φ, ρ, g) = measured / true N d depends only on the solid fraction φ, the widening ratio
ρ = d_eff/d, and the minimum gap g between stem surfaces (in diameters; centres at least d(1 + g) apart). No
closed form was accurate, so R is computed by a deterministic 2D Monte Carlo simulation (random sequential
addition of 1500 non-overlapping discs, 20 000 lines from random free points, fixed seed), tabulated on demand
and interpolated (e.g. R(0.1, 2.25, 0) = 1.22; R(0.3, 2.25, 0) = 2.15). For each voxel, the neighbourhood's
measured density m and path-weighted mean ρ are used to solve m = a·R(π d a/4, ρ, g) for a.

Limits: the solid fraction of the exclusion discs is capped at 0.4 (below the 0.547 random packing limit), so
N d ≤ 1.6/(π d (1 + g)²), i.e. 51 m⁻¹ (5100 stems/m²) for 1 cm stems with no gap. Above this, voxels are clamped
and reported. The correction needs the mean stem diameter, and becomes sensitive to it at high densities. At
typical Avicennia densities (100–400 stems/m², 0.5–1 cm) φ is 0.01–0.03 and the correction is a few percent,
mostly due to the footprint widening.

The minimum gap matters at high density: stems placed one diameter apart and corrected assuming no gap were
overestimated by 23% at 1000 stems/m² (10% at 300 stems/m²); assuming the right gap removed this. It should be
set from field observations of how closely stems pack.

## 6. Validation (synthetic, pooled estimator, 5 cm voxels, 1 × 1 × 0.1 m cells)

| Case | Estimate / truth, lowest 0.1 m (and 0.1–0.2 m) |
|---|---|
| Random, 30 / 100 / 300 / 1000 / 1500 stems/m² | 1.06 / 1.04 / 1.03 / 0.95 / 0.93 (0.1–0.2 m: 1.10, 1.08, 1.08, 1.02, 0.99) |
| Without the occlusion correction, 300 / 1500 stems/m² | 1.10 / 1.28 |
| Patches of 0.5 m and 0.2 m radius, 300 stems/m² within | 0.99 (1.04) and 0.97 (1.01) |
| Cable roots, 3 cm lateral spread | about 0.87 |
| Sloped ground (10%), heights above a ground mesh (`raydensity --ground`) | within 0.4% of flat ground |
| Same, but flattened with `raytranslate subtract` | 0.5–0.7 (it shifts each whole ray by the ground at its end) |
| Tiled (`raysplit tiles` + `raydensity`) versus a single run | equal to float rounding |

Residual biases: about +10 to +15% in the layer holding the stem tops (top hits are not modelled; using only
shallower rays reduces this); about −7% in the lowest 0.1 m above 1000 stems/m²; +6% in very sparse areas
(30 stems/m², from few hits per neighbourhood; a wider pool radius reduces it). Each scene is a single random seed
on a 6 × 6 m plot, so individual figures carry a few percent of sampling noise.

## 7. Sparse trunks in the bottom 0.3 m

Scenario: 300 stems/m² with trunks at 0.1 per m² (mean diameter 0.2 m), so a cell holding a trunk has 0.21 m⁻¹
of trunk frontal area, about 7% of its pneumatophore N d. Results for 1 × 1 m cells (three trunks in the plot):

| | Cells with a trunk: estimate / pneumatophore truth, 0–0.1 / 0.1–0.2 / 0.2–0.3 m | Cells without |
|---|---|---|
| A: trunks not separated | 1.11 / 1.14 / 1.29 (1.04 / 1.05 / 1.12 of pneumatophores + trunk) | 1.03 / 1.08 / 1.16 |
| B: trunk returns removed (perfect labels) | 1.07 / 1.09 / 1.18 | 1.03 / 1.08 / 1.15 |
| C: trunks detected above the pneumatophores and removed | 1.06 / 1.08 / 1.18 | 1.03 / 1.08 / 1.15 |

**Impact.** Unseparated, a trunk cell's estimate is close to the combined pneumatophore and trunk frontal area,
but this is a coincidence: alone, trunks are counted at about 1.24× their frontal area (1.40× with the trunk's own
diameter in the footprint model, 3.1× without footprint correction), because a 20 cm trunk spans several voxels,
so violates the stem-scale assumptions (thin elements, footprint widening of d, and the occlusion correction for
d). Trunks should therefore be separated, and their frontal area added from their fitted diameters. A trunk also
shadows the volume behind it, reducing coverage (but not biasing) nearby.

**Separability.** In simulation, trunks were cleanly separable in the bottom 0.3 m. All 8 trunks near the plot were
found by clustering returns in a band above the tallest pneumatophores (0.69–1.29 m here) and fitting a horizontal
circle to each cluster: centres within 0.1 cm, diameters 1.1–1.5 cm too large (the footprint widening, d_eff − D,
which the footprint model can remove). Removing returns within each fitted circle plus 1.5 cm, at all heights,
removed 99.9% of trunk returns and 0.05% of pneumatophore returns in the bottom 0.3 m. The pneumatophore estimate in
trunk cells then exceeded trunk-free cells by 0–3%, consistent with the trunk's 3% share of the cell's ground area:
the estimate is per unit of free (observable) volume, so it should be multiplied by (1 − trunk cross-sectional
area / cell area) for a per-cell value.

Recommended procedure: (1) detect trunks above the pneumatophore layer (a circle fit per cluster, or
`rayextract trunks`), and extrapolate them vertically to the ground; (2) set returns within each trunk's fitted
outline (plus a small margin) to zero alpha, so rays stopped by trunks count as unobstructed up to the trunk;
(3) estimate pneumatophores; (4) add each trunk's frontal area D·Δz/V to its cell(s), with D corrected for the
footprint, and scale the pneumatophore density by the cell's free area fraction.

**Caveats for real stands.** The simulation used perfectly vertical, circular trunks of uniform cross-section and
nothing else above the pneumatophores. Separation will be harder where trunks lean, branch low, flare at the base,
or are multi-stemmed; where saplings or seedlings (1–5 cm) fall between the pneumatophore and trunk scales (those
taller than the pneumatophores would be detected and handled as trunks, shorter ones would be counted as
pneumatophores of the wrong diameter); and where low branches, prop-like roots, debris or epiphytes return in the
detection band. Vertical extrapolation assumes a trunk is visible above the pneumatophores at the same position as
at its base; checking the fitted circle at several heights would detect lean.

## 8. Main assumptions

Stems are vertical cylinders of one (mean) diameter, randomly placed or with a known minimum gap; ground returns are
identified and given zero alpha; heights are measured above a ground mesh rather than by flattening the cloud; the
beam is Gaussian with a known divergence and a calibrated detection threshold; and every voxel used is observed by
enough rays (cells report the observed fraction).
