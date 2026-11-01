# I2EC_3d — imposed earthquake cycles on a 3D megathrust

Goal: a MATLAB analogue of [seqeas](https://github.com/tobiscode/seqeas-public) built on
`unicycle`. Earthquakes are imposed on a triangular fault mesh, either periodically on a
locked part or as a sequence of events with any slip distribution, and the rest of the
mesh evolves as an initial value problem governed by steady-state rate-dependent friction.

## Contents

**Reference** · [Conventions](#conventions) · [Physics](#physics) · [Code layout](#code-layout) · [Tests](#tests)

**Tasks**

| # | task | status | checklist |
|---|---|---|---|
| 1 | [Mesh](#task-1) | ✅ Done | 1/1 |
| 2 | [Kernels](#task-2) | ✅ Done | 2/2 |
| 3 | [Interseismic ODE (rate-dependent friction)](#task-3) | ✅ Done | 4/4 |
| 4 | [Interseismic verification](#task-4) | ✅ Done | 5/5 |
| 5 | [Imposed earthquakes](#task-5) | ✅ Done | 8/8 |
| 6 | [Cycle driver and spin-up](#task-6) | ✅ Done | 5/6 |
| 7 | [Outputs](#task-7) | 🔄 In progress | 1/2 |
| 7b | [Instrumental era](#task-7b) | ⬜ Not started | 0/1 |
| 8 | [Later / optional](#task-8) | 💤 Backlog | 0/5 |
| 9 | [Linear viscous friction: analytic cycles and inversion](#task-9) | 🔄 In progress | 3/7 |
| 10 | [Time-windowed effective viscosity: inversion design](#task-10) | 📝 Planning | 0/5 |
| 11 | [Inverse package: design](#task-11) | 📝 Planning | 1/5 |
| 12 | [Experiment: effective viscosity of an RSF cycle, by time window](#task-12) | 📝 Planning | 0/6 |
| 13 | [Experiment: stochastic event sequence without a locked zone](#task-13) | ✅ Done | 3/3 |

Status: ✅ done · 🔄 in progress · 📝 planning · ⬜ not started · 💤 backlog

---

## Conventions

- **Units:** m, s, MPa. Velocities are in m/s and time is in seconds; convert to years only for plotting.
  - With `G = 30e3` (MPa) the traction kernels are in MPa/m, independent of the time unit.
  - The `triangleReceiver` constructor fills in defaults, some in m/yr (`Vs = 3e3*3.1536e7`). These defaults are not used.
- **Coordinates:** the axis order differs between the files and the loaded object.
  - The `.ned` columns are `index, x1, x2, x3` = north, east, depth (positive down).
  - On reading, `triangleReceiver` sets `rcv.x = [x2, x1, -x3]`, so `rcv.x` and `rcv.xc` are **[east, north, up]**. The first two columns are swapped and depth is negated.
  - When writing a mesh (`make_fault_mesh.m`), write `[north, east, depth]`. When reading one, use `rcv.xc(:,1)` = east, `rcv.xc(:,2)` = north and `rcv.xc(:,3)` = −depth.
  - Checked on `generic_sz15`: `rcv.xc(:,1)` spans ±247 km along strike, `rcv.xc(:,2)` spans 3–287 km (north, down-dip) and `rcv.xc(:,3)` spans −1 to −77 km.
- **Sign conventions** (checked numerically):
  - Positive dip slip is thrust.
  - Positive strike slip is left-lateral.
  - The diagonals of `Kss` and `Kdd` are negative.
- **Kernel:** `bigK = [Kss,Kds;Ksd,Kdd]`, so `[τs;τd] = bigK·[s_s;s_d]`.
  - `Kds` is the strike traction due to dip slip.
- **Mesh object:** a `unicycle.geometry.triangleReceiver` (a handle class), referred to as `rcv` or `mesh`. It already carries everything the model needs:

  | property | meaning |
  |---|---|
  | `a`, `b` | rate-and-state parameters |
  | `sigma` | effective normal stress |
  | `mu0` | reference friction coefficient (f0) |
  | `Vo` | reference slip rate (v0) |
  | `Vs` | shear-wave speed, for radiation damping `η = G/(2Vs)` (`G` from `earthModel.G`) |
  | `Vpl`, `Vrake` | long-term slip rate magnitude and rake (degrees) |
  | `pinnedPosition` | indices of patches held at v = 0 |
  | `earthModel` | elastic properties `G` and `nu` |

  **Workflow:** the user reads the mesh, then writes all of these properties into the mesh object in m/s and MPa. The model code only reads them from `rcv`; it never sets physical parameters or relies on the constructor defaults.

  Don't build a parallel parameter struct. Only the derived quantities go in a separate struct: partitioned kernels, the loading rate and the initial state.
  - Before solving, check each property's size (N × 1, or a scalar expanded to N × 1).
  - Check `a > b` on creeping patches.
  - Warn if `Vs` still holds the m/yr constructor default.

<sub>[↑ Contents](#contents)</sub>

---

## Physics

- **Friction** (steady state, with v collinear with τ; collinearity is an assumption). Two options, chosen by `law` in `assembleODE`:
  - `'log'`: `|τ| = σ (f0 + (a−b) ln(|v|/v0))`
  - `'asinh'` (default): `|τ| = (a−b)σ asinh(|v|/(2v0) · exp(f0/(a−b)))`
  - They are identical unless `|v| ≲ v0·exp(−f0/(a−b))`, where the log law gives a negative |τ|; the asinh form stays positive as v → 0.
  - **(f0, v0) enter only through one speed, `v* = v0·exp(−f0/(a−b))`.**
    - log: `|τ| = (a−b)σ ln(v/v*)`; asinh: `|τ| = (a−b)σ asinh(v/(2v*))`.
    - Any (f0, v0) pair with the same v* is the same model. Changing v0 just shifts f0 by `(a−b) ln(v0'/v0)`.
    - The physically meaningful number is the friction at the operating speed, `f(Vpl) = f0 + (a−b) ln(Vpl/v0)`, which gives `v* = Vpl·exp(−f(Vpl)/(a−b))`.
    - Prefer `v0 = Vpl`, so that f0 is the friction at plate rate.
  - **The log law fails at small friction, and no choice of v0 fixes it.** It needs |τ| < 0 when v < v*, i.e. slip against its own traction, which is unphysical.
    - For f(Vpl) = 0.01 and a−b = 0.005, v* ≈ 0.14·Vpl. Any patch slowed by about 7× (e.g. in an asperity's stress shadow) breaks the law.
    - With v0 = Vpl the log law is unsafe for f0 ≲ 0.02–0.05 next to a locked asperity, and fine for f0 ≳ 0.1.
  - **Where the asinh law comes from:** thermally activated slip with forward and backward jumps, `v = v*(e^{τ/A} − e^{−τ/A}) = 2v* sinh(τ/A)` with A = (a−b)σ. The log law keeps only the forward jumps.
    - For v ≫ v*, asinh is classical rate-and-state friction (identical to log).
    - For v ≲ v* it becomes **linear viscous**: `v ≈ (2v*/A)·τ`, with viscosity A/(2v*). Traction goes to 0 with v and is never negative.
    - It does not by itself reverse slip. sinh is odd, so v ∥ τ always; slip reverses only if the traction vector reverses, which asinh handles smoothly and log cannot.
  - **⚠ Low-f runs are not classical RSF.** When `f(Vpl)/(a−b)` is O(1–5), the asinh fault sits in or near its linear-viscous regime over the operating speeds. Its dynamics are a different style: `dF/dζ ≪ (a−b)σ`, velocity responds nearly linearly to stress, and the speed time scale is shortened.
    - Interpret such runs as a regularized, viscous-like fault, not as rate-strengthening friction.
    - The first f0 sweep (`studies/f0ConvergenceStudy.m`) used v0 = 1e-6 m/s with Vpl ≈ 1.7e-9 m/s, so f(Vpl) = f0 − 0.032. Its f0 = 0.01 and 0.02 runs have negative friction at Vpl under log (v*/Vpl ≈ 80 and 11), i.e. they are fully in the viscous regime. The sweep was rerun with v0 = Vpl and f0 ≥ 0.1 (see [Task 4](#task-4)).
- **Decisions (2026-10-30).**
  - **asinh only from here on.** The log law stays in the code (and its tests) for reference, but all experiments use `'asinh'`: the log law cannot represent a reduced or reversed traction (|τ| < 0 below v*), which imposed events next to or inside creeping regions produce.
  - **f0 ≥ 0.1 only (with v0 = Vpl).** We do not think f0 < 0.1 is stable or mechanically meaningful enough to explore through experiments:
    - f(Vpl)/(a−b) = O(1–10) puts the fault in or near the linear-viscous regime of asinh, so it is no longer rate-and-state friction (see the ⚠ note above).
    - |τ| ≈ σ f0 is then only a few MPa, comparable to coseismic stress changes (Δτ ~ 3–12 MPa near the rupture). Events then reverse tractions and rotate slip by up to 180° ([Task 5](#task-5)), and the (ζ, θ) state approaches its singularity.
    - The system becomes stiff for `ode45` (the f0 = 0.01 sweep with v0 = 1e-6 did not finish one run in ~2 h).
  - f0 = 0.1 is the lower end of the explored range and is checked in every new experiment.
- **State per creeping patch:** `[slip_s, slip_d, ζ = ln(|v|/v0), θ = rake]`.
  - This state is singular at v = 0: ζ → −∞, and `θ̇ = τ̇⊥/(|τ| + η|v|)` becomes very large as |τ| → 0. That is stiff but finite unless the traction passes exactly through zero.
  - Option, if strong stress shadows or slip reversal matter (especially in the viscous regime): use the traction components (τs, τd) as the state, with v from the asinh law. That has no singularity.
- **Elasticity:** `τ̇ = K (v − Vpl)` with **v = 0 on pinned patches**. On creeping patches this is `τ̇_c = K_cc (v_c − Vpl_c) + K_cp (0 − Vpl_p) = K_cc v_c − K(c,:)·Vpl`, so the pinned columns enter only through the constant `τ̇_load`. Checked against the full `K (v − Vpl)` with `v_p = 0`: relative error 4e-15.
- **Projecting τ̇ onto the slip direction and its normal:**
  - Parallel: `τ̇·v̂ = (d|τ|/dζ + η|v|) ζ̇`, with `d|τ|/dζ = (a−b)σ` (exact for log, and to round-off for asinh). This sets acceleration and deceleration; f0 and v0 drop out, exactly as in 2D.
  - Perpendicular: `τ̇·v̂⊥ = (|τ| + η|v|) θ̇`. This sets rotation; because `|τ|/(a−b)σ = f0/(a−b) + ζ`, the static friction level matters.

<sub>[↑ Contents](#contents)</sub>

---

## Code layout

| function | role |
|---|---|
| `unicycle.greens.computeKernels(rcv,name,kdir,nCheck)` | load or compute cached traction kernels |
| `unicycle.ode.assembleODE(rcv,bigK,law)` | one-time setup: partitioned kernel (`Kcc`, `Kcp`), loading, `y0` |
| `unicycle.ode.steadyState(sys)` | interseismic steady state `v = K_cc \ (−τ̇_load)` (asperity locked forever) as a state vector; reference for single-event experiments, not the cycle-spin-up start (see [Task 6](#task-6)) |
| `unicycle.ode.eventStressChange(sys,Trec,weight,weightC)` | coseismic slip `weight·Vpl·Trec` on pinned patches (and `weightC·Vpl·Trec` on creeping patches inside the rupture); `Δτ = K_cp slip_P + K_cc slip_C` |
| `unicycle.ode.ruptureSlipForStressChange(sys,Trec,inRuptureC,dtauTarget,direction)` | slip on creeping patches inside the rupture that gives each of them a prescribed Δτ·d (e.g. a 3 MPa stress drop) |
| `unicycle.ode.periodicEvent(sys,Trec,weightC)` | event that releases exactly Vpl·Trec on the locked patches (closes their slip budget), plus optional slip on creeping rupture patches |
| `unicycle.ode.sequenceEvent(sys,t,slip)` | one event of a sequence at time t from its N × 2 slip on any patches: `slipP`, `slipC` and `Δτ = K_cp slip_P + K_cc slip_C`; collect them in `seq = struct('T',T,'events',events)` for `spinUp`, `runCycles` or `linvisSpinUp` |
| `unicycle.ode.spinUp(sys,event,opts)` / `spinUp(sys,seq,opts)` | spin-up only: repeats the characteristic event, or a sequence (`T`, `events`, as for `linvisSpinUp`), until the (super)cycle is periodic and returns the primed state a time `tSince` after its start (default: just before the next one), plus convergence info and the state before/after each event. `opts.repeat = false` runs a non-repeating sequence once (e.g. a stochastic catalogue, [Task 13](#task-13)), with optional dense output per interval |
| `unicycle.ode.runCycles(sys,event,opts)` / `runCycles(sys,seq,opts)` | `spinUp` + one stored full (super)cycle with dense output after every event |
| `unicycle.greens.computeSurfaceKernels(rcv,xObs)` | surface displacement kernels at stations, (3D × 2N) |
| `unicycle.ode.surfaceDeformation(G,sys,t,Y,opts)` | surface velocity and displacement time series, optionally relative to the long-term motion |
| `unicycle.ode.stateTraction(sys,y)` | traction vector implied by the state, `|T| = F(v) + ηv` along θ |
| `unicycle.ode.applyEvent(sys,y,dtau,slipC)` | coseismic state jump `T⁺ = T⁻ + Δτ`, damped inverse for v⁺, θ⁺ from T⁺; adds `slipC` to the slip state |
| `unicycle.ode.odeRateDependent(t,y,sys)` | ODE right-hand side |
| `unicycle.ode.frictionLaw(sys,zeta)` | `F = |τ|` and `dF/dζ` for `'log'` or `'asinh'` |
| `unicycle.ode.frictionInverse(sys,tauMag,useDamping)` | `ζ` from the traction, with radiation damping by default |

`+unicycle/Contents.m` lists an older `unicycle/ode/*` (evolution, rateandstate, …) that is not in this copy, so there are no name clashes.

<sub>[↑ Contents](#contents)</sub>

---

## Tests

`tests/`, run with `matlab -batch "run('tests/runAllTests.m')"`. There are 85 tests and all pass; the run takes about 1 min.

| file | covers |
|---|---|
| `testGreensConventions` | `.ned` → `rcv.x` axis swap; unit vectors; positive dip slip = thrust; positive strike slip = left-lateral; negative self-stiffness |
| `testComputeKernels` | block ordering; reload; changed G, nu or geometry → new file; tampered or foreign cache caught |
| `testFrictionLaw` | log law; asinh = log for realistic f0/(a−b); asinh regularization at low v; derivative; inverse round trip; no overflow; damped inverse |
| `testAssembleODE` | partitioning; scalar expansion; Vpl/Vrake components; loading `K(v−Vpl)` with `v_p = 0`; steady state; `y0`; validation and warning |
| `testImposedEvent` | `Kcp` partition; `Δτ` matches the full kernel and is linear in the weight; nearest neighbours loaded in the slip direction; `applyEvent` ∘ `stateTraction` is the identity for Δτ = 0; traction balance and direction after rotating/amplifying/reducing steps; damping caps v⁺; spring-slider afterslip → Δτ/k; `steadyState` is a fixed point; log-law negative traction raises an error; rupture onto creeping patches (Δτ from the full kernel, stress drop, slip jump); traction reversal with asinh gives backward slip ∥ T⁺; slip solved for a prescribed stress change hits the target exactly |
| `testRunCycles` | locked-budget closure gives Δs* = Vpl·Trec on creeping patches (with and without creeping rupture slip); spring-slider cycles converge to Vpl·T + Δτ/k; small-mesh spin-up at f0 = 0.1 closes the slip budget and repeats; warning when the locked budget does not close; `spinUp` matches the `runCycles` cycle at tSince = Trec, mid-cycle and 0; primed state independent of the starting state (contraction < 1); tSince outside [0, Trec] rejected; sequences: one event = sequence of one event at t = 0 (dtau from slip when omitted); one pass through a repeated event = the spin-up loop; periodic supercycle of three events on/off the locked zone closes the budget (Vpl·T on creeping patches), stored supercycle ordered with each event twice, warning when an event is dropped; no locked zone, one pass: jumps = `applyEvent`, slip accounting, dense output joins the event states, stopping at tSince truncates consistently |
| `testSurfaceDeformation` | surface kernel layout and slip signs; long-term removal; back slip is landward; displacement = ∫ velocity; coseismic step |
| `testLinvis` | linear viscous model: partition and v* vs RSF set-up; eigendecomposition; analytic v, s vs ode45 of `odeLinearViscous`; closed-form periodic cycle vs numerically spun-up cycle; fixed point; η-time scaling; event jump; exact periodic state vs iterated cycle map; slip budget independent of η; surface output; event sequences: closed form vs numerical spin-up, superposition of shifted single-event cycles, locked-budget check |
| `testKLBasis` | KL basis: orthonormality (area weighted), covariance reconstruction, zero-mean modes, default ℓ, long-ℓ limits, Matérn forms, mesh independence, patch subsets |
| `testOdeRateDependent` | projection formulas; f0-independent acceleration; `θ̇ ∝ 1/|τ|`; spring-slider logistic solution; rake relaxation `exp(−k Vpl t/|τ|)`; radiation damping; damped traction balance conserved; small mesh converges to BEM |

Helpers: `makeTestMesh` (60 × 30 km plane, 36 triangles), `setTestProperties` and `oneSliderSystem`.

<sub>[↑ Contents](#contents)</sub>

---

## Tasks

<a id="task-1"></a>

### 1. Mesh

> **Status:** ✅ Done · checklist 1/1

- [x] `make_fault_mesh.m`: planar 500 × 300 km megathrust dipping 15°, 3000 triangles of about 10 km.
  - Writes `fault/generic_sz15.{ned,tri}`.

<sub>[↑ Contents](#contents)</sub>

---

<a id="task-2"></a>

### 2. Kernels

> **Status:** ✅ Done · checklist 2/2

- [x] `unicycle.greens.computeKernels(rcv,name,kdir,nCheck)` returns `bigK` (2N × 2N shear) and `bigKn = [Ksn,Kdn]`. It caches them in `kernels/<name>_<hash>.mat`.
  - The hash is a SHA-256 of `rcv.x`, `rcv.vertices`, the Green's function class, G and nu. A different geometry or elastic model gets a different file, so stale kernels can't be picked up.
  - On load it compares the full stored fingerprint, then recomputes `nCheck` (default 2) random source columns and compares them with the cache. That takes about 1 s for 3000 patches, against about 260 s for the full computation.
  - Tested: reload, changed G, moved vertex, a tampered kernel entry, and a cache file renamed from another mesh are all handled or caught.
  - The file is 415 MB because it also stores `Ksn`/`Kdn` (6 blocks instead of 4).
- [x] Cached for `generic_sz15`: `kernels/generic_sz15_0f17a6d68309.mat`.

<sub>[↑ Contents](#contents)</sub>

---

<a id="task-3"></a>

### 3. Interseismic ODE (rate-dependent friction)

> **Status:** ✅ Done · checklist 4/4

- [x] `sys = unicycle.ode.assembleODE(rcv,bigK,law)`: the one-time setup per run, so each RHS call is a single mat-vec with `Kcc`.
  - Reads every parameter from `rcv` and expands scalars to N × 1.
  - Validates sizes, checks `a > b` on creeping patches and that `Vpl` is not all zero, and warns if `Vs` is still the m/yr default.
  - Maps `Vpl`/`Vrake` to `[Vpl·cosd(Vrake), Vpl·sind(Vrake)]` and computes `η = G/(2Vs)`.
  - Returns `iC`, `Kcc`, `taudotLoad`, `y0` and the creeping-patch friction properties.
- [x] `unicycle.ode.odeRateDependent(t,y,sys)` uses `frictionLaw` for `|τ|` and `d|τ|/dζ`. `prepareMesh.m` removed.
- [x] Friction-law options `'asinh'` and `'log'`, with `frictionInverse` for the coseismic update.
- [x] Radiation damping in the inverse. Given a traction `|T|` (e.g. after an imposed earthquake), v solves `|T| = F(v) + η v` instead of `|T| = F(v)`.
  - `g(ζ) = F(ζ) + η v0 e^ζ − |T|` is increasing and convex for both laws, so the root is unique.
  - Newton's method starts from the undamped solution, or from the cap `v = |T|/η` when that is closer and still has g ≥ 0. It converges monotonically from above.
  - `frictionInverse(sys,tau,false)` gives the undamped value.
  - Tested: the residual is 1e-12; damped v is no larger than the undamped v or `|T|/η`; at creep rates the shift matches `ηv/(dF/dζ)`; and in a spring slider started from a damped state, the ODE keeps `F(v) + ηv` equal to the elastic traction.
  - With f0/(a−b) = 120, the two laws agree to 7e-15 MPa over ζ ∈ [−20, 20], and a 300-yr pinned run differs by about 1e-13.
  - The asinh derivative matches finite differences, both inverses round-trip to about 1e-14, and f0/(a−b) = 1200 stays free of overflow.

<sub>[↑ Contents](#contents)</sub>

---

<a id="task-4"></a>

### 4. Interseismic verification

> **Status:** ✅ Done · checklist 5/5

- [x] On the 36-patch test mesh with one pinned patch, 1e4 yr of integration converges to the BEM solution within 1e-3·Vpl (`testMeshSteadyStateMatchesBEM`).
  - f0 was set to 0.01 deliberately, to make rake rotation fast.
  - This shows the loading and the ODE reach the correct steady state. It does not test the slow-rotation hypothesis below, because the same mesh was not run with f0 = 0.6 for comparison.
- [x] With no pinned patches, steady sliding at Vpl is a fixed point in SI units: |ζ̇|, |θ̇| ~ 1e-23 1/s.
- [x] Long-time convergence to the BEM steady state `v_c = K_cc \ (−τ̇_load)`: f0 sweep, `studies/f0ConvergenceStudy.m`.
  - Setup: `generic_sz15_dx20` (750 triangles). Locked zone: rectangle 200 km along strike × 100 km down dip (x ∈ [−110, 90] km, down-dip distance ∈ [40, 140] km, depth ≈ 10–36 km; 100 triangles), edges on grid lines so no lone triangles stay unlocked (`studies/rectangularAsperity.m`). An earlier run with an elliptical asperity (75 triangles, jagged edge) gave the same conclusions. asinh law, v0 = Vpl, a−b = 0.005, σ = 50 MPa, Vpl = 5.4 cm/yr at rake 68°, f0 = 0.1, 0.2, 0.4, 0.6, 1e4 yr from v = Vpl at the plate rake everywhere. `ode45` at 1e-8 takes 3–5 s per run.
  - The first attempt (full 3000-triangle mesh, v0 = 1e-6, f0 from 0.01) was killed after 1 h 44 min still on f0 = 0.01: in the viscous regime dF/dζ ≈ (a−b)σ/160, so the explicit integrator crawls. With v0 = Vpl and f0 ≥ 0.1 the regime is classical RSF (dF/dζ = (a−b)σ).
  - **Result: the hypothesis holds.**
    - Magnitude: the curves for all f0 are identical up to ~100 yr, where the max error first drops to 0.1·Vpl. The speed time scale does not depend on f0. (With the rectangle, f0 = 0.4 and 0.6 rise back above 0.1·Vpl afterwards because of the unconverged rake, and settle below it only at 630 and 1000 yr.)
    - Rake: the error stays at its initial ~52° max (~8.5° rms) for hundreds of years, then decays on a time scale proportional to f0. Max rake error at 1e4 yr: 0.6°, 5.4°, 19°, 30° for f0 = 0.1, 0.2, 0.4, 0.6 (rectangle; ellipse 0.6°, 5.8°, 20°, 30°). Fitting `52°·exp(−t/T)` gives T ≈ 2.2, 4.6, 10.5, 18 kyr, roughly 25–30 kyr × f0.
    - After ~100 yr the residual magnitude error is controlled by the rake error, so it also converges more slowly for larger f0.
    - Consistency: the speed and rake time scales imply the same effective stiffness. (a−b)σ/(k Vpl) = 112 yr and |τ|/(k Vpl) ≈ 2.2 kyr at f0 = 0.1 both give k ≈ 0.04 MPa/m, i.e. a fault-scale mode with wavelength ~ G/k ~ 700 km.
    - The final states for different f0 differ only by the unconverged rake.
    - Where (`studies/output/f0ConvergenceStudyVectors.png`): the rake errors are concentrated between the asperity and the trench, at its two along-strike ends. There the BEM steady state rotates strongly away from the plate rake (≈120° west of the asperity, ≈35° east of it), because the locked zone's stress shadow pushes the creep around it. The runs still sit between the plate rake (68°) and that rotated rake: Δrake < 0 on the west side and > 0 on the east side. Directly up-dip of the asperity |v| < 0.2·Vpl, so the rake there is ill-defined. Elsewhere all runs agree with BEM to within a few degrees.
  - **Plate rake 90° (pure thrust), same setup otherwise** (`plateRake = 90; run('studies/f0ConvergenceStudy.m')`; outputs tagged `_rake68` / `_rake90`):
    - The BEM steady state is nearly antisymmetric about the asperity's centre line: rake ≈ 138° west of it and ≈ 42° east of it next to the trench. It is not exactly symmetric because the box is offset 10 km west.
    - Same qualitative result: magnitude curves are identical for all f0 for ~100 yr, then rake controls everything. Max rake error at 1e4 yr: 0.2°, 3.5°, 14.5°, 22.9° for f0 = 0.1, 0.2, 0.4, 0.6 (68°: 0.6°, 5.4°, 19°, 30°). Fitted e-folding times ≈ 1.8, 3.8, 8.4, 13.5 kyr, about 20–25 % faster than at 68°.
    - Remaining errors are again at the two trench-side ends of the asperity, now antisymmetric (west Δrake < 0, east > 0).
    - The "time until max |v| error < 0.1·Vpl" metric is fragile: after the drop at ~100 yr the rake-driven rebound peaks at ~0.12·Vpl, just above the threshold, so the reported times (0.5–3 kyr) mostly reflect the rake.
  - **Implications.**
    - With realistic f0 ≈ 0.6 the e-folding time of rake errors is ~20 kyr, i.e. many earthquake cycles. Rake rotations from imposed earthquakes will persist across cycles, and a spin-up that starts from v = Vpl at the plate rake would take ~1e5 yr.
    - ~~Start cycle spin-up from the BEM steady state.~~ Superseded by [Task 6](#task-6): with periodic events the cycle-averaged slip is Vpl along the plate rake, so v = Vpl at the plate rake is the better start (2.5× fewer cycles at f0 = 0.6).
    - The periodicity check needs a rake criterion, not just magnitude.
- [x] Compare log and asinh on the same sweep: plate rake 90°, `plateRake = 90; law = 'log'; run('studies/f0ConvergenceStudy.m')` (outputs tagged `_rake90_log`; asinh keeps the untagged-law names).
  - They agree to round-off. At 1e4 yr the final states differ by at most 3e-14·Vpl in |v| and 3e-13° in rake, and the error-vs-time curves differ by at most 1e-13°. The ratio |τ|/(dF/dζ) at Vpl is the same (20, 40, 80, 120), and dF/dζ = (a−b)σ exactly.
  - So every rake-90° number above holds for the log law too: max rake error 0.2°, 3.5°, 14.5°, 22.9° at 1e4 yr.
  - This is expected. The laws differ only for v ≲ v0·exp(−f0/(a−b)) = Vpl·e^−20 or smaller, and the slowest creeping patch runs at 0.06·Vpl.
- [x] Integrator: **`ode45`**, `RelTol = AbsTol = 1e-8` (decided 2026-10-30). One 150-yr postseismic cycle on `generic_sz15_dx20` takes ~0.3 s with v0 = Vpl. It becomes very slow only when the asinh law is in its viscous regime (f(Vpl)/(a−b) ~ O(1)); revisit (`ode15s` with a Jacobian, or the traction-state formulation) only if spin-up runs become slow.

<sub>[↑ Contents](#contents)</sub>

---

<a id="task-5"></a>

### 5. Imposed earthquakes

> **Status:** ✅ Done · checklist 8/8

- [x] Define the asperity with simple geometry and set it as `rcv.pinnedPosition`: `studies/rectangularAsperity.m` (box in along-strike x and down-dip distance, edges on grid lines).
- [x] Coseismic slip on the asperity: `weight·Vpl·T_rec` along the plate rake (`eventStressChange`). `weight` (scalar or per pinned patch) allows an edge taper; with weight < 1 the deficit is only partly released.
- [x] Coseismic stress change on creeping patches: `Δτ = K_cp · slip_P` (two components). `K_cp` is stored by `assembleODE`.
- [x] Coseismic state update in traction space (`applyEvent`):
  - `T⁻` is reconstructed from the state as `|T⁻| = F(v⁻) + ηv⁻` along θ⁻ (`stateTraction`), so it is consistent with the damped ODE.
  - `T⁺ = T⁻ + Δτ`, `θ⁺ = atan2(T⁺_d, T⁺_s)`, and `F(v⁺) + ηv⁺ = |T⁺|` is solved with the damped `frictionInverse`.
  - Slip on creeping patches is continuous across the event.
  - With the log law, `F < 0` below v* makes T⁻ undefined; `stateTraction` raises `unicycle:ode:stateTraction:negativeTraction`.
- [x] Radiation damping makes the postseismic start well posed (no velocity cap needed).
- [x] Example: `studies/postseismicExample.m` → `studies/output/postseismicExample.png`. `generic_sz15_dx20`, rectangular asperity, f0 = 0.6 at v0 = Vpl, a−b = 0.005, rake 68°, T_rec = 150 yr (8.1 m of slip). Starts from `steadyState`, applies one event, integrates 150 yr (0.3 s with `ode45`).
  - Δτ along the plate rake is up to 7 MPa on the patches bordering the asperity and positive around it.
  - v⁺ reaches 0.36 m/s ≈ 2e8 Vpl next to the asperity; damping keeps it far below the undamped value (|T⁺|/η ≈ 7 m/s there).
  - Patches one row away first accelerate as the afterslip front arrives (up-dip centre peaks at ~3e5 Vpl after ~1e-5 yr, i.e. ~5 min; the down-dip probe at ~4e2 Vpl after ~1e-3 yr), then decay roughly as 1/t to below Vpl after ~10 yr.
  - Afterslip after 150 yr is up to 6.8 m (84 % of the coseismic slip) on the ring of patches around the asperity, largest up-dip, between the asperity and the trench.
  - Coseismic rotation of the slip direction is only ≤ 6.5°: the rotation is ≈ Δτ⊥/|T| and |T| ≈ σ f0 = 30 MPa ≫ Δτ. So with f0 = 0.6, f0 matters for rake much more through the slow interseismic relaxation ([Task 4](#task-4)) than through the coseismic jump. Smaller f0 would give larger coseismic rotations.
  - After one cycle the state is far from where it started (max |v/v_ss − 1| = 0.68, max Δrake 5°): the steady state is not the periodic cycle, so spin-up ([Task 6](#task-6)) is needed.

- [x] f0 and rupture-size comparison: `studies/postseismicExample.m` with `f0` ∈ {0.1, 0.6} and `ruptureExtend` ∈ {0, 20 km} (rupture 200 × 100 km = asperity, or 240 × 100 km = 120 % of its area, the extra 20 creeping patches slip 8.1 m coseismically). Outputs `studies/output/postseismicExample_f0_<f0>_ext<km>km.{png,mat}`; comparison `studies/postseismicCompare.m` → `postseismicCompare.png`. asinh law only.

  | rupture | f0 | \|T⁻\| (MPa) | min \|T⁺\| (MPa) | reversed patches | max coseismic rake rotation | runtime | max traction mismatch |
  |---|---|---|---|---|---|---|---|
  | asperity | 0.6 | 29.9 | 29.9 | 0 | 6.5° | 0.3 s | 2.5e-8 MPa |
  | asperity | 0.1 | 4.9 | 4.9 | 0 | 22.5° | 0.2 s | 2.4e-8 MPa |
  | +20 % | 0.6 | 29.9 | 17.8 | 0 | 6.5° | 0.2 s | 2.8e-8 MPa |
  | +20 % | 0.1 | 4.9 | 0.78 | 14 | 180° | 0.3 s | 2.7e-8 MPa |

  - **f0 does not change the speeds.** |v|(t) at every probe and the afterslip after 150 yr (up to ~7 m, ~85 % of the coseismic slip) are nearly identical for f0 = 0.1 and 0.6. Speed is controlled by Δτ/((a−b)σ).
  - **f0 controls the direction.** Coseismic rake rotation ≈ Δτ⊥/|T|, and |T| ≈ σ·f0 is 6× smaller at f0 = 0.1: up to 22.5° vs 6.5° (8° vs 2° along the mid-asperity profile). The rotated rake then persists; after one cycle the rake differs from the steady state by up to 16° (f0 = 0.1) vs 5° (f0 = 0.6).
  - **Rupture onto creeping patches.** The 20 creeping patches inside the rupture get a stress drop of up to 12 MPa along the plate rake.
    - f0 = 0.6: |T| drops from 30 to ≥ 17.8 MPa, never reverses. These patches are effectively locked after the event (v ~ Vpl·e^{−12/0.25}) and only creep again after ~10–100 yr of reloading.
    - f0 = 0.1: the drop exceeds the 4.9 MPa traction, so on 14 patches the traction reverses (|T⁺| down to 0.78 MPa, pointing backward). With asinh this is the linear-viscous regime: backward slip at ~1e-7 Vpl (negligible slip). As the patch is reloaded, the traction rotates through 180° instead of passing through zero (min |T| over the run 0.11 MPa) and returns to forward creep at ~Vpl by ~150 yr.
  - **Numerically stable in all four cases.** `ode45` at 1e-8 needs < 0.3 s, no warnings. The traction implied by the state, `(F(v) + ηv)·[cos θ, sin θ]`, matches the elastic traction `T⁺ + K_cc(slip − slip⁺) + τ̇_load t` to ~1e-7 (a−b)σ over the whole run, including the reversal and the 180° rotation.
  - **Caveat.** The (ζ, θ) state is singular only if the traction passes exactly through zero, i.e. stays exactly anti-parallel while it reloads. Here a perpendicular traction component (oblique rake, rupture ends) keeps it away from zero. A symmetric setup (e.g. a pure-thrust 2D-like rupture with reversal on the symmetry line) could hit the singularity; then switch to the traction-state formulation (Physics section). The log law cannot represent a reversed traction at all (`stateTraction` errors).

- [x] Extended rupture with a prescribed 3 MPa stress drop: same 240 × 100 km rupture, but the slip on the 20 extension patches is solved for (`unicycle.ode.ruptureSlipForStressChange`, a 20 × 20 linear system) so that every one of them gets Δτ·T̂⁻ = −3 MPa (projected on its pre-event traction direction). `dtauTarget = -3` in `postseismicExample.m`; outputs tagged `_drop3`; third row of `postseismicCompare.png`, whose last column now shows |τ|(t).
  - (A first attempt with a cos² taper scaled to keep −3 ≤ Δτ ≤ 10 MPa only made one patch reach −3 MPa; the outer row was loaded by up to +3 MPa. Replaced.)
  - Slip on the extension: 3.5–6.6 m (vs 8.1 m on the asperity). Elsewhere Δτ·T̂⁻ stays in [0, 7] MPa. The slip is the same for both f0 (it depends only on the kernel and the steady-state traction direction).
  - f0 = 0.1: |τ| on the extension drops from ~4.9 to ~1.9 MPa (min |T⁺| 1.68 MPa with the perpendicular component); no reversal, no backward slip. f0 = 0.6: 29.9 → 26.9 MPa.
  - Speed: identical for both f0, v⁺ ≈ Vpl·e^{−3/0.25} ≈ 1e-5 Vpl. The patches are reloaded by the afterslip of their neighbours and by plate loading; v is back to ~Vpl after ~30–100 yr and |τ| back to its pre-event value by ~150 yr.
  - Direction: at f0 = 0.1 the stress drop rotates the traction (and slip) on the extension by up to 28° (west-edge probe: 73° → 96°), relaxing back toward 73° over ~100 yr; at f0 = 0.6 the rotation is ≤ 2°.
  - Same stability: < 0.3 s, traction mismatch ~1e-7 (a−b)σ.

<sub>[↑ Contents](#contents)</sub>

---

<a id="task-6"></a>

### 6. Cycle driver and spin-up

> **Status:** ✅ Done · checklist 5/6

- [x] Loop over cycles (in `unicycle.ode.spinUp`; `runCycles` calls it and then stores one cycle): apply the event to the pre-event state, integrate Trec with `ode45`, repeat. The slip state is reset each cycle, so its end value is the slip of that cycle. Since [Task 13](#task-13) the same loop runs a sequence of events per (super)cycle, a single event being a sequence of one.
- [x] **Slip budget.** `unicycle.ode.periodicEvent` makes the locked patches release exactly Vpl·Trec along the plate rake in every event (`spinUp` warns `unicycle:ode:spinUp:lockedBudget` otherwise). Over one cycle the creeping traction changes by `ΔT = Δτ_event + K_cc Δs_inter + τ̇_load Trec`, so the periodic cycle has a unique slip per cycle `Δs* = slip_C − K_cc \ (Δτ_event + τ̇_load Trec)`, which equals Vpl·Trec on every creeping patch when the locked budget closes (also with coseismic slip on creeping patches). The slip budget therefore closes everywhere once the cycle is periodic.
- [x] Periodicity check, on successive pre-event states: `max |Δs − Δs*|/(Vpl Trec) < 1e-3`, `max |Δζ| < 1e-3`, `max |Δrake| < 1e-3°`, for 3 cycles in a row. Then one full cycle is stored (161 log-spaced outputs from 1e-8 yr to Trec).
- [x] Study: `studies/cycleSpinup.m` (parameters `f0`, `ruptureExtend`, `dtauTarget`, `maxCycles`, `startFrom`) → `studies/output/cycleSpinup_f0_<f0>_ext<km>km[_drop3].{png,mat}`. Same setup as [Task 5](#task-5): `generic_sz15_dx20`, 200 × 100 km asperity, rake 68°, a−b = 0.005, v0 = Vpl, Trec = 150 yr (8.1 m per event), asinh.

  | rupture | f0 | cycles to periodic (from plate rate) | from BEM steady state | wall time | slip/(Vpl T) in stored cycle | min \|T⁺\| over spin-up | reversed |
  |---|---|---|---|---|---|---|---|
  | asperity | 0.6 | 78 (11.7 kyr) | 194 (29 kyr) | 17 s | 0.99996–1.00005 | 29.9 MPa | 0 |
  | asperity | 0.1 | 31 (4.7 kyr) | 43 (6.5 kyr) | 6 s | 1.00000 | 4.94 MPa | 0 |
  | +20 %, 3 MPa drop | 0.6 | 77 (11.6 kyr) | 206 (31 kyr) | 18 s | 0.99996–1.00004 | 26.6 MPa | 0 |
  | +20 %, 3 MPa drop | 0.1 | 31 (4.7 kyr) | 44 (6.6 kyr) | 7 s | 1.00000 | 1.59 MPa | 0 |

  - **The slip budget closes** to < 5e-5 of Vpl·Trec on every creeping patch, and the converged state repeats (Δζ ~ 1e-6, Δrake ~ 1e-3° per cycle).
  - **Convergence is geometric** and set by the slow rake mode: ~4.5× fewer cycles at f0 = 0.1 than at 0.6, consistent with the rake time scale ∝ |τ| ∝ f0 ([Task 4](#task-4)).
  - **The periodic cycle is not the BEM steady state.** With the asperity released every cycle, every creeping patch slips Vpl·Trec along (close to) the plate rake, whereas the steady state around a permanently locked asperity has rakes rotated by up to −58°/+34° and speeds 0.85–1.8 × v_ss. In the periodic cycle the rake stays within ±0.5° of 68° (f0 = 0.6) or within a few degrees except right after the event (f0 = 0.1, e.g. 84° → 68° at the west edge of the extended rupture). Hence the spin-up default is now `sys.y0` (v = Vpl at the plate rake), 1.4–2.5× fewer cycles than starting from the steady state; both reach the same periodic cycle.
  - **f0 = 0.1 is stable through the whole spin-up**, including the extended rupture with a 3 MPa drop: no reversals, min |T⁺| = 1.59 MPa (above the 1 MPa threshold), each cycle ~0.2 s.
  - Late-interseismic coupling (1 − |v|/Vpl just before the event) is high on the ring of creeping patches around the asperity and decays over ~100 km.
- [x] **Spun-up initial condition** (`unicycle.ode.spinUp`). Idea: the largest (characteristic) earthquake dominates the relaxation of the system, so repeating it until periodic primes slip rates and tractions; the cycle map is contractive (creep diffuses the memory of the initial condition). Recorded instrumental-era events and creep are then imposed on this state with `applyEvent` and `odeRateDependent` using the same `sys`.
  - Returns the state `tSince` after the last characteristic event, or after the start of the last supercycle for a sequence (slip zeroed), `info.slipSince`, `info.yPre`, `info.history`, and `info.contraction` (geometric decrease of the slip error per cycle).
  - 20 km mesh, rupture = asperity, Trec = 150 yr: contraction per cycle 0.79 (f0 = 0.1, 31 cycles) and 0.96 (f0 = 0.6, 78 cycles). Starting from the locked-asperity steady state instead gives the same primed state (at tSince = 120 yr: max Δζ 1e-5 / 2e-4, max Δrake 0.003° / 0.04° for f0 = 0.1 / 0.6) after 43 / 194 cycles.
  - The stopping rule bounds the change per cycle, not the remaining error; with contraction r the remaining error is ≈ change·r/(1−r), i.e. ~25× the last change at f0 = 0.6. Tighten `tolRake`/`tolZeta` if the primed rake must be accurate to better than a few hundredths of a degree.
- [ ] Save `t`, slip and velocity on **all** patches (pinned patches: v = 0 between events, coseismic slip at events) — move to [Task 7](#task-7) together with surface displacements.

<sub>[↑ Contents](#contents)</sub>

---

<a id="task-7"></a>

### 7. Outputs

> **Status:** 🔄 In progress · checklist 1/2

- [x] Surface velocity and displacement time series at stations.
  - `G = unicycle.greens.computeSurfaceKernels(rcv,xObs)`: `rcv.displacementKernels(xObs,3)` as one (3D × 2N) matrix `[G_strike, G_dip]`, rows `[e1 n1 u1 e2 ...]`; xObs in `[east, north, up]` like `rcv.xc`.
  - `[vSurf,uSurf] = unicycle.ode.surfaceDeformation(G,sys,t,Y,opts)`: velocity `G_C v_C(t)` (pinned patches v = 0), displacement `G_C slip_C(t) + G_P slip_P`, both (nt × D × 3) [east, north, up].
    - The displacement is evaluated from the slip state, which `ode45` integrates together with ζ and θ, so it is the exact time integral of the velocity. A trapezoidal integral of the velocity is not reliable right after an event (v ~ 1e8 Vpl decaying within seconds): it was off by ~5 % on a grid starting at 1e-6 yr.
    - `opts.removeLongTerm` subtracts `G Vpl` (and `G Vpl t`), giving the signal relative to the long-term motion, i.e. `G (v − Vpl)`.
    - `opts.pinnedSlip` (e.g. `out.event.slipP`) adds the locked-zone coseismic slip; `opts.reference = 'first'` makes the displacement relative to the first time.
  - Tests (`testSurfaceDeformation`): kernel layout, thrust moves the hanging wall toward the trench and up, positive strike slip moves it along the strike vector; steady sliding at Vpl is removed exactly by `removeLongTerm`; a station above a locked patch moves landward relative to the long-term motion (back slip); displacement = time integral of velocity after an event; coseismic step = G·slip.
  - Study `studies/surfaceTimeseries.m` (trench-normal profile at east = −10 km, f0 = 0.1 and 0.6, periodic cycle): relative to the long-term motion the displacement is a sawtooth — a coseismic step (up to 3.6 m), fast postseismic transient over ~1–10 yr, then linear interseismic recovery — and returns to zero after one cycle to 2e-4 m (f0 = 0.6) / 3e-5 m (f0 = 0.1), i.e. the slip budget closes at the surface too. The two f0 give nearly identical surface signals (f0 mainly changes rake, which the surface sees only weakly). Stations exactly on the fault trace (north = 0) are singular; keep them off it.
- [ ] Plots: map views of v/Vpl and rake at selected times, time series at selected patches, and the cumulative slip budget (total slip per cycle ≈ Vpl·T_rec everywhere).

<sub>[↑ Contents](#contents)</sub>

---

<a id="task-7b"></a>

### 7b. Instrumental era

> **Status:** ⬜ Not started · checklist 0/1

- [ ] Driver that starts from `spinUp(..., tSince = time since the last characteristic event)` and imposes the recorded sequence (the event loop now exists: `spinUp(sys,seq,struct('repeat',false,'y0',y))`, see [Task 13](#task-13)): coseismic slip of recorded events via `applyEvent` (traction change from their slip distributions, which may include creeping patches), creep in between with `odeRateDependent`, loading continuing with the locked patches pinned.

<sub>[↑ Contents](#contents)</sub>

---

<a id="task-8"></a>

### 8. Later / optional

> **Status:** 💤 Backlog · checklist 0/5

- [ ] Multiple asperities with different recurrence times and non-uniform coseismic slip.
- [ ] Alternative laws behind the same interface: the regularized sinh form and power-law viscous `v = τⁿ/η`.
- [ ] Apply the model to the Mexico mesh (`fault/mexico_mt_xyz`) with Euler-pole loading. This needs the outline and lat/lon files that are not in this folder.
- [ ] Performance for larger meshes: compressed kernels such as H-matrices, or a GPU.
- [ ] Benchmark against seqeas in an equivalent configuration.

<sub>[↑ Contents](#contents)</sub>

---

<a id="task-9"></a>

### 9. Linear viscous friction: analytic cycles and inversion

> **Status:** 🔄 In progress · checklist 3/7

Goal: an effective linear viscous friction law on creeping patches, so the cycle ODE is linear and has a closed-form solution. That gives exact spin-up and time evolution, and makes the effective viscosity structure the unknown of an inverse problem fitted to geodetic time series.

**Formulation** (creeping patches; pinned patches keep v = 0 between events)
- Friction: `τ = τ0 + D v`. Here `D` is a 2Nc × 2Nc positive (block-)diagonal matrix of effective viscosities, in units of MPa·s/m (traction per unit slip rate, not Pa·s).
  - **Decision: isotropic only.** `D = diag(η, η)`, one scalar η per patch acting on both slip components.
  - Use a different symbol from the radiation damping `sys.eta = G/(2Vs)`. Damping simply adds to D and is negligible: 5 MPa·s/m against η ~ 1e8 (below).
- Elasticity, as now: `τ̇ = K_cc v + τ̇_load`, where `τ̇_load = −K(c,:)·Vpl` keeps the pinned patches at v = 0.
- Combined: `D v̇ = K_cc v + τ̇_load`, i.e. `v̇ = A v + D⁻¹τ̇_load` with `A = D⁻¹ K_cc`.
- τ0 drops out of the rate form entirely. There is no counterpart of f0, and no collinearity assumption is needed beyond the choice of D.
- Steady state: `v* = −K_cc⁻¹ τ̇_load`. This is the existing BEM steady state (`steadyState`), and it is independent of D.
- Transient: `u = v − v*` obeys `u̇ = A u`, so `u(t) = exp(A t) u0`. Slip has a closed form too: `s(t) = v* t + A⁻¹ (exp(A t) − I) u0`. Surface displacement is `G_surf·s(t)`, using the existing kernels.
- Event with traction change Δτ: instantaneous velocity jump `u⁺ = u⁻ + D⁻¹Δτ`.
  - `Δτ = K_cp slip_P + K_cc slip_C`, the same as `eventStressChange`.
  - Slip on creeping rupture patches, if any, jumps by `slip_C`.

#### Exact periodic spin-up

- Let `E = exp(A T)` and `b = D⁻¹Δτ`. With `p` the pre-event state, the periodic condition `p = E (p + b)` gives `p = (I − E)⁻¹ E b`.
  - In the eigenbasis of A this is `p_i = b_i e^{λ_i T}/(1 − e^{λ_i T})`.
  - It exists whenever all `Re λ_i < 0`.
- State at `tSince` after the last event: `u = exp(A tSince)(p + b)`.
- **Slip budget closes exactly, independently of D.** The slip per cycle is `v* T − A⁻¹ b + slip_C = −K_cc⁻¹(Δτ + τ̇_load T) + slip_C`, the same Δs* as in `spinUp`, which equals Vpl·T when the locked budget closes. So long-term slip and coupling carry no information on D. Only the transients do: the timing of postseismic and interseismic rates.
- This is a closed-form counterpart of `spinUp` / `runCycles`. It needs no time stepping and no convergence loop.

#### Is K_cc negative definite?

- In the continuum, the elastic energy of any nonzero slip on a bounded fault in a half-space is positive. So the operator is negative definite, not just semi-definite. It would be only semi-definite with a zero mode for an infinitely long fault, or for rigid translation of an unbounded one. Its smallest |eigenvalue| scales as G/L_fault: [Task 4](#task-4) found an effective k ≈ 0.04 MPa/m for the slowest mode of the 20 km mesh.
- Removing pinned patches (K_cc is a principal submatrix) keeps definiteness.
- **Caveat: the discrete K from centroid collocation is not exactly symmetric.** Reciprocity holds for patch-integrated quantities, not for centroid values. So `A = D⁻¹K_cc` may have complex eigenvalues and be non-normal. For symmetric K, A is similar to `D^{-1/2} K D^{-1/2}`, whose eigenvalues are real and negative.

#### Relation to rate-and-state (what the linear model can and cannot do)

- **Tangent viscosity of asinh RSF** around a reference speed: `dτ/dv = (a−b)σ / v_ref` along the slip direction.
  - At v_ref = Vpl: η ≈ 0.25 MPa / 1.7e-9 m/s ≈ 1.5e8 MPa·s/m.
  - The resulting relaxation time is η/k ≈ (a−b)σ/(k Vpl) ≈ 110 yr, matching the RSF speed time scale in [Task 4](#task-4).
- **Anisotropy.** RSF resists rotation of slip with `|τ|/v`, not `(a−b)σ/v`. This is the f0 effect on rake ([Task 4](#task-4)). A faithful linearization therefore uses a 2 × 2 tensor per patch in the frame of the reference slip direction: `η_∥ = (a−b)σ/v_ref` and `η_⊥ = |τ_ref|/v_ref = η_∥ (f0/(a−b) + ln(v_ref/v0))`. That is still linear (block-diagonal D), and it brings back the rake/f0 behaviour. An isotropic η is equivalent to `|τ| = (a−b)σ`, i.e. the low-f regime. **We use isotropic η only**, so the linear model does not carry the f0 effect on rake; keep this in mind when comparing with RSF runs at f0 = 0.6 ([Task 9.4](#task-9)).
- **Limits.** A linear model relaxes as a sum of exponentials. RSF after a large stress step (Δτ ≫ (a−b)σ) gives ~1/t (Omori-like) decay from v ~ 1e8 Vpl. A single constant D per patch cannot reproduce the early postseismic phase of large events. For the inverse problem, decide whether D is a fixed per-patch property, or is allowed to vary with time window or to be re-linearized.

#### Tasks

- [x] 9.1 Properties of K (both generic meshes, rectangular asperity pinned):

  | mesh | Nc | ‖K−Kᵀ‖/‖K‖ | eig(K_cc) real parts | max\|Im\|/max\|Re\| | cond(V) | eig time |
  |---|---|---|---|---|---|---|
  | `generic_sz15_dx20` | 650 | 2.4 % | −3.6 … −0.12 MPa/m | 4e-4 | 64 | 0.4 s |
  | `generic_sz15` | 2600 | 2.1 % | −7.2 … −0.13 MPa/m | 6e-4 | 400 | 13 s |

  - All eigenvalues have negative real parts; the symmetric part is negative definite over the same range. The slowest mode, k ≈ 0.12 MPa/m, is fault-scale and nearly the same on both meshes; the fastest scales with 1/cell size.
  - The imaginary parts are negligible and V is well conditioned, so **eigendecomposition** is used (no symmetrization). `linvisAssemble` warns if Re λ ≥ 0 or cond(V) > 1e10.
  - With uniform η = (a−b)σ/Vpl ≈ 1.5e8 MPa·s/m, the relaxation times −1/Re λ span ~1–40 yr (0.6–37 yr on the 10 km mesh).
- [x] 9.2 Forward model, analytic (`+unicycle/+ode/`, same conventions as the RSF code: [strike; dip] blocks on creeping patches):
  - `lin = linvisAssemble(rcv,bigK,eta)`: partition (`iC`, `iP`, `Kcc`, `Kcp`, `taudotLoad`, `Vpl`, so `periodicEvent`/`eventStressChange`/`surfaceDeformation` accept `lin`), `vStar`, eigendecomposition of `A = D⁻¹K_cc`. η scalar, Nc × 1 or N × 1. Needs only `Vpl`, `Vrake`, `pinnedPosition` on `rcv`. **Decision:** η is an argument of `linvisAssemble`, not a property on `rcv` (the receiver class has no viscosity property, and passing it keeps inversion loops over trial η simple). This is the one exception to the rule that physical properties live on the mesh object.
  - `[v,s] = linvisEvolve(lin,v0,t)`: `v = v* + V(e^{λt}.*c)`, `s = v* t + V((e^{λt}−1)./λ.*c)`, `c = V⁻¹(v0 − v*)`; `e^z − 1` evaluated as `2e^{z/2} sinh(z/2)` to avoid cancellation.
  - `v = linvisEvent(lin,v,dtau)`: `v + Δτ/η`.
  - `[v,info] = linvisSpinUp(lin,event,tSince)`: exact periodic cycle, `p_i = b_i e^{λ_i T}/(1 − e^{λ_i T})`; returns the slip rate at tSince, `vPre`, `vPost`, `slipCycle`, `slipSince`, and `spectralRadius` = max |e^{λT}|.
  - `[v,info] = linvisSpinUp(lin,seq,tOut,opts)`: the same function, given a sequence (struct with `T` and `events`) instead of an event, gives the exact periodic cycle for a sequence of events (times t_k in [0,T), slip on locked and/or creeping patches) repeating every supercycle T. With p = u(0⁻) = v(0⁻) − v*, `p = (I − E(T))⁻¹ Σ_k E(T − t_k) b_k`, i.e. the sum of single-event periodic cycles of period T shifted to t_k (periodic Green's function `G_T(s) = E(s mod T)(I − E(T))⁻¹`). The constant loading is carried by v*: slip is `v* t + ∫u + Σ slip_C`. Returns v at tOut (event at t_k included at t_k; t = T is the state before the next supercycle), `vPre`/`vPost` at each event, slip, `slipCycle`, `budgetErrP`, `budgetErrC`. **Slip budget:** Δs_C = Vpl_C T − K_cc⁻¹K_cp(Σ slip_P − Vpl_P T), so closing the locked budget closes it on every creeping patch, wherever the events rupture; the locked budget is checked (error above 1e-6 |Vpl| T; `opts.checkBudget = false` to skip). **Decision:** one function for both; `linvisSpinUp` recognises a sequence by its `events` field and treats an event (with `Trec`) as a sequence of one event at t = 0 (budget checked when `event.slipP` is given). `sequenceEvent(sys,t,slip)` builds an event from N × 2 slip (Δτ = K_cp slip_P + K_cc slip_C). `linvisSpinUpNumerical` also accepts a sequence.
  - Surface: `surfaceDeformation(G,lin,t,struct('v',v,'s',s),opts)` (it now also accepts slip rate and slip directly).
- [x] 9.3 Verification (`tests/testLinvis.m`, 8 tests): partition and v* match `assembleODE`/`steadyState`; eigendecomposition reproduces A; analytic v and s match `ode45` for depth-dependent η; v* is a fixed point; η → cη equals t → t/c; event jump; closed-form periodic state returns to itself after one cycle and equals 200 iterations of the cycle map; slip per cycle = Vpl·T for two different η, with and without slip on creeping rupture patches; output feeds `surfaceDeformation` and the relative displacement closes after a cycle.
  - **Analytic vs numerically spun-up cycle.** `odeLinearViscous(t,y,lin)` is the time-stepping ODE (state `[slip_s; slip_d; v_s; v_d]`, `η v̇ = K_cc v + τ̇_load`, `ṡ = v`), and `linvisSpinUpNumerical(lin,event,opts)` repeats event + `ode45` over Trec until periodic (slip error and change in pre-event v < 1e-8, 3 cycles in a row). `testAnalyticSpinUpMatchesNumericalSpinUp` compares pre/post-event slip rate, slip per cycle and the whole periodic trajectory for depth-dependent η, with coseismic slip on creeping patches, and for 30× larger η: agreement to 1e-7 Vpl. On the 20 km mesh (`studies/linvisCycle.m`) the numerical spin-up needs 8 cycles (~1 s) and matches the closed form to 2e-11 Vpl (uniform η) and 2e-9 Vpl (η ×10 with depth) in pre-event slip rate, and to ~1e-11 Vpl·T in slip per cycle.
  - **Extended rupture** (`testExtendedRuptureAnalyticMatchesNumerical`): rupture 20 % larger in area than the locked zone, uniform slip Vpl·Trec on the locked patches and 0.5·Vpl·Trec on the creeping patches of the extension. Analytic and numerical periodic cycles agree to 1e-7 Vpl (depth-dependent η and 10× that), and the slip per cycle is still Vpl·Trec everywhere. On the 20 km mesh (20 km extension on each side, exactly 20 % of the locked area): Δτ along slip on the extension is −5.3 to +0.9 MPa (with half the slip, loading by the adjacent locked-zone slip wins on a few patches), the numerical spin-up takes 8 cycles and matches the closed form to 2e-11 Vpl (uniform η) and 8e-10 Vpl (η ×10 with depth), and slip/(Vpl·T) = 1.000000.
  - **Event sequences** (`testSequenceAnalyticMatchesNumerical`, `testSequenceIsSumOfShiftedSingleEventCycles`, `testSequenceSlipBudgetCheck`): a 120 yr supercycle of four events, given out of order — west half of the locked zone with Vpl·T extending onto the nearest creeping patches, east half in two events (0.6 and 0.4 Vpl·T, the second extending onto creeping patches), and an event on deep creeping patches only. Closed form vs numerical spin-up agree to 1e-7 Vpl in v before and after every event and in slip per supercycle (depth-dependent η and 10× that); slip per supercycle = Vpl·T on every creeping patch. The sequence solution equals the sum of shifted single-event cycles, and one event at t_k equals the `linvisSpinUp` cycle shifted by t_k. Dropping the west-half event raises the budget error, and with the check off the creeping slip matches `Vpl_C T − K_cc⁻¹K_cp(Σ slip_P − Vpl_P T)`. For η → ∞ the slip rate is constant, Vpl − Σ slip_C/T.
  - Viscosity limits (`testViscosityLimits`): η → ∞ gives steady creep at Vpl through the whole cycle (`v₀ − v* = (I − E)⁻¹b → −K_cc⁻¹Δτ/T`, the events are averaged out); η → 0 gives instant relaxation to v*. Both close the slip budget. The spin-up only needs `I − E` invertible (all Re λ < 0); `e^{λT} − 1` is computed by `linvisExpm1` (sinh form for |z| < 1, direct otherwise), because `1 − exp(λT)` loses all digits for |λT| ≲ 1e-8 and the sinh form gives 0·∞ for λT → −∞.
  - Only the loading term needs A invertible (v* = −A⁻¹c). If near-zero modes ever appear, use `v₀ = (I − E)⁻¹[T φ₁(AT) c + b]`; a periodic solution then requires `T c̃ᵢ + b̃ᵢ = 0` in that mode (slip-budget closure per mode).
  - Not yet done: the small-perturbation limit against `odeRateDependent` with the tangent viscosity.
- [~] 9.4 Linear vs RSF, first look (`studies/linvisCycle.m` → `studies/output/linvisCycle.png`; 20 km mesh, rectangular asperity, rake 68°, Trec = 150 yr; RSF asinh f0 = 0.1; linear with uniform η0 = (a−b)σ/Vpl and with η increasing ×10 with depth around η0):
  - Cost: RSF spin-up 31 cycles in 5.8 s; linear assemble 0.7 s (eigendecomposition), then exact spin-up + 162 output times in ~0.015 s.
  - Slip budget: slip/(Vpl·T) = 1.000000 on every creeping patch for both η.
  - Memory: the linear cycle map contracts by 0.02–0.04 per cycle (all relaxation times ≪ Trec), so the linear model forgets its initial condition within one cycle. RSF contracts by only 0.79 per cycle at f0 = 0.1 because of its slow rake mode, which the isotropic linear model does not have.
  - Early postseismic: the linear jump is Δτ/η ≈ 7 MPa / 1.5e8 ≈ 30 Vpl, versus 1e5–1e8 Vpl for RSF. With η tuned to plate rate the linear model is far too slow in the first year and spreads afterslip over ~10 yr; the surface displacements differ most in the first ~10–20 yr near the asperity, and agree closely in the late interseismic period (where both are controlled by v* and the same budget).
  - Next: fit η (uniform first) to the RSF synthetics over different windows (e.g. 1–10 yr, 10–150 yr) to quantify which windows a linear model can represent. Planned as [Task 12](#task-12).
- [ ] 9.5 Inverse problem design:
  - **Data:** station displacement time series over a window, relative to long-term motion.
  - **Unknowns:** the coefficients of a low-rank viscosity model (see 9.5a), possibly tSince or T. η is not free at every patch.
  - **Known:** Δτ from the characteristic event, locked region, Vpl.
  - **Structure:** given D, the prediction is linear in b (and in u0), so amplitude-type unknowns can be eliminated by linear least squares (variable projection).
  - **Derivatives with respect to D:** eigenvalue/eigenvector perturbation, or the Fréchet derivative of `expm`; or derivative-free/MCMC for a few parameters.
- [ ] 9.5a Low-rank viscosity parameterization.
  - Work with `m = log η` (keeps η > 0; viscosity varies over orders of magnitude).
  - Prior covariance on patch centroids, e.g. Matérn, `C_ij = σ_m² · Matérn_ν(r_ij/ℓ)`; ν = 1/2 is the exponential kernel `σ_m² exp(−r_ij/ℓ)`. r_ij is the distance between centroids: Euclidean in 3D for now, along-fault distance if the mesh is strongly curved.
  - Truncated Karhunen–Loève expansion: `m = m̄ + Φ_r Λ_r^{1/2} c`, where (Φ_r, Λ_r) are the r leading eigenpairs of C, and `c ~ N(0, I)` are the unknowns. r is set by the fraction of prior variance kept, e.g. 95–99 %; for ℓ ~ 50–100 km on a 500 × 300 km fault, expect tens of modes, not hundreds.
  - Hyperparameters: mean m̄ (e.g. log of the RSF tangent value (a−b)σ/Vpl), σ_m, ℓ, ν, r. Fix them first, then test their sensitivity or treat them hierarchically.
  - The modes depend only on the mesh and hyperparameters: compute once and cache them, like the kernels. The covariance is on all creeping patches (Nc × Nc), so its eigendecomposition is cheap on both meshes.
  - Alternative basis to compare: smooth eigenmodes of the elastic operator K_cc or of a mesh Laplacian.
  - Gradients: `∂d/∂c = (∂d/∂m) · Φ_r Λ_r^{1/2}`, so only r directional derivatives of the forward model are needed. With r small, finite differences or complex-step on c are affordable; otherwise use the eigen-perturbation/Fréchet derivative of expm.
- [ ] 9.6 Identifiability on synthetics. Generate data with the linear model and with the RSF model, then invert. Which KL modes (and hence which parts of η) are resolved by a given station geometry and window length (near-asperity vs distant, shallow vs deep), and how do trade-offs with tSince, T and Δτ look?
- [ ] 9.7 Cost and scaling: one decomposition per D, O((2Nc)³). That is trivial for 750 patches (1300 × 1300) and a few seconds for 3000 patches (~5400 × 5400); it governs the cost of each inversion iteration.

#### Open questions

- ~~Isotropic η or the RSF-derived tensor?~~ Isotropic (decided).
- ~~One D for the whole cycle?~~ Different D per time window (decided, see [Task 10](#task-10)).
- ~~Parameterization of D~~: low-rank KL modes of a Matérn covariance on log η (decided). Still open: values of ℓ, ν, σ_m and r, and whether to estimate them.
- Which data: GNSS displacement time series only, or also InSAR / long-term rates?

<sub>[↑ Contents](#contents)</sub>

---

<a id="task-10"></a>

### 10. Time-windowed effective viscosity: inversion design

> **Status:** 📝 Planning · checklist 0/5

Idea: effective viscosity is not one constant. It is large during the late interseismic period and small right after an earthquake, increasing with time since the event. The plan is a piecewise-constant η in time, with slip rate continuous across window boundaries. Every window is solved with the [Task 9](#task-9) closed form.

#### Windows

- **W0, interseismic:** `η_0` (large), used in the exact spin-up (`linvisSpinUp`) with the characteristic event (Δτ_char, Trec). It gives the slip rate at `tSince` = time from the last characteristic event to the modern event. This is fitted to late-interseismic surface velocities (relative to long-term motion).
- **Coseismic:** the slip-rate jump uses the late-interseismic slip rate from W0 and the coseismic traction change of the modern event, `Δτ_obs = K_cp slip_P + K_cc slip_C` (from its slip distribution, which may include creeping patches): `v(0⁺) = v_W0(tSince) + Δτ_obs/η_1`. The jump is divided by the viscosity of the **first postseismic window**, η_1.
- **W1 … WK, postseismic:** windows `[t_{k−1}, t_k]` after the modern event (e.g. 0–1 month, 1–3 months, 3–12 months, 1–3 yr, …), each with its own `η_k`. Within a window, `v(t) = v* + V_k e^{Λ_k (t − t_{k−1})} V_k⁻¹ (v(t_{k−1}) − v*)`. The slip-rate field is continuous on the whole mesh at every boundary, `v(t_k⁻) = v(t_k⁺)`.
- **Shared steady state:** `v* = −K_cc⁻¹τ̇_load` is the same for all windows (independent of η), so only the decaying part changes from window to window. Slip and surface displacement per window are closed form (`linvisEvolve`) and accumulate across windows; add the coseismic step `G_P slip_P + G_C slip_C`.

#### Physics notes

- Continuity of v across a change in η means the traction offset τ0 jumps (`τ0_k = τ − η_k v`). That is acceptable for an effective model: τ0 never enters the rate equation.
- RSF tangent viscosity `η ≈ (a−b)σ/v` predicts exactly this behaviour: small η when v is high right after the event, growing as v decays, and largest late in the cycle. It could serve as a physics-based prior or consistency check: `η_k · v̄_k ≈` const (≈ (a−b)σ) on each patch.
- **Interseismic η sets the late-interseismic coupling pattern.** Its limits ([Task 9](#task-9)) are η → ∞, steady creep at Vpl around the locked zone (no stress shadow), and η → 0, the full stress shadow v* right after relaxation. The fitted η_0 places the late-interseismic state between these.
- **η_0 is a target of the inversion, not a nuisance parameter (decided).** We do not observe a full earthquake cycle, only a quasi-steady late-interseismic velocity field. η_0 is defined as the effective viscosity that, used throughout a periodic characteristic cycle (`linvisSpinUp` with η_0 only), reproduces that quasi-steady state at `tSince`. It is a value we interpret, so the spin-up deliberately uses one η for the whole cycle rather than a window schedule.
- **What η_0 can be identified against.** For uniform η the model has the scaling `v(η, Trec, tSince) = v(cη, cTrec, c·tSince)` (`testViscosityScalesTime`). So a single velocity snapshot constrains η_0 only relative to the assumed Trec and tSince: their values (paleoseismic/historical) must be fixed or carried as uncertain, and the reported η_0 is conditional on them. Spatial variations of η_0 (the KL modes) are constrained by the shape of the coupling pattern, between the limits v* (small η) and Vpl (large η).

#### Inverse problem

- **Unknowns:** log η per window as truncated KL fields ([Task 9.5a](#task-9)), `m_k = m̄_k + Φ_r Λ_r^{1/2} c_k`, k = 0…K, with **r = 2–3 modes per window, including the interseismic window W0** (decided). The truncation is the regularization. W0 needs several modes because which parts of the mesh sit in the locked zone's stress shadow (small η, v → v*) and which are stiff enough to creep at full plate rate (large η, v → Vpl) is not known a priori. So the unknown vector is `[c_0; c_1; …; c_K]`, about 3(K+1) numbers; possibly tSince/Trec as well.
- **Data:** late-interseismic GNSS velocities, which constrain η_0 through the spin-up (W0), and postseismic GNSS displacement time series (W1…WK), all relative to long-term motion. The W0 fit can be done first and on its own; it fixes the pre-event slip rate that the postseismic windows start from.
- **Forward cost:** one eigendecomposition per window (0.7 s on the 20 km mesh, ~13 s on the 10 km mesh), then cheap closed-form propagation.
- **Jacobian structure ("triangular in time").** Rows are grouped by data window, columns by the coefficient block c_k of each window. Data in window k depend on c_0…c_k (earlier windows enter through the slip rate carried across window edges and, for W1, through the pre-event rate and η_1 in the coseismic jump), but never on c_{k+1}…c_K. So `∂d_k/∂c_j = 0` for j > k and the Jacobian is block lower-triangular. In practice:
  - perturbing c_k needs a new eigendecomposition for window k only, then re-propagation from window k onward; windows before k are unchanged;
  - with ~3(K+1) unknowns, finite differences (or complex step) on c are affordable, and so is sampling (MCMC);
  - a sequential fit (W0 alone, then W1 given W0, …) is a natural starting model for the joint inversion.
- **Regularization:** spatial smoothness from the KL prior; in time, penalize differences of log η between consecutive windows, and optionally enforce η_1 ≤ η_2 ≤ … ≤ η_K ≤ η_0 (monotone increase).

#### Tasks

- [ ] 10.1 Forward chain: given (η_0, Δτ_char, Trec, tSince), (Δτ_obs, slip_C) and (η_1…η_K, window edges), return slip rate and slip on the fault and surface velocity/displacement at stations over all windows. Reuse `linvisAssemble`/`linvisSpinUp`/`linvisEvolve`/`linvisEvent`/`surfaceDeformation`; one `lin` per window, sharing kernels and v*.
- [ ] 10.2 Tests: continuity of v at window edges; with all η_k equal, identical to a single `linvisEvolve`; each window matches `ode45` of `odeLinearViscous` with that window's η; slip and displacement accumulate correctly across windows.
- [ ] 10.3 (Planned in detail as [Task 12](#task-12).) Synthetic test against RSF: generate postseismic time series with the asinh RSF model ([Task 5](#task-5) setup), fit η_k per window with the linear chain, and compare the fitted η_k with the RSF tangent `(a−b)σ/v̄_k`. How many windows does a log-like RSF decay need?
- [ ] 10.3b Interseismic-only inversion (W0): fit η_0 (KL field) to late-interseismic velocities on synthetics; show the η_0–Trec–tSince trade-off and how well the spatial pattern of η_0 is recovered.
- [ ] 10.4 Inversion prototype on synthetics with KL-parameterized η_k ([Task 9.5a](#task-9)), including resolution/identifiability of each window's η.

#### Open questions

- Window edges: fixed and log-spaced in time since the event, or part of the inversion?
- ~~Spin-up with η_0 only, or with a window schedule?~~ η_0 only: η_0 is the effective viscosity of the quasi-steady late-interseismic state and is itself a result (decided).
- Trec and tSince for the characteristic event: fixed, or inverted with priors (η_0 trades off with them)?
- Independent η fields per window (own c_k, r = 2–3 each), or a shared spatial pattern scaled per window? Current plan: independent c_k per window.
- Does the modern event replace the characteristic event in the cycle (tSince measured from it afterwards), or is it superposed on the characteristic cycle?

<sub>[↑ Contents](#contents)</sub>

---

<a id="task-11"></a>

### 11. Inverse package: design

> **Status:** 📝 Planning · checklist 1/5

**Location.** All inversion code goes in a new package folder **`+unicycle/+inverse/`** (decided), called as `unicycle.inverse.*`. **Plain functions, no classes (decided):** each function takes the objects and data it needs (`rcv`/`lin`, kernels, KL basis, data structure, windows) and returns plain structures. This follows the convention that machinery lives in `+unicycle` packages (ODE in `+ode`, kernels in `+greens`). The forward pieces it calls stay in `+ode` (linvis*, surfaceDeformation). Tests in `tests/` as for every module.

#### User-supplied

- Time windows: edges after the modern event (and tSince, Trec, Δτ_char for W0).
- Data and their uncertainties (GNSS velocities for W0; displacement time series for W1…WK), station coordinates.
- Optionally the correlation length ℓ and the KL kernel (Matérn ν; default ν = 1/2, exponential).

#### Data with gaps and mixed sampling

- Time series will have gaps and different sampling rates between stations (and possibly within a station). Store observations in long format, one row per observation: station index, component (E/N/U), time, value, σ (and optionally a data covariance per station). There is no common time grid: each station keeps its own times.
- The forward model is closed form, so it is evaluated at exactly the observed times; no interpolation and no time stepping. Per window, collect the unique observation times in that window, evaluate slip there with `linvisEvolve`, and map to each station with its kernel rows.
- Displacement time series have an arbitrary reference position per station and component, so add one offset per station-component as a nuisance parameter. It enters linearly and is eliminated by least squares inside the misfit (variable projection). If a station has equipment offsets or long gaps, allow one offset per segment.
- Interseismic data (W0) are velocities, not positions, and need no offsets.
- Windows contain whatever observations fall in them; a window with few or no data is simply less constrained (fewer modes selected, see below).

#### Parameters per window: r + 1 (decided)

- Each window k has **r_k + 1 unknowns**: the uniform level `m̄_k = log η̄_k` (log of a uniform viscosity) and r_k KL coefficients `c_k`, so `log η_k = m̄_k + Phi(:,1:r_k) c_k`. The KL modes are zero-mean (`klBasis` removes the mean), so they describe log-space spatial variations *about* the uniform viscosity, and m̄_k does not trade off with them.
- **Workflow for every window:**
  1. First estimate the best-fit **uniform** viscosity (r_k = 0, one parameter). This is a 1-D problem (a line search on log η̄), gives the starting model, and is the r = 0 baseline in the model selection.
  2. Then estimate m̄_k **jointly** with the r_k mode coefficients (r_k + 1 parameters), starting from the uniform fit with c_k = 0. m̄_k is re-estimated, not held at the uniform value, because spatial variations change the best uniform level (in log space, its mean under the fitted pattern).
- The mean sets the relaxation time scale and is the best-constrained number, so it is always a free parameter, not a prior mean.
- KL modes `Φ_r` of the covariance `exp(−r_ij/ℓ)` (or Matérn) on creeping-patch centroids, normalized by `Λ_r^{1/2}`; computed once per (mesh, ℓ, ν) and cached.

#### Correlation length ℓ: one value for all windows (decided)

- ℓ only defines a common basis: the same KL modes `Φ` serve every window. Windows differ in how many of these nested modes they use and in their coefficients.
- Default when not supplied: **ℓ = L/2**, with L the largest extent of the creeping-patch centroids. Use the largest principal extent (PCA of the centroids), so curved or rotated meshes work too. For the generic meshes that is ~500 km along strike → ℓ ≈ 250 km.
- With only 2–3 modes, ℓ mainly sets the shape of those modes. For ℓ ≳ L the leading modes tend to a near-constant field and smooth gradients along the long and short dimensions, which is what "stress shadow near the asperity vs full plate rate far away" needs. For ℓ ≪ L the variance spreads over many modes and 2–3 modes capture little of the prior. So ℓ ~ L/2 is a sensible default.
- Check the sensitivity by treating ℓ ∈ {L/4, L/2, L} as candidate models in the same information-criterion comparison as r.

#### Number of KL modes: information criterion

- For each candidate r (and ℓ), fit by least squares and compute `BIC = p ln n + χ²_min` (Gaussian errors, so −2 ln L̂ = χ² + const), with p the number of free parameters (means + modes over all windows, plus tSince/Trec if free) and n the number of data. Pick the minimum; also report AIC/AICc (less conservative, better for small n).
- Search strategy: modes are nested (r = 0, 1, 2, 3 …), with r = 0 the uniform-viscosity model (1 parameter) and r modes giving r + 1 parameters per window. Select r_0 for W0 on the interseismic data alone, then r_k window by window given the earlier windows (the triangular structure), then a final joint fit. A full combinatorial search over all windows is not needed. The intended outcome: data-rich windows get more modes, poorly sampled windows fewer (down to r = 0, only the mean level).
- Caveats: the criterion is only meaningful with realistic data uncertainties. GNSS time series have temporally correlated noise, so n overstates the independent information. Use a data covariance or an effective n (e.g. from the noise model), or confirm with cross-validation (leave stations out).

#### Tasks

- [x] 11.1 KL modes: `basis = unicycle.inverse.klBasis(rcv,iC,opts)`.
  - Matérn covariance between centroids (ν = 0.5 exponential by default; 1.5, 2.5, Inf, or general ν via the Bessel form evaluated in log space), `σ_m` default 1.
  - **Area-weighted** eigenproblem `C W φ = λ φ`, `φ' W φ = I` (the discretized integral operator), so modes and eigenvalues are nearly independent of mesh resolution (leading λ within 3 % between 10 km and 5 km meshes).
  - **Mean removed by default:** the covariance is centred on the area-weighted mean, so every mode has zero mean and the mean level m̄ is a separate parameter (no trade-off with a near-constant first mode).
  - Returns `Phi = φ diag(√λ)` (so `m = m̄ + Phi(:,1:r) c`, c ~ N(0, I)), `phi`, `lambda`, cumulative `varFraction`, `L` (largest principal extent of the centroids) and the options used. Default ℓ = L/2. Sign convention: largest-magnitude entry of each mode positive (signs can differ between meshes; irrelevant for the inversion).
  - No caching: the dense eigendecomposition takes < 1 s for 650 patches and a few seconds for 2600, so it is recomputed.
  - Tests (`testKLBasis`, 7): orthonormality, eigen-equation and exact covariance reconstruction; zero-mean modes; default ℓ = L/2; ℓ ≫ L gives a constant (uncentred) or linear-gradient (centred) leading mode; Matérn forms agree (including large ν without overflow); mesh independence; subset of patches.
  - `studies/klModes.m` → `klModes.png`: on the creeping patches of the 20 km mesh (L ≈ 490 km), the leading modes are an along-strike gradient, an along-dip (shallow vs deep) gradient, centre vs ends along strike, and a quadrupole. Their shapes hardly change between ℓ = L/2 and L. The exponential kernel spreads variance slowly: 2 / 3 / 5 modes capture 42 / 50 / 61 % of the prior variance at ℓ = L/2 (30 / 37 / 48 % at L/4, 50 / 58 / 68 % at L). That is fine for a truncated basis, but a smoother kernel (ν = 1.5 or 2.5) would concentrate more variance in the first few modes.
- [ ] 11.2 Data structure (long format, per-station times, σ, offsets per station-component/segment) and the window forward chain ([Task 10.1](#task-10)) evaluated at the observed times.
- [ ] 11.3 Misfit, Jacobian (finite differences on [m̄_k, c_k], exploiting the triangular structure), and a least-squares solver (Gauss–Newton / Levenberg–Marquardt); later MCMC.
- [ ] 11.4 Model selection over r (and ℓ) with BIC/AIC, sequentially by window then jointly.
- [ ] 11.5 Synthetic tests: recover known η fields (linear-model data), then RSF-generated data ([Task 10.3](#task-10)).

#### Open questions

- ~~Class or plain functions?~~ Plain functions (decided).
- ~~One ℓ for all windows?~~ One ℓ, a common basis (decided).
- Default ℓ = L/2 from the mesh, or always user-supplied? Current plan: user value if given, else L/2.
- Noise model for the time series (white σ only, or white + flicker/random-walk) and how to set the effective n in the information criterion.

<sub>[↑ Contents](#contents)</sub>

---

<a id="task-12"></a>

### 12. Experiment: effective viscosity of an RSF cycle, by time window

> **Status:** 📝 Planning · checklist 0/6

Goal: generate synthetic geodetic data from the asinh RSF cycle, then fit them with the linear viscous model (Tasks 9–11) over a range of time windows. Questions: which windows a linear model can represent and to what misfit; how the fitted η_k evolve with time since the event; whether they follow the RSF tangent viscosity `η_tan = (a−b)σ/v`; and how many windows a log-like RSF decay needs. This merges the "Next" item of 9.4 with 10.3, 10.3b and 11.5 into one experiment.

#### Setup (reference case)

- Same set-up as `studies/linvisCycle.m` / `surfaceTimeseries.m`: `generic_sz15_dx20`, 200 × 100 km rectangular asperity, rake 68°, Vpl = 54 mm/yr, Trec = 150 yr; asinh, a−b = 0.005, σ = 50 MPa, v0 = Vpl, f0 = 0.1.
- RSF truth: spin up with `runCycles` to the periodic cycle. The **modern event is the next characteristic event** at the end of the spun-up cycle, so the RSF cycle supplies both the late-interseismic state before it (t = Trec − τ) and the postseismic response after it (t = 0⁺ … Trec). The model is then exactly consistent with the W0 assumptions (characteristic cycle, tSince = Trec) and only the rheology differs.
- Stations: a land-side grid (north > 0) plus the trench-normal profile of `surfaceTimeseries.m`, through `computeSurfaceKernels`. Use E, N, U displacement relative to long-term motion (`removeLongTerm`, `pinnedSlip`), including the coseismic step.
- Sampling: daily for the first year, weekly up to 10 yr, then monthly to Trec. Dense output from `runCycles` (`nOut`, `tFirst` ≈ 1 day); the closed-form linear model is evaluated at the same times.
- Noise: noise-free first; then white noise of 1 mm horizontal and 3 mm vertical (and later flicker noise, see [Task 11](#task-11) open questions).

#### Data products (saved once, `studies/output/rsfSynthetics.mat`)

- Station displacements and times; for diagnostics only (not used in the fit), the fault slip rate and slip on creeping patches at the same times, and the patch-wise tangent `η_tan(t) = (a−b)σ/|v|` (asinh: `∂τ/∂v` at the current state, which tends to (a−b)σ/v for v ≫ v0 e^{−f0/(a−b)}).
- Δτ_char, slip_P, slip_C of the event, Trec, and the late-interseismic velocities (from the last few years before the event).

#### Windows

- W0: interseismic velocity field over the last 5–10 yr before the event, fitted with `linvisSpinUp(lin_0, periodicEvent, tSince)`.
- Postseismic, log-spaced edges in time since the event, reference set: 0–0.1, 0.1–1, 1–3, 3–10, 10–30, 30–150 yr. Also try coarser sets (2–3 windows) and finer ones (two per decade), and a few contiguous windows of increasing length starting at 0 (0–1, 0–10, 0–150 yr) to show how a single η fails over long windows.

#### Fitting stages

1. **Uniform η, oracle start (isolates the rheology).** For each window, fix the slip rate at the window start to the RSF fault slip rate (not observable, a diagnostic) and line-search one uniform log η over the window's data. Output: misfit vs η curves per window, best η_k, and residual misfit. The first window uses the RSF rate just after the event; the linear jump Δτ/η is not used here.
2. **Uniform η, chained (realistic).** The [Task 10](#task-10) chain: fit η_0 to W0 alone; then the coseismic jump `v(0⁺) = v_W0 + Δτ/η_1`; then η_1, η_2, … sequentially (each given the earlier ones), then jointly. Compare with stage 1: what the carried slip rate and the linear coseismic jump cost.
3. **Spatial structure.** Add r = 1–3 KL modes per window (`klBasis`, ℓ = L/2), select r by BIC per window ([Task 11.4](#task-11)), and compare the fitted log η_k fields with the area-averaged log η_tan over each window.
4. **Window design.** Vary the number and placement of the edges; report the total misfit and BIC against the number of windows. Find the smallest window set that reaches the noise level.
5. **Robustness.** Noise on/off; dense vs sparse stations (e.g. 10 stations); f0 = 0.1 vs 0.6 (from `surfaceTimeseries.m`); a modern event that differs from the characteristic one (rupture extending onto creeping patches, or part of a sequence via `linvisSpinUp(lin,seq)`).

#### Metrics

- Misfit per window (RMS and χ²/n) and the fraction of signal explained.
- η_k against η_tan averaged over the window (area-weighted, over patches with significant afterslip), and against `η_k · v̄_k ≈ (a−b)σ`.
- Prediction (not fitted): fault afterslip per window vs the RSF afterslip, and the late-interseismic coupling pattern from η_0 vs RSF.
- Expected: for log-like RSF decay, v ~ 1/t, so η_tan ∝ t and fitted η_k should grow roughly linearly with the window's mid-time; the first window will need very small η (RSF jumps by 1e5–1e8 Vpl vs Δτ/η ≈ 30 Vpl at η_tan(Vpl), [Task 9.4](#task-9)).
- Known mismatch: RSF has a slow rake-rotation mode that the isotropic linear model lacks ([Task 9.4](#task-9)); check how much misfit it leaves late in the cycle.

#### Deliverables

- `studies/rsfSynthetics.m` (generate and save the data), `studies/rsfToLinvis.m` (stages 1–4), figures in `studies/output/`: η_k vs time since event with η_tan overlaid; misfit vs η per window; data fit at representative stations; maps of fitted log η_k vs log η_tan.

#### Prerequisites

- [Task 10.1](#task-10)–10.2 forward chain over windows (minimum: uniform η per window; then the KL form of [Task 11.2](#task-11)).
- [Task 11.3](#task-11) misfit and solver (a 1-D line search is enough for stages 1–2).

#### Tasks

- [ ] 12.1 RSF synthetic data set (`studies/rsfSynthetics.m`) with the diagnostics above.
- [ ] 12.2 Window forward chain ([Task 10.1](#task-10)) and its tests (10.2).
- [ ] 12.3 Stage 1: uniform η per window with an oracle start; η_k vs η_tan.
- [ ] 12.4 Stage 2: chained uniform fit, W0 first.
- [ ] 12.5 Stage 3: KL modes per window with BIC.
- [ ] 12.6 Stages 4–5: window design and robustness.

#### Open questions

- Modern event = next characteristic event (proposed), or a different event on top of the spun-up cycle?
- Is the oracle start (stage 1) worth keeping in the final figures, or only as a diagnostic?
- Should W0 data be velocities over the last few years (proposed), or displacement time series like the postseismic windows?
- Station layout: generic grid, or a layout modelled on a real network?

<sub>[↑ Contents](#contents)</sub>

---

<a id="task-13"></a>

### 13. Experiment: stochastic event sequence without a locked zone

> **Status:** ✅ Done · checklist 3/3

Question: on a fully velocity-strengthening fault (asinh RSF, no pinned patches), impose events whose slip varies randomly from cycle to cycle. How does the system spin up, and how variable is the late-interseismic coupling? Follows the no-locked-zone runs (periodic 8.1 m events on the 200 × 100 km rectangle, where locking emerges from the stress drop alone).

#### Code

- `unicycle.ode.spinUp` takes a single event or a sequence (`T`, `events`), like `linvisSpinUp` (decided: one entry point, no separate sequence driver). With `opts.repeat = false` it runs the sequence once from `opts.y0` without a convergence check, returning the state before/after every event (`zetaPre`, `thetaPre`, `zetaPost`, `thetaPost`, `slipPre`), the slip of the pass (`slipSince`) and optional dense output per interval (`opts.dense`). `runCycles` accepts sequences too and stores the periodic supercycle.
- `studies/stochasticSequence.m` (one case: `seed`, `nCycles`, `f0`) and `studies/stochasticSequenceCompare.m` (figures).

#### Set-up

- `generic_sz15_dx20`, asinh, a−b = 0.005, σ = 50 MPa, **f0 = 0.3**, v0 = Vpl, Vpl = 54 mm/yr at rake 68°, `pinnedPosition = []`.
- Every 150 yr cycle: three events, one per segment (west, centre, east) at **0, 20 and 30 yr** into the cycle.
- Each event: a square of **uniform slip** along the plate rake. Nominal centres at along-strike −150, 0, +150 km and down-dip distance 90 km, **jittered** per event (σ = 25 km along strike, 15 km down dip). Edge length N(100, 25) km truncated to [40, 160] km. **Slip N(8.1, 1.5) m** (mean = Vpl·150 yr, truncated to > 0), so the events alone close the slip budget on average at the nominal centres (decided; a first run used U(4, 8) m, see below). Patches by centroid in the square (`rectangularAsperity`); edges snap to the 20 km grid.
- Patches near the jittered edges and in the gaps are not ruptured every cycle, so creep still supplies part of their budget.
- Reference: the deterministic sequence of mean events (nominal centres, 100 km, 8.1 m), `seed = 0`.
- 300 cycles (45 kyr) per run, seeds 1–3, start from v = Vpl at the plate rake. About 0.8 s per cycle alone, ~2 s with four runs in parallel.

**Diagnostics** (coupling C = 1 − |v|/Vpl at the end of each cycle, just before the next west event)
- Maps: reference C; stochastic mean and std of C after spin-up; mean − reference; seismic fraction of the slip budget.
- Spin-up: segment-averaged C vs cycle (seeds vs reference).
- Variability: stationary distribution of segment-averaged C.
- Memory: correlation of C with the coseismic slip a patch received k cycles earlier (k = −1, the next cycle, as a null check) and with the slip deficit Vpl t − s.

#### Tasks

- [x] 13.1 `spinUp`/`runCycles` accept sequences; one-pass mode; tests in `testRunCycles`.
- [x] 13.2 `studies/stochasticSequence.m`, `studies/stochasticSequenceCompare.m`.
- [x] 13.3 Run reference + 3 seeds, figures (`studies/output/stochasticSequenceMaps_f0_0.3.png`, `stochasticSequenceStats_f0_0.3.png`).
  - Moving the event loop into `spinUp` reproduces the earlier `runSequence` runs bit for bit (coupling and slip identical for all four cases).

#### Results, slip N(8.1, 1.5) m (f0 = 0.3, 300 cycles, statistics over cycles 51–300 and 3 seeds)

- **The mean-event reference locks its squares completely, but never becomes periodic.** C ≈ 1.00 on all three segments (v down to 1e-10–1e-14 Vpl at the centres), with seismic fraction 1. C stops changing within ~20 cycles (|ΔC| = 3e-5 between the last two), but only because it saturates at 1. The cumulative deficit Vpl t − s at the centre patch keeps falling: −3.4, −7.4, −10.6, −13.9, −16.0 m after 1, 5, 20, 100, 300 cycles. This is the knife edge: the events supply exactly Vpl·T, so any creep is a surplus nothing pays back, and the patch locks ever harder (drift roughly logarithmic in time).
- **The stochastic sequence is less locked on average**, mostly at the edges. Mean segment-averaged C: west 0.79, centre 0.87, east 0.84 (reference ≈ 1.0). Mean C − C_ref reaches −0.4 along the square edges and is ≈ −0.05 in the interiors. Jitter blurs the edges: patches there are inside the rupture in some cycles and creep in others. Segment-averaged seismic fraction 0.68 (reference 1.0); creep supplies the rest at the edges and in the gaps.
- **Spin-up of the coupling is fast:** segment-averaged C reaches its stationary distribution within ~5–10 cycles from v = Vpl.
- **The slip budget does not run away.** At the segment centres the cumulative slip settles 7–25 m ahead of Vpl·t, i.e. those patches sit at a lower stress than the uniform-creep start. This offset wanders by ±5–10 m between 1–2 kyr windows but shows no trend over 45 kyr in any seed (centre patch, seed 1: −15.6 ± 3.8 m over cycles 51–300). So the seismic fraction = 1 knife edge of the reference is avoided: the jittered squares miss the centre patches in some cycles, which gives a small creep share and with it a restoring drift. D settles where the average creep, whose rate depends on D through the traction, fills that shortfall.
- **Why the deficit is not zero.** Budget closure means D/t → 0, not D → 0. D is a traction offset: with backslip loading, τ ≈ τ_start − K D, and the start (D = 0) is steady state at Vpl. Locking a patch late in the cycle needs (a−b)σ ln(v/Vpl) ≈ −3 to −6 MPa for v/Vpl = 1e-5–1e-10, i.e. D ≈ −10 to −20 m with K ≈ μ/L ≈ 0.3 MPa/m for a 100 km square. Even a strictly periodic cycle has a constant nonzero D at a given phase; D → 0 only on a patch back at v = Vpl (C = 0).
- **Variability of late-interseismic coupling.** Segment-averaged 5–95 % ranges: west [0.56, 0.96], centre [0.71, 0.97], east [0.65, 0.97], skewed with a tail to low C. The std map is ≈ 0.05 in the square interiors (locked in almost every cycle) and 0.2–0.25 along the jittered edges, especially the outer along-strike ends of the west and east segments.
- **Memory.** C correlates with the current cycle's coseismic slip (0.75), weakly with the previous cycle's (0.15), and with nothing after ~4 cycles. It does not correlate with the next event (0, as it must by construction). Its correlation with the slip deficit is weaker than with U(4, 8) slip (median −0.56), because C saturates at 1 on over-slipped patches.
- Implication: when events close the budget on average, the interiors of the rupture squares look fully locked in almost every cycle, so the cycle-to-cycle variability of a coupling map sits at the rupture edges. A single snapshot can underestimate the locked area by a jitter length (~25 km) on any side.

<details>
<summary><b>First run, slip U(4, 8) m (mean 6 m, outputs in <code>studies/output/slipUniform4to8/</code>), f0 = 0.3, cycles 51–300, 3 seeds</b> (superseded, click to expand)</summary>

- **Spin-up is fast and unremarkable.** Segment-averaged coupling reaches its stationary range within a few cycles from v = Vpl. The mean-event reference becomes periodic just as quickly (|ΔC| = 4e-9 between its last two cycles). After that, randomness, not spin-up, dominates.
- **The slip budget closes without a locked zone.** The slip deficit Vpl t − s at the segment centres, in the gaps, up-dip and down-dip fluctuates by a few metres but does not drift over 45 kyr. Creep absorbs the shortfall: the seismic fraction is ≈ 0.75 at the segment centres and ≈ 0.5 averaged over the nominal 100 km boxes (0.74 for the reference), because the jittered squares spread the coseismic slip.
- **The mean coupling matches the deterministic one.** The stochastic mean of C is within ±0.06 of the mean-event sequence everywhere and within 0.03 on segment averages: west 0.39 vs 0.36, centre 0.58 vs 0.56, east 0.54 vs 0.54. The west segment is less coupled because its event comes first in the cycle, so it has relaxed for 150 yr at the snapshot against 130/120 yr for centre and east.
- **Late-interseismic coupling varies a lot from cycle to cycle.** Segment-averaged C has a std of 0.12–0.15 and a 5–95 % range of about [0.2, 0.6] (west) and [0.3–0.4, 0.8] (centre, east), i.e. ±50 % around the mean. The std map peaks (≈ 0.17) at the jittered segment edges and in the gaps between segments, where a patch is inside the rupture in some cycles and not in others.
- **Coupling records the past, not the future.** The correlation of C with the coseismic slip a patch received is 0.55 for the current cycle's event, 0.15 one cycle earlier and < 0.05 after three cycles. It is 0 for the next event (a null check, since events are independent). C is strongly tied to the accumulated slip budget: over-slipped patches (s > Vpl t) are under-stressed and stay more locked. The median correlation of C with Vpl t − s is −0.80 over 237 ruptured patches.
- Implication for geodesy: under random event slip, one late-interseismic coupling map reflects the slip-budget state left by the last one to two events. A single snapshot of a single segment can differ from the long-term mean pattern by ±0.2–0.3 in C.

</details>

<sub>[↑ Contents](#contents)</sub>
