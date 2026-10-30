# I2EC_3d — imposed earthquake cycles on a 3D megathrust

Goal: a MATLAB analogue of [seqeas](https://github.com/tobiscode/seqeas-public) built on
`unicycle`. Earthquakes are imposed periodically on a chosen part of a triangular fault
mesh, and the rest of the mesh evolves as an initial value problem governed by
steady-state rate-dependent friction.

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
    - The first f0 sweep (`studies/f0ConvergenceStudy.m`) used v0 = 1e-6 m/s with Vpl ≈ 1.7e-9 m/s, so f(Vpl) = f0 − 0.032. Its f0 = 0.01 and 0.02 runs have negative friction at Vpl under log (v*/Vpl ≈ 80 and 11), i.e. they are fully in the viscous regime. The sweep was rerun with v0 = Vpl and f0 ≥ 0.1 (see Task 4).
- **Decisions (2026-10-30).**
  - **asinh only from here on.** The log law stays in the code (and its tests) for reference, but all experiments use `'asinh'`: the log law cannot represent a reduced or reversed traction (|τ| < 0 below v*), which imposed events next to or inside creeping regions produce.
  - **f0 ≥ 0.1 only (with v0 = Vpl).** We do not think f0 < 0.1 is stable or mechanically meaningful enough to explore through experiments:
    - f(Vpl)/(a−b) = O(1–10) puts the fault in or near the linear-viscous regime of asinh, so it is no longer rate-and-state friction (see the ⚠ note above).
    - |τ| ≈ σ f0 is then only a few MPa, comparable to coseismic stress changes (Δτ ~ 3–12 MPa near the rupture). Events then reverse tractions and rotate slip by up to 180° (Task 5), and the (ζ, θ) state approaches its singularity.
    - The system becomes stiff for `ode45` (the f0 = 0.01 sweep with v0 = 1e-6 did not finish one run in ~2 h).
  - f0 = 0.1 is the lower end of the explored range and is checked in every new experiment.
- **State per creeping patch:** `[slip_s, slip_d, ζ = ln(|v|/v0), θ = rake]`.
  - This state is singular at v = 0: ζ → −∞, and `θ̇ = τ̇⊥/(|τ| + η|v|)` becomes very large as |τ| → 0. That is stiff but finite unless the traction passes exactly through zero.
  - Option, if strong stress shadows or slip reversal matter (especially in the viscous regime): use the traction components (τs, τd) as the state, with v from the asinh law. That has no singularity.
- **Elasticity:** `τ̇ = K (v − Vpl)` with **v = 0 on pinned patches**. On creeping patches this is `τ̇_c = K_cc (v_c − Vpl_c) + K_cp (0 − Vpl_p) = K_cc v_c − K(c,:)·Vpl`, so the pinned columns enter only through the constant `τ̇_load`. Checked against the full `K (v − Vpl)` with `v_p = 0`: relative error 4e-15.
- **Projecting τ̇ onto the slip direction and its normal:**
  - Parallel: `τ̇·v̂ = (d|τ|/dζ + η|v|) ζ̇`, with `d|τ|/dζ = (a−b)σ` (exact for log, and to round-off for asinh). This sets acceleration and deceleration; f0 and v0 drop out, exactly as in 2D.
  - Perpendicular: `τ̇·v̂⊥ = (|τ| + η|v|) θ̇`. This sets rotation; because `|τ|/(a−b)σ = f0/(a−b) + ζ`, the static friction level matters.

## Code layout

| function | role |
|---|---|
| `unicycle.greens.computeKernels(rcv,name,kdir,nCheck)` | load or compute cached traction kernels |
| `unicycle.ode.assembleODE(rcv,bigK,law)` | one-time setup: partitioned kernel (`Kcc`, `Kcp`), loading, `y0` |
| `unicycle.ode.steadyState(sys)` | interseismic steady state `v = K_cc \ (−τ̇_load)` (asperity locked forever) as a state vector; reference for single-event experiments, not the cycle-spin-up start (see Task 6) |
| `unicycle.ode.eventStressChange(sys,Trec,weight,weightC)` | coseismic slip `weight·Vpl·Trec` on pinned patches (and `weightC·Vpl·Trec` on creeping patches inside the rupture); `Δτ = K_cp slip_P + K_cc slip_C` |
| `unicycle.ode.ruptureSlipForStressChange(sys,Trec,inRuptureC,dtauTarget,direction)` | slip on creeping patches inside the rupture that gives each of them a prescribed Δτ·d (e.g. a 3 MPa stress drop) |
| `unicycle.ode.periodicEvent(sys,Trec,weightC)` | event that releases exactly Vpl·Trec on the locked patches (closes their slip budget), plus optional slip on creeping rupture patches |
| `unicycle.ode.runCycles(sys,event,opts)` | spin-up: event + `ode45` over Trec per cycle until periodic (slip per cycle = Δs*, ζ and rake repeat); stores one full cycle |
| `unicycle.ode.stateTraction(sys,y)` | traction vector implied by the state, `|T| = F(v) + ηv` along θ |
| `unicycle.ode.applyEvent(sys,y,dtau,slipC)` | coseismic state jump `T⁺ = T⁻ + Δτ`, damped inverse for v⁺, θ⁺ from T⁺; adds `slipC` to the slip state |
| `unicycle.ode.odeRateDependent(t,y,sys)` | ODE right-hand side |
| `unicycle.ode.frictionLaw(sys,zeta)` | `F = |τ|` and `dF/dζ` for `'log'` or `'asinh'` |
| `unicycle.ode.frictionInverse(sys,tauMag,useDamping)` | `ζ` from the traction, with radiation damping by default |

`+unicycle/Contents.m` lists an older `unicycle/ode/*` (evolution, rateandstate, …) that is not in this copy, so there are no name clashes.

**Tests:** `tests/`, run with `matlab -batch "run('tests/runAllTests.m')"`. There are 52 tests and all pass; the run takes about 1 min.

| file | covers |
|---|---|
| `testGreensConventions` | `.ned` → `rcv.x` axis swap; unit vectors; positive dip slip = thrust; positive strike slip = left-lateral; negative self-stiffness |
| `testComputeKernels` | block ordering; reload; changed G, nu or geometry → new file; tampered or foreign cache caught |
| `testFrictionLaw` | log law; asinh = log for realistic f0/(a−b); asinh regularization at low v; derivative; inverse round trip; no overflow; damped inverse |
| `testAssembleODE` | partitioning; scalar expansion; Vpl/Vrake components; loading `K(v−Vpl)` with `v_p = 0`; steady state; `y0`; validation and warning |
| `testImposedEvent` | `Kcp` partition; `Δτ` matches the full kernel and is linear in the weight; nearest neighbours loaded in the slip direction; `applyEvent` ∘ `stateTraction` is the identity for Δτ = 0; traction balance and direction after rotating/amplifying/reducing steps; damping caps v⁺; spring-slider afterslip → Δτ/k; `steadyState` is a fixed point; log-law negative traction raises an error; rupture onto creeping patches (Δτ from the full kernel, stress drop, slip jump); traction reversal with asinh gives backward slip ∥ T⁺; slip solved for a prescribed stress change hits the target exactly |
| `testRunCycles` | locked-budget closure gives Δs* = Vpl·Trec on creeping patches (with and without creeping rupture slip); spring-slider cycles converge to Vpl·T + Δτ/k; small-mesh spin-up at f0 = 0.1 closes the slip budget and repeats; warning when the locked budget does not close |
| `testOdeRateDependent` | projection formulas; f0-independent acceleration; `θ̇ ∝ 1/|τ|`; spring-slider logistic solution; rake relaxation `exp(−k Vpl t/|τ|)`; radiation damping; damped traction balance conserved; small mesh converges to BEM |

Helpers: `makeTestMesh` (60 × 30 km plane, 36 triangles), `setTestProperties` and `oneSliderSystem`.

## Tasks

### 1. Mesh — done
- [x] `make_fault_mesh.m`: planar 500 × 300 km megathrust dipping 15°, 3000 triangles of about 10 km.
  - Writes `fault/generic_sz15.{ned,tri}`.

### 2. Kernels
- [x] `unicycle.greens.computeKernels(rcv,name,kdir,nCheck)` returns `bigK` (2N × 2N shear) and `bigKn = [Ksn,Kdn]`. It caches them in `kernels/<name>_<hash>.mat`.
  - The hash is a SHA-256 of `rcv.x`, `rcv.vertices`, the Green's function class, G and nu. A different geometry or elastic model gets a different file, so stale kernels can't be picked up.
  - On load it compares the full stored fingerprint, then recomputes `nCheck` (default 2) random source columns and compares them with the cache. That takes about 1 s for 3000 patches, against about 260 s for the full computation.
  - Tested: reload, changed G, moved vertex, a tampered kernel entry, and a cache file renamed from another mesh are all handled or caught.
  - The file is 415 MB because it also stores `Ksn`/`Kdn` (6 blocks instead of 4).
- [x] Cached for `generic_sz15`: `kernels/generic_sz15_0f17a6d68309.mat`.

### 3. Interseismic ODE (rate-dependent friction) — done
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

### 4. Interseismic verification
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
    - ~~Start cycle spin-up from the BEM steady state.~~ Superseded by Task 6: with periodic events the cycle-averaged slip is Vpl along the plate rake, so v = Vpl at the plate rake is the better start (2.5× fewer cycles at f0 = 0.6).
    - The periodicity check needs a rake criterion, not just magnitude.
- [x] Compare log and asinh on the same sweep: plate rake 90°, `plateRake = 90; law = 'log'; run('studies/f0ConvergenceStudy.m')` (outputs tagged `_rake90_log`; asinh keeps the untagged-law names).
  - They agree to round-off. At 1e4 yr the final states differ by at most 3e-14·Vpl in |v| and 3e-13° in rake, and the error-vs-time curves differ by at most 1e-13°. The ratio |τ|/(dF/dζ) at Vpl is the same (20, 40, 80, 120), and dF/dζ = (a−b)σ exactly.
  - So every rake-90° number above holds for the log law too: max rake error 0.2°, 3.5°, 14.5°, 22.9° at 1e4 yr.
  - This is expected. The laws differ only for v ≲ v0·exp(−f0/(a−b)) = Vpl·e^−20 or smaller, and the slowest creeping patch runs at 0.06·Vpl.
- [x] Integrator: **`ode45`**, `RelTol = AbsTol = 1e-8` (decided 2026-10-30). One 150-yr postseismic cycle on `generic_sz15_dx20` takes ~0.3 s with v0 = Vpl. It becomes very slow only when the asinh law is in its viscous regime (f(Vpl)/(a−b) ~ O(1)); revisit (`ode15s` with a Jacobian, or the traction-state formulation) only if spin-up runs become slow.

### 5. Imposed earthquakes — done
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
  - Coseismic rotation of the slip direction is only ≤ 6.5°: the rotation is ≈ Δτ⊥/|T| and |T| ≈ σ f0 = 30 MPa ≫ Δτ. So with f0 = 0.6, f0 matters for rake much more through the slow interseismic relaxation (Task 4) than through the coseismic jump. Smaller f0 would give larger coseismic rotations.
  - After one cycle the state is far from where it started (max |v/v_ss − 1| = 0.68, max Δrake 5°): the steady state is not the periodic cycle, so spin-up (Task 6) is needed.

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

### 6. Cycle driver and spin-up — done
- [x] Loop over cycles (`unicycle.ode.runCycles`): apply the event to the pre-event state, integrate Trec with `ode45`, repeat. The slip state is reset each cycle, so its end value is the slip of that cycle.
- [x] **Slip budget.** `unicycle.ode.periodicEvent` makes the locked patches release exactly Vpl·Trec along the plate rake in every event (`runCycles` warns `unicycle:ode:runCycles:lockedBudget` otherwise). Over one cycle the creeping traction changes by `ΔT = Δτ_event + K_cc Δs_inter + τ̇_load Trec`, so the periodic cycle has a unique slip per cycle `Δs* = slip_C − K_cc \ (Δτ_event + τ̇_load Trec)`, which equals Vpl·Trec on every creeping patch when the locked budget closes (also with coseismic slip on creeping patches). The slip budget therefore closes everywhere once the cycle is periodic.
- [x] Periodicity check, on successive pre-event states: `max |Δs − Δs*|/(Vpl Trec) < 1e-3`, `max |Δζ| < 1e-3`, `max |Δrake| < 1e-3°`, for 3 cycles in a row. Then one full cycle is stored (161 log-spaced outputs from 1e-8 yr to Trec).
- [x] Study: `studies/cycleSpinup.m` (parameters `f0`, `ruptureExtend`, `dtauTarget`, `maxCycles`, `startFrom`) → `studies/output/cycleSpinup_f0_<f0>_ext<km>km[_drop3].{png,mat}`. Same setup as Task 5: `generic_sz15_dx20`, 200 × 100 km asperity, rake 68°, a−b = 0.005, v0 = Vpl, Trec = 150 yr (8.1 m per event), asinh.

  | rupture | f0 | cycles to periodic (from plate rate) | from BEM steady state | wall time | slip/(Vpl T) in stored cycle | min \|T⁺\| over spin-up | reversed |
  |---|---|---|---|---|---|---|---|
  | asperity | 0.6 | 78 (11.7 kyr) | 194 (29 kyr) | 17 s | 0.99996–1.00005 | 29.9 MPa | 0 |
  | asperity | 0.1 | 31 (4.7 kyr) | 43 (6.5 kyr) | 6 s | 1.00000 | 4.94 MPa | 0 |
  | +20 %, 3 MPa drop | 0.6 | 77 (11.6 kyr) | 206 (31 kyr) | 18 s | 0.99996–1.00004 | 26.6 MPa | 0 |
  | +20 %, 3 MPa drop | 0.1 | 31 (4.7 kyr) | 44 (6.6 kyr) | 7 s | 1.00000 | 1.59 MPa | 0 |

  - **The slip budget closes** to < 5e-5 of Vpl·Trec on every creeping patch, and the converged state repeats (Δζ ~ 1e-6, Δrake ~ 1e-3° per cycle).
  - **Convergence is geometric** and set by the slow rake mode: ~4.5× fewer cycles at f0 = 0.1 than at 0.6, consistent with the rake time scale ∝ |τ| ∝ f0 (Task 4).
  - **The periodic cycle is not the BEM steady state.** With the asperity released every cycle, every creeping patch slips Vpl·Trec along (close to) the plate rake, whereas the steady state around a permanently locked asperity has rakes rotated by up to −58°/+34° and speeds 0.85–1.8 × v_ss. In the periodic cycle the rake stays within ±0.5° of 68° (f0 = 0.6) or within a few degrees except right after the event (f0 = 0.1, e.g. 84° → 68° at the west edge of the extended rupture). Hence the spin-up default is now `sys.y0` (v = Vpl at the plate rake), 1.4–2.5× fewer cycles than starting from the steady state; both reach the same periodic cycle.
  - **f0 = 0.1 is stable through the whole spin-up**, including the extended rupture with a 3 MPa drop: no reversals, min |T⁺| = 1.59 MPa (above the 1 MPa threshold), each cycle ~0.2 s.
  - Late-interseismic coupling (1 − |v|/Vpl just before the event) is high on the ring of creeping patches around the asperity and decays over ~100 km.
- [ ] Save `t`, slip and velocity on **all** patches (pinned patches: v = 0 between events, coseismic slip at events) — move to Task 7 together with surface displacements.

### 7. Outputs
- [ ] Surface displacements and velocities at observation points via `rcv.displacementKernels` (relative to Vpl, as in `BEM_interseismicsolve.m`).
- [ ] Plots: map views of v/Vpl and rake at selected times, time series at selected patches, and the cumulative slip budget (total slip per cycle ≈ Vpl·T_rec everywhere).

### 8. Later / optional
- [ ] Multiple asperities with different recurrence times and non-uniform coseismic slip.
- [ ] Alternative laws behind the same interface: the regularized sinh form and power-law viscous `v = τⁿ/η`.
- [ ] Apply the model to the Mexico mesh (`fault/mexico_mt_xyz`) with Euler-pole loading. This needs the outline and lat/lon files that are not in this folder.
- [ ] Performance for larger meshes: compressed kernels such as H-matrices, or a GPU.
- [ ] Benchmark against seqeas in an equivalent configuration.
