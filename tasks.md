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

- **Friction:** `|τ| = σ (f0 + (a−b) ln(|v|/v0))`, with v collinear with τ. Collinearity is an assumption.
- **State per creeping patch:** `[slip_s, slip_d, ζ = ln(|v|/v0), θ = rake]`.
- **Elasticity:** `τ̇ = K_cc v_c − K(c,:)·Vpl`. This is backslip loading with v = 0 on pinned patches.
- **Projecting τ̇ onto the slip direction and its normal:**
  - Parallel: `τ̇·v̂ = ((a−b)σ + η|v|) ζ̇`. This sets acceleration and deceleration; f0 and v0 drop out, exactly as in 2D.
  - Perpendicular: `τ̇·v̂⊥ = (|τ| + η|v|) θ̇`. This sets rotation; because `|τ|/(a−b)σ = f0/(a−b) + ζ`, the static friction level matters.

## Tasks

### 1. Mesh — done
- [x] `make_fault_mesh.m`: planar 500 × 300 km megathrust dipping 15°, 3000 triangles of about 10 km.
  - Writes `fault/generic_sz15.{ned,tri}`.

### 2. Kernels
- [x] Compute traction kernels for `generic_sz15` and cache them in `kernels/generic_sz15_traction.mat` (about 4 min).
- [ ] `computeKernels.m`: a helper that computes or loads cached kernels for any mesh basename and returns `bigK`.

### 3. Interseismic ODE (rate-dependent friction)
- [x] First draft: `prepareMesh.m` and `odeRateDependent.m`.
- [ ] Refactor `prepareMesh(rcv,bigK)` to read every parameter from the `triangleReceiver` properties above, so callers don't pass a duplicate struct.
  - It returns only the derived quantities: `iC`, `Kcc`, `taudotLoad`, `y0`, and per-patch `(a−b)σ`, `η`.
  - Map `Vpl`/`Vrake` to the strike and dip components `[Vpl·cosd(Vrake), Vpl·sind(Vrake)]`.
- [ ] Update `odeRateDependent(t,y,ode)` to match.

### 4. Interseismic verification
- [x] With no pinned patches, steady sliding at Vpl is a fixed point: |ẏ| ~ 1e-16.
- [ ] Long-time convergence to the BEM steady state `v_c = K_cc \ (−τ̇_load)`.
  - The first run was not converged after 3000 yr: max error was 0.23·Vpl and the rake was off by up to 12°.
  - Hypothesis: the magnitude relaxes quickly, but the rake relaxes about `f0/(a−b)` times slower.
  - To test: separate the magnitude and rake errors, and vary f0.
- [ ] Choose an integrator and tolerances (`ode45` vs `ode23s`/`ode15s`) and record the runtime.

### 5. Imposed earthquakes
- [ ] Define the asperity with simple geometry (an ellipse or box in x and depth) and set it as `rcv.pinnedPosition`.
- [ ] Coseismic slip on the asperity: `Vpl·T_rec` along the plate rake, with an optional edge taper.
- [ ] Coseismic stress change on creeping patches: `Δτ = K_cA · slip_A` (two components).
- [ ] Coseismic state update from the new traction vector `T⁺ = T⁻ + Δτ`:
  - `θ⁺ = atan2(T⁺_d, T⁺_s)`
  - `ζ⁺ = (|T⁺|/σ − f0)/(a−b)`
  - The jump is exact in traction space; this is where f0 enters coseismically too.
- [ ] Cap velocity or apply radiation damping so the postseismic start is well posed.

### 6. Cycle driver and spin-up
- [ ] Loop over cycles: integrate the interseismic period, apply the event, and repeat.
- [ ] Periodicity check: compare the state at the start of successive cycles against `atol`/`rtol`. Stop once it is converged, then store one full cycle.
- [ ] Save `t`, slip and velocity on all patches. Pinned patches get v = 0 between events and their coseismic slip at events.

### 7. Outputs
- [ ] Surface displacements and velocities at observation points via `rcv.displacementKernels` (relative to Vpl, as in `BEM_interseismicsolve.m`).
- [ ] Plots: map views of v/Vpl and rake at selected times, time series at selected patches, and the cumulative slip budget (total slip per cycle ≈ Vpl·T_rec everywhere).

### 8. Later / optional
- [ ] Multiple asperities with different recurrence times and non-uniform coseismic slip.
- [ ] Alternative laws behind the same interface: the regularized sinh form and power-law viscous `v = τⁿ/η`.
- [ ] Apply the model to the Mexico mesh (`fault/mexico_mt_xyz`) with Euler-pole loading. This needs the outline and lat/lon files that are not in this folder.
- [ ] Performance for larger meshes: compressed kernels such as H-matrices, or a GPU.
- [ ] Benchmark against seqeas in an equivalent configuration.
