# I2EC_3d: imposed earthquake cycles on a 3D fault

MATLAB code for earthquake-cycle models in which earthquakes are **imposed** periodically on a locked part of a 3D triangular fault mesh, and the rest of the fault **creeps** according to rate-strengthening friction. It follows the approach of [seqeas](https://github.com/tobiscode/seqeas-public), extended to 3D with two-component (strike and dip) slip. It is built on the [unicycle](https://bitbucket.org/sbarbot/unicycle/src/master/) library, included under `+unicycle/`.

The model:
- **Elasticity:** traction kernels between triangular patches in a half-space (Nikkhoo & Walter 2015), `τ̇ = K (v − Vpl)`, with `v = 0` on locked patches between events.
- **Friction:** steady-state rate-strengthening friction on creeping patches, `|τ| = (a−b)σ asinh(|v|/(2v0) · exp(f0/(a−b)))`, with radiation damping `G/(2Vs)`. The slip-rate vector is assumed collinear with the traction vector; the friction law sets their magnitudes.
- **Events:** each event releases the slip deficit of the locked region, `Vpl·Trec` along the plate rake, so its slip budget closes. Optionally the rupture extends onto creeping patches with a prescribed stress drop.
- **Spin-up:** the characteristic event is repeated until the cycle is periodic. The resulting state is the primed initial condition on which recorded (instrumental-era) events and creep are imposed.

Status, decisions and results of all experiments so far are in [tasks.md](tasks.md).

## Requirements

- MATLAB (developed and tested with R2025a). No toolboxes beyond core MATLAB.
- Run everything from the repository root, so that the `+unicycle` package is on the path.

## Quick start

```matlab
yr = 3.1536e7;
addpath('studies')                                    % rectangularAsperity

% mesh and kernels (cached in kernels/ after the first run)
eM   = unicycle.greens.nikkhoo15(30e3,0.25);          % shear modulus (MPa), Poisson's ratio
rcv  = unicycle.geometry.triangleReceiver('fault/generic_sz15_dx20',eM);
bigK = unicycle.greens.computeKernels(rcv,'generic_sz15_dx20');

% properties, written by the user into the mesh object (m, s, MPa)
rcv.a = 0.015;  rcv.b = 0.010;  rcv.sigma = 50;  rcv.mu0 = 0.6;
rcv.Vpl = 0.054/yr;  rcv.Vo = rcv.Vpl;  rcv.Vrake = 68;  rcv.Vs = 3e3;
rcv.pinnedPosition = rectangularAsperity(rcv,[-110e3,90e3],[40e3,140e3],15);

% periodic characteristic event every 150 yr
sys   = unicycle.ode.assembleODE(rcv,bigK,'asinh');
event = unicycle.ode.periodicEvent(sys,150*yr);

% spun-up initial condition, 120 yr after the last characteristic event
[y0,info] = unicycle.ode.spinUp(sys,event,struct('tSince',120*yr));

% or: spin up and store one full periodic cycle
out = unicycle.ode.runCycles(sys,event);
```

On the 750-triangle mesh, one cycle takes about 0.2 s. Spin-up takes about 30 cycles at f0 = 0.1 and about 80 at f0 = 0.6.

## Workflow

| step | call |
|---|---|
| Mesh | `rcv = unicycle.geometry.triangleReceiver('fault/<name>', eM)` reads `<name>.ned` / `<name>.tri` |
| Kernels | `bigK = unicycle.greens.computeKernels(rcv, '<name>')` loads cached kernels if they match the mesh geometry and elastic constants, otherwise computes and caches them |
| Properties | set `a, b, sigma, mu0, Vo, Vpl, Vrake, Vs, pinnedPosition` on `rcv` (scalar or one value per patch) |
| System | `sys = unicycle.ode.assembleODE(rcv, bigK, 'asinh')` (rebuild after changing any property) |
| Event | `event = unicycle.ode.periodicEvent(sys, Trec, weightC)` |
| Spin-up | `[y0, info] = unicycle.ode.spinUp(sys, event, opts)` or `out = unicycle.ode.runCycles(sys, event, opts)` |
| Single events | `unicycle.ode.applyEvent(sys, y, dtau, slipC)`, then `ode45(@(t,y) unicycle.ode.odeRateDependent(t,y,sys), ...)` |
| Linear viscous, analytic | `lin = unicycle.ode.linvisAssemble(rcv, bigK, eta)`, `[v, info] = unicycle.ode.linvisSpinUp(lin, unicycle.ode.periodicEvent(lin, Trec), tSince)`, `[v, info] = unicycle.ode.linvisSpinUp(lin, seq, tOut)`, `[v, s] = unicycle.ode.linvisEvolve(lin, v0, t)` |
| Surface stations | `G = unicycle.greens.computeSurfaceKernels(rcv, xObs)` once, then `[vSurf, uSurf] = unicycle.ode.surfaceDeformation(G, sys, t, Y, opts)` |

The ODE state covers the creeping patches only (`sys.Nc` of them, mesh indices `sys.iC`; locked patches are `sys.iP`):

```
y = [slip_strike; slip_dip; zeta; theta]     zeta = log(|v|/v0), theta = rake (rad)
```

Slip rate is `v = sys.v0 .* exp(zeta)` along `[cos(theta), sin(theta)]` in (strike, dip). In `runCycles` output, `out.t` is the time since the event and `out.y` holds the state of one periodic cycle, with slip relative to the start of the cycle. Locked patches have `v = 0` between events and slip `out.event.slipP` at each event.

### Surface time series

```matlab
xObs = [east, north, zeros(D,1)];                       % stations (m), same frame as rcv.xc
G    = unicycle.greens.computeSurfaceKernels(rcv,xObs);  % (3D x 2N)
o    = struct('removeLongTerm',true,'pinnedSlip',out.event.slipP);
[vSurf,uSurf] = unicycle.ode.surfaceDeformation(G,sys,out.t,out.y,o);
% vSurf (m/s), uSurf (m): nt x D x 3, components [east, north, up]
```

- `removeLongTerm` subtracts the deformation from steady slip at Vpl on the whole mesh, `G·Vpl` (and `G·Vpl·t` from the displacement), leaving the signal relative to the long-term motion.
- `pinnedSlip` adds the coseismic slip of the locked patches, so the displacement includes the coseismic step at t = 0.
- The displacement comes from the slip state, which `ode45` integrates exactly along with the slip rate, rather than from numerical quadrature of the velocity.
- Keep stations off the fault trace, where the kernels are singular.

### Friction parameters to vary

| property | meaning | typical |
|---|---|---|
| `a`, `b` | rate-and-state parameters, `a > b` on creeping patches | 0.015, 0.010 |
| `sigma` | effective normal stress (MPa) | 50 |
| `mu0` | f0, friction coefficient at v0 | 0.1–0.6 |
| `Vo` | v0; set it to `Vpl` so that f0 is the friction at plate rate | `rcv.Vpl` |
| `Vs` | shear-wave speed (m/s), for radiation damping | 3e3 |

## Conventions

- **Units:** m, s and MPa. Velocities are in m/s, even though the unicycle constructor defaults are in m/yr, so always overwrite `Vs` and `Vo`.
- **Coordinates:** `.ned` files are `[index, north, east, depth (positive down)]`. After loading, `rcv.x` / `rcv.xc` are `[east, north, up]`.
- **Slip signs:** positive dip slip is thrust; positive strike slip is left-lateral. Kernels are `bigK = [Kss, Kds; Ksd, Kdd]`.
- **Who sets properties:** the user writes all physical properties into the mesh object after reading it; library code only reads them.

## Linear viscous friction (analytic)

With `τ = τ0 + η v` on creeping patches (isotropic η per patch, MPa·s/m), the cycle ODE is linear, `v̇ = D⁻¹(K_cc v + τ̇_load)`, and is solved in closed form from the eigendecomposition of `A = D⁻¹K_cc`. That includes the exact periodic cycle for a characteristic event, or for any sequence of events that repeats every supercycle T: no time stepping, no spin-up loop. A sequence must release the slip deficit of the locked patches (Σ slip_P = Vpl·T, checked); events may rupture locked and creeping patches. The slip per cycle is Vpl·Trec for any η, so only transients constrain η. This is the forward model for inverting effective viscosity structure; see Task 9 in tasks.md.

```matlab
lin   = unicycle.ode.linvisAssemble(rcv,bigK,eta);      % eta: scalar, Nc x 1 or N x 1
event = unicycle.ode.periodicEvent(lin,Trec);
[vNow,info] = unicycle.ode.linvisSpinUp(lin,event,tSince);   % exact periodic state
[v,s] = unicycle.ode.linvisEvolve(lin,info.vPost,t);          % 2Nc x nt, closed form

% sequence of events (slip on any patches, N x 2), repeating every T
events(1) = unicycle.ode.sequenceEvent(lin,0,slip1);
events(2) = unicycle.ode.sequenceEvent(lin,t2,slip2);
[v,info] = unicycle.ode.linvisSpinUp(lin,struct('T',T,'events',events),tOut);   % same function
```

## Modelling choices

- **asinh friction only** for experiments. The log law is kept in the code for reference, but it cannot represent the reduced or reversed tractions that imposed events produce next to and inside creeping regions.
- **f0 ≥ 0.1 with v0 = Vpl.** Below this the fault is in the near-linear-viscous regime of the asinh law, tractions are comparable to coseismic stress changes, and the system becomes stiff. We do not consider it mechanically meaningful to explore.
- **`ode45`** with `RelTol = AbsTol = 1e-8`.

## Repository layout

```
+unicycle/
  +greens/computeKernels.m       cached traction kernels with geometry check
  +greens/computeSurfaceKernels.m  surface displacement kernels at stations
  +ode/                          imposed-earthquake-cycle machinery
    assembleODE.m                partition locked/creeping patches, loading, friction terms
    odeRateDependent.m           ODE right-hand side
    frictionLaw.m, frictionInverse.m   asinh/log laws, radiation-damped inverse
    steadyState.m                interseismic steady state with the asperity locked forever
    eventStressChange.m          coseismic slip and traction change
    ruptureSlipForStressChange.m slip on creeping rupture patches for a prescribed stress drop
    stateTraction.m, applyEvent.m      traction from state; coseismic state jump
    surfaceDeformation.m         surface velocity and displacement time series
    linvisAssemble.m, linvisEvolve.m, linvisEvent.m, linvisSpinUp.m, linvisExpm1.m
                                 linear viscous friction: closed-form evolution and exact periodic cycle
    sequenceEvent.m              one event of a sequence (linvisSpinUp also takes repeating sequences)
    odeLinearViscous.m, linvisSpinUpNumerical.m
                                 linear viscous ODE and time-stepped spin-up, to verify the closed form
    periodicEvent.m              characteristic event that closes the locked slip budget
    spinUp.m, runCycles.m        spin-up to the periodic cycle
  +inverse/                      inversion machinery (plain functions)
    klBasis.m                    Karhunen-Loève basis of a Matérn covariance on patches (low-rank log-viscosity fields)
  (rest of unicycle: geometry, Green's functions, utilities)
fault/                           meshes (.ned/.tri); make_fault_mesh.m writes the generic ones
kernels/                         kernel cache (not tracked)
studies/                         experiment scripts; outputs in studies/output/ (not tracked)
tests/                           matlab.unittest suite
tasks.md                         plan, decisions and results
```

Meshes:
- **`generic_sz15`:** planar fault, 500 km along strike × 300 km down dip, dip 15°, 10 km triangles (3000).
- **`generic_sz15_dx20`:** the same fault with 20 km triangles (750). Use this one for fast experiments.

## Studies

| script | what it does |
|---|---|
| `studies/f0ConvergenceStudy.m` | interseismic relaxation to the steady state for several f0; rake relaxes on a time scale ∝ f0 |
| `studies/postseismicExample.m` | one imposed event and 150 yr of postseismic response (`f0`, `ruptureExtend`, `dtauTarget`) |
| `studies/postseismicCompare.m` | compares the `postseismicExample` cases |
| `studies/cycleSpinup.m` | spin-up to the periodic cycle, slip budget and convergence figures |
| `studies/surfaceTimeseries.m` | surface velocity and displacement over the periodic cycle along a trench-normal profile |
| `studies/linvisCycle.m` | analytic linear viscous cycle vs the RSF cycle |
| `studies/klModes.m` | leading KL modes on the creeping patches for several correlation lengths |
| `studies/plotLockedAsperity.m` | plots the mesh and locked region |

## Tests

```bash
/Applications/MATLAB_R2025a.app/bin/matlab -batch "run('tests/runAllTests.m')"
```

78 tests cover:
- sign and axis conventions;
- kernel caching;
- the friction laws and the damped inverse;
- loading with locked patches;
- spring-slider solutions of the ODE;
- convergence to the boundary-element steady state;
- imposed events (traction balance, reversal, stress-drop solve);
- cycle spin-up (slip-budget closure, independence of the starting state);
- surface deformation (kernel signs, long-term removal, displacement as the time integral of velocity);
- the analytic linear viscous model (against `ode45`, a numerically spun-up cycle, the iterated cycle map and the slip budget);
- the KL basis for low-rank viscosity fields.

The suite takes about a minute.
