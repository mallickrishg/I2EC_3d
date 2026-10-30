function out = runCycles(sys,event,opts)
% RUNCYCLES spin up periodic imposed earthquake cycles until the cycle is
% periodic, then store one full cycle
%
%   out = unicycle.ode.runCycles(sys,event,opts)
%
% Runs unicycle.ode.spinUp (see there for the cycle loop, the slip budget
% and the convergence criteria) and then integrates one more cycle from
% the periodic pre-event state with output times. Use spinUp directly when
% only the spun-up initial condition is needed.
%
% INPUT:
%   sys   - structure from unicycle.ode.assembleODE
%   event - structure from unicycle.ode.periodicEvent
%   opts  - (optional) structure with fields (defaults in brackets)
%           y0         - initial pre-event state [sys.y0, v = Vpl at the
%                        plate rake]. In the periodic cycle every creeping
%                        patch slips Vpl Trec along (close to) the plate
%                        rake, so this converges faster than the
%                        interseismic steady state, whose rake is rotated
%                        by up to ~60 deg around a permanently locked asperity
%           maxCycles  - [2000]
%           tolSlip    - [1e-3]
%           tolZeta    - [1e-3]
%           tolRake    - [1e-3] (deg)
%           nConverged - [3]
%           odeOpts    - [odeset('RelTol',1e-8,'AbsTol',1e-8)]
%           nOut       - number of log-spaced output times in the stored
%                        cycle [161]
%           tFirst     - first output time after the event (s) [1e-8 yr]
%           verbose    - print progress every this many cycles [50], 0 = off
%           (y0, maxCycles, tol*, nConverged, odeOpts, verbose are passed
%           to unicycle.ode.spinUp)
%
% OUTPUT: structure out with
%   converged, nCycles
%   history  - per-cycle table: slipErr, dzeta, drakeDeg, minTplus,
%              nReversed (traction reversed by the event), runtime
%   yStart   - (4Nc x nCycles+1) pre-event state at the start of each cycle
%   t,y      - one stored cycle starting from the converged pre-event
%              state: t (nt x 1) from the event, y (nt x 4Nc), with slip
%              relative to the start of the cycle
%   slipCycle - (Nc x 2) slip of the stored cycle
%   slipTarget - (Nc x 2) slip per cycle of the periodic solution, Δs*
%   event    - the event structure
%
% SEE ALSO: unicycle.ode.spinUp, unicycle.ode.periodicEvent, unicycle.ode.applyEvent,
%           unicycle.ode.odeRateDependent

if nargin < 3, opts = struct(); end
yr = 3.1536e7;
if ~isfield(opts,'nOut') || isempty(opts.nOut), opts.nOut = 161; end
if ~isfield(opts,'tFirst') || isempty(opts.tFirst), opts.tFirst = 1e-8*yr; end
Nc = sys.Nc;
Trec = event.Trec;

spinOpts = rmfield(opts,intersect(fieldnames(opts),{'nOut','tFirst','tSince'}));
spinOpts.keepStates = true;
[~,info] = unicycle.ode.spinUp(sys,event,spinOpts);
if ~isfield(opts,'odeOpts') || isempty(opts.odeOpts)
    opts.odeOpts = odeset('RelTol',1e-8,'AbsTol',1e-8);
end

% store one full cycle from the periodic pre-event state
tout = [0,logspace(log10(opts.tFirst),log10(Trec),opts.nOut)];
yPlus = unicycle.ode.applyEvent(sys,info.yPre,info.event.dtau,info.event.slipC);
[t,Y] = ode45(@(t,y) unicycle.ode.odeRateDependent(t,y,sys),tout,yPlus,opts.odeOpts);

out.converged = info.converged;
out.nCycles = info.nCycles;
out.history = info.history;
out.contraction = info.contraction;
out.yStart = info.yStart;
out.t = t;
out.y = Y;
out.slipCycle = [Y(end,1:Nc)',Y(end,Nc+1:2*Nc)'];
out.slipTarget = info.slipTarget;
out.event = info.event;

end
