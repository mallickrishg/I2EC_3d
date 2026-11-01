function out = runCycles(sys,event,opts)
% RUNCYCLES spin up periodic imposed earthquake cycles until the cycle is
% periodic, then store one full cycle
%
%   out = unicycle.ode.runCycles(sys,event,opts)   % one event
%   out = unicycle.ode.runCycles(sys,seq,opts)     % a sequence
%
% Runs unicycle.ode.spinUp (see there for the cycle loop, the slip budget,
% the convergence criteria and the event/sequence input) and then
% integrates one more (super)cycle from the periodic state with dense
% output after every event. Use spinUp directly when only the spun-up
% initial condition is needed, or for one pass through a sequence that
% does not repeat (opts.repeat = false).
%
% INPUT:
%   sys   - structure from unicycle.ode.assembleODE
%   event - structure from unicycle.ode.periodicEvent, or a sequence
%           (T, events), see unicycle.ode.spinUp
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
%           nOut       - number of log-spaced output times after each
%                        event of the stored cycle [161]
%           tFirst     - first output time after an event (s) [1e-8 yr]
%           verbose    - print progress every this many cycles [50], 0 = off
%           (y0, maxCycles, tol*, nConverged, odeOpts, verbose are passed
%           to unicycle.ode.spinUp)
%
% OUTPUT: structure out with
%   converged, nCycles
%   history  - per-cycle table: slipErr, dzeta, drakeDeg, minTplus,
%              nReversed (traction reversed by an event), runtime
%   yStart   - (4Nc x nCycles+1) state at the start of each cycle
%   t,y      - one stored cycle starting from the converged state: t
%              (nt x 1) from the start of the cycle (for one event, the
%              time since the event), y (nt x 4Nc), with slip relative to
%              the start of the cycle. Each event appears at its time
%              twice, before and after the jump (except an event at t = 0,
%              whose pre-event state is the end of the cycle)
%   slipCycle - (Nc x 2) slip of the stored cycle
%   slipTarget - (Nc x 2) slip per cycle of the periodic solution, Δs*
%   tEvents  - (M x 1) event times, sorted
%   zetaPre, thetaPre, zetaPost, thetaPost - (Nc x M) state just before
%              and after each event of the stored cycle
%   event    - the event or sequence as given
%   seq      - the sequence (events sorted by time)
%
% SEE ALSO: unicycle.ode.spinUp, unicycle.ode.periodicEvent,
%           unicycle.ode.sequenceEvent, unicycle.ode.applyEvent,
%           unicycle.ode.odeRateDependent

if nargin < 3, opts = struct(); end
yr = 3.1536e7;
if ~isfield(opts,'nOut') || isempty(opts.nOut), opts.nOut = 161; end
if ~isfield(opts,'tFirst') || isempty(opts.tFirst), opts.tFirst = 1e-8*yr; end
Nc = sys.Nc;

spinOpts = rmfield(opts,intersect(fieldnames(opts),{'nOut','tFirst','tSince','repeat','dense'}));
spinOpts.keepStates = true;
[~,info] = unicycle.ode.spinUp(sys,event,spinOpts);
M = numel(info.seq.events);

% store one full cycle from the periodic state, dense after every event
storeOpts = struct('y0',info.yPre,'repeat',false,'dense',0:M,'nOut',opts.nOut,...
    'tFirst',opts.tFirst,'verbose',0);
if isfield(opts,'odeOpts'), storeOpts.odeOpts = opts.odeOpts; end
[~,cyc] = unicycle.ode.spinUp(sys,info.seq,storeOpts);
d = cyc.dense;

out.converged = info.converged;
out.nCycles = info.nCycles;
out.history = info.history;
out.contraction = info.contraction;
out.yStart = info.yStart;
out.t = vertcat(d.t);
out.y = vertcat(d.y);
out.slipCycle = cyc.slipSince;
out.slipTarget = info.slipTarget;
out.tEvents = cyc.t;
out.zetaPre = cyc.zetaPre;
out.thetaPre = cyc.thetaPre;
out.zetaPost = cyc.zetaPost;
out.thetaPost = cyc.thetaPost;
out.event = event;
out.seq = info.seq;
assert(size(out.y,2) == 4*Nc)

end
