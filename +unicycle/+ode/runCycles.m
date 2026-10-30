function out = runCycles(sys,event,opts)
% RUNCYCLES spin up periodic imposed earthquake cycles until the cycle is
% periodic, then store one full cycle
%
%   out = unicycle.ode.runCycles(sys,event,opts)
%
% Each cycle: apply the event to the pre-event state
% (unicycle.ode.applyEvent), integrate the creeping patches for Trec with
% ode45, and compare with the previous cycle. The slip state is reset at
% the start of every cycle, so its value at the end is the slip of that
% cycle (coseismic + interseismic).
%
% Slip budget: the traction change on creeping patches over one cycle is
%   ΔT = Δτ_event + K_cc Δs_inter + τ̇_load Trec
% so a periodic cycle (ΔT = 0) has a unique slip per cycle
%   Δs* = slip_C - K_cc \ (Δτ_event + τ̇_load Trec)
% With the pinned patches slipping Vpl Trec in every event
% (unicycle.ode.periodicEvent) this is exactly Δs* = Vpl Trec on every
% creeping patch: the slip budget closes everywhere. The cycle is
% converged when
%   max |Δs_c - Δs*| / (|Vpl| Trec)  < opts.tolSlip
%   max |Δzeta| between successive pre-event states < opts.tolZeta
%   max |Δrake| between successive pre-event states < opts.tolRake (deg)
% for opts.nConverged cycles in a row. The rake converges slowest, on the
% time scale |τ|/(k Vpl) (~1e4 yr for f0 = 0.6, a-b = 0.005).
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
% SEE ALSO: unicycle.ode.periodicEvent, unicycle.ode.applyEvent,
%           unicycle.ode.odeRateDependent

if nargin < 3, opts = struct(); end
yr = 3.1536e7;
def = struct('y0',[],'maxCycles',2000,'tolSlip',1e-3,'tolZeta',1e-3,'tolRake',1e-3,...
    'nConverged',3,'odeOpts',odeset('RelTol',1e-8,'AbsTol',1e-8),'nOut',161,...
    'tFirst',1e-8*yr,'verbose',50);
names = fieldnames(def);
for k = 1:numel(names)
    if ~isfield(opts,names{k}) || isempty(opts.(names{k}))
        opts.(names{k}) = def.(names{k});
    end
end

Nc = sys.Nc;
Trec = event.Trec;
if isempty(opts.y0)
    opts.y0 = sys.y0;
end

% the locked region must release Vpl Trec in every event
if isfield(sys,'iP') && ~isempty(sys.iP)
    budget = sys.Vpl(sys.iP,:)*Trec;
    if max(abs(event.slipP(:) - budget(:))) > 1e-9*max(abs(budget(:)))
        warning('unicycle:ode:runCycles:lockedBudget',...
            'runCycles: event slip on pinned patches differs from Vpl Trec; the slip budget of the locked region does not close')
    end
end
if ~isfield(event,'slipC') || isempty(event.slipC)
    event.slipC = zeros(Nc,2);
end

% slip per cycle of the periodic solution
sTarget = event.slipC(:) - sys.Kcc \ (event.dtau(:) + sys.taudotLoad*Trec);
sTarget = [sTarget(1:Nc),sTarget(Nc+1:end)];
VplMag = hypot(sys.Vpl(sys.iC,1),sys.Vpl(sys.iC,2))*Trec;
rhs = @(t,y) unicycle.ode.odeRateDependent(t,y,sys);
iZ = 2*Nc+1:3*Nc;
iT = 3*Nc+1:4*Nc;

y = opts.y0(:);
y(1:2*Nc) = 0;
yStart = zeros(4*Nc,opts.maxCycles+1);
yStart(:,1) = y;
H = zeros(opts.maxCycles,6);
nOk = 0;
converged = false;
for k = 1:opts.maxCycles
    tic
    yMinus = y;
    yPlus = unicycle.ode.applyEvent(sys,yMinus,event.dtau,event.slipC);
    Tm = unicycle.ode.stateTraction(sys,yMinus);
    Tp = unicycle.ode.stateTraction(sys,yPlus);
    [~,Y] = ode45(rhs,[0,Trec/2,Trec],yPlus,opts.odeOpts);
    yEnd = Y(end,:)';

    ds = [yEnd(1:Nc),yEnd(Nc+1:2*Nc)];
    slipErr = max(hypot(ds(:,1) - sTarget(:,1),ds(:,2) - sTarget(:,2))./VplMag);
    dzeta = max(abs(yEnd(iZ) - yMinus(iZ)));
    drake = max(abs(rad2deg(angle(exp(1i*(yEnd(iT) - yMinus(iT)))))));
    H(k,:) = [slipErr,dzeta,drake,min(hypot(Tp(:,1),Tp(:,2))),sum(sum(Tp.*Tm,2) < 0),toc];

    y = yEnd;
    y(1:2*Nc) = 0;
    y(iT) = angle(exp(1i*y(iT)));
    yStart(:,k+1) = y;

    if slipErr < opts.tolSlip && dzeta < opts.tolZeta && drake < opts.tolRake
        nOk = nOk + 1;
    else
        nOk = 0;
    end
    if opts.verbose > 0 && (mod(k,opts.verbose) == 0 || k == 1 || nOk >= opts.nConverged)
        fprintf('cycle %4d: slip err %.2e, dzeta %.2e, drake %.2e deg, min|T+| %.3g MPa, %d reversed, %.2f s\n',...
            k,H(k,1),H(k,2),H(k,3),H(k,4),H(k,5),H(k,6))
    end
    if nOk >= opts.nConverged
        converged = true;
        break
    end
end
nCycles = k;
if ~converged
    warning('unicycle:ode:runCycles:notConverged',...
        'runCycles: not periodic after %d cycles (slip err %.2e, dzeta %.2e, drake %.2e deg)',...
        nCycles,H(k,1),H(k,2),H(k,3))
end

% store one full cycle from the final pre-event state
tout = [0,logspace(log10(opts.tFirst),log10(Trec),opts.nOut)];
yPlus = unicycle.ode.applyEvent(sys,y,event.dtau,event.slipC);
[t,Y] = ode45(rhs,tout,yPlus,opts.odeOpts);

out.converged = converged;
out.nCycles = nCycles;
out.history = array2table(H(1:nCycles,:),'VariableNames',...
    {'slipErr','dzeta','drakeDeg','minTplus','nReversed','runtime'});
out.yStart = yStart(:,1:nCycles+1);
out.t = t;
out.y = Y;
out.slipCycle = [Y(end,1:Nc)',Y(end,Nc+1:2*Nc)'];
out.slipTarget = sTarget;
out.event = event;

end
