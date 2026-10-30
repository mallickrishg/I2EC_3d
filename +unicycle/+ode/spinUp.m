function [y,info] = spinUp(sys,event,opts)
% SPINUP initial condition primed by a periodic characteristic earthquake
%
%   [y,info] = unicycle.ode.spinUp(sys,event,opts)
%
% Repeats the characteristic event (unicycle.ode.periodicEvent) every Trec
% until the cycle is periodic, and returns the state a time tSince after
% the last characteristic event. The cycle map is contractive: creep
% diffuses the memory of the initial condition, so the periodic state does
% not depend on where the spin-up starts, and the error shrinks
% geometrically from cycle to cycle (info.contraction). The returned state
% is the initial condition for the instrumental era, on which recorded
% events (unicycle.ode.applyEvent) and creep (unicycle.ode.odeRateDependent
% with the same sys) are then imposed.
%
% Each cycle: apply the event to the pre-event state, integrate the
% creeping patches for Trec with ode45. The slip state is reset at the
% start of every cycle, so its value at the end is the slip of that cycle.
%
% Slip budget: the traction change on creeping patches over one cycle is
%   ΔT = Δτ_event + K_cc Δs_inter + τ̇_load Trec
% so a periodic cycle (ΔT = 0) has a unique slip per cycle
%   Δs* = slip_C - K_cc \ (Δτ_event + τ̇_load Trec)
% With the pinned patches slipping Vpl Trec in every event this is exactly
% Δs* = Vpl Trec on every creeping patch: the slip budget closes
% everywhere. The cycle is converged when
%   max |Δs_c - Δs*| / (|Vpl| Trec)  < opts.tolSlip
%   max |Δzeta| between successive pre-event states < opts.tolZeta
%   max |Δrake| between successive pre-event states < opts.tolRake (deg)
% for opts.nConverged cycles in a row. The rake converges slowest, on the
% time scale |τ|/(k Vpl).
%
% INPUT:
%   sys   - structure from unicycle.ode.assembleODE
%   event - structure from unicycle.ode.periodicEvent
%   opts  - (optional) structure with fields (defaults in brackets)
%           tSince     - time since the last characteristic event at which
%                        to return the state (s), 0 <= tSince <= Trec
%                        [Trec, i.e. just before the next event]
%           y0         - initial pre-event state [sys.y0, v = Vpl at the
%                        plate rake, close to the cycle-averaged state]
%           maxCycles  - [2000]
%           tolSlip    - [1e-3]
%           tolZeta    - [1e-3]
%           tolRake    - [1e-3] (deg)
%           nConverged - [3]
%           odeOpts    - [odeset('RelTol',1e-8,'AbsTol',1e-8)]
%           keepStates - also return the pre-event state of every cycle [false]
%           verbose    - print progress every this many cycles [50], 0 = off
%
% OUTPUT:
%   y    - (4Nc x 1) state [slip_s; slip_d; zeta; theta] a time tSince
%          after the last characteristic event, with slip set to zero
%   info - structure with
%          converged, nCycles, tSince
%          history     - per-cycle table: slipErr, dzeta, drakeDeg,
%                        minTplus, nReversed (traction reversed by the
%                        event), runtime
%          contraction - geometric decrease of slipErr per cycle over the
%                        last 10 cycles (NaN if fewer)
%          slipSince   - (Nc x 2) slip from the last event to tSince,
%                        including coseismic slip on creeping patches
%          yPre        - (4Nc x 1) periodic pre-event state (slip zero)
%          slipTarget  - (Nc x 2) slip per cycle of the periodic solution
%          yStart      - (4Nc x nCycles+1) pre-event states, if keepStates
%          event       - the event structure
%
% SEE ALSO: unicycle.ode.periodicEvent, unicycle.ode.runCycles,
%           unicycle.ode.applyEvent, unicycle.ode.odeRateDependent

if nargin < 3, opts = struct(); end
def = struct('tSince',event.Trec,'y0',[],'maxCycles',2000,'tolSlip',1e-3,'tolZeta',1e-3,...
    'tolRake',1e-3,'nConverged',3,'odeOpts',odeset('RelTol',1e-8,'AbsTol',1e-8),...
    'keepStates',false,'verbose',50);
names = fieldnames(def);
for k = 1:numel(names)
    if ~isfield(opts,names{k}) || isempty(opts.(names{k}))
        opts.(names{k}) = def.(names{k});
    end
end

Nc = sys.Nc;
Trec = event.Trec;
assert(opts.tSince >= 0 && opts.tSince <= Trec,'spinUp: tSince must be in [0, Trec]')
if isempty(opts.y0)
    opts.y0 = sys.y0;
end

% the locked region must release Vpl Trec in every event
if isfield(sys,'iP') && ~isempty(sys.iP)
    budget = sys.Vpl(sys.iP,:)*Trec;
    if max(abs(event.slipP(:) - budget(:))) > 1e-9*max(abs(budget(:)))
        warning('unicycle:ode:spinUp:lockedBudget',...
            'spinUp: event slip on pinned patches differs from Vpl Trec; the slip budget of the locked region does not close')
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
if opts.keepStates
    yStart = zeros(4*Nc,opts.maxCycles+1);
    yStart(:,1) = y;
end
H = zeros(opts.maxCycles,6);
nOk = 0;
converged = false;
for k = 1:opts.maxCycles
    tCycle = tic;
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
    H(k,:) = [slipErr,dzeta,drake,min(hypot(Tp(:,1),Tp(:,2))),sum(sum(Tp.*Tm,2) < 0),toc(tCycle)];

    y = yEnd;
    y(1:2*Nc) = 0;
    y(iT) = angle(exp(1i*y(iT)));
    if opts.keepStates
        yStart(:,k+1) = y;
    end

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
    warning('unicycle:ode:spinUp:notConverged',...
        'spinUp: not periodic after %d cycles (slip err %.2e, dzeta %.2e, drake %.2e deg)',...
        nCycles,H(k,1),H(k,2),H(k,3))
end
yPre = y;

% state a time tSince after the last characteristic event
if opts.tSince < Trec
    yPlus = unicycle.ode.applyEvent(sys,yPre,event.dtau,event.slipC);
    if opts.tSince > 0
        [~,Y] = ode45(rhs,[0,opts.tSince/2,opts.tSince],yPlus,opts.odeOpts);
        y = Y(end,:)';
    else
        y = yPlus;
    end
    slipSince = [y(1:Nc),y(Nc+1:2*Nc)];
    y(1:2*Nc) = 0;
else
    slipSince = ds;
end

info.converged = converged;
info.nCycles = nCycles;
info.tSince = opts.tSince;
info.history = array2table(H(1:nCycles,:),'VariableNames',...
    {'slipErr','dzeta','drakeDeg','minTplus','nReversed','runtime'});
info.contraction = NaN;
if nCycles >= 10
    info.contraction = exp(mean(diff(log(H(nCycles-9:nCycles,1)))));
end
info.slipSince = slipSince;
info.yPre = yPre;
info.slipTarget = sTarget;
if opts.keepStates
    info.yStart = yStart(:,1:nCycles+1);
end
info.event = event;

end
