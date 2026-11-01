function [y,info] = spinUp(sys,event,opts)
% SPINUP initial condition primed by a periodic characteristic earthquake
% or a periodic sequence of earthquakes; or one pass through a sequence
%
%   [y,info] = unicycle.ode.spinUp(sys,event,opts)   % one event
%   [y,info] = unicycle.ode.spinUp(sys,seq,opts)     % a sequence
%
% One characteristic event (unicycle.ode.periodicEvent) is a sequence of
% one event at t = 0 with T = Trec; a sequence has T and events (struct
% array, e.g. from unicycle.ode.sequenceEvent), as for
% unicycle.ode.linvisSpinUp.
%
% Repeats the sequence every T until the (super)cycle is periodic, and
% returns the state a time tSince after the start of a periodic
% (super)cycle, i.e. after the start of the last events. The cycle map is contractive: creep diffuses the memory of
% the initial condition, so the periodic state does not depend on where
% the spin-up starts, and the error shrinks geometrically from cycle to
% cycle (info.contraction). The returned state is the initial condition
% for the instrumental era, on which recorded events
% (unicycle.ode.applyEvent) and creep (unicycle.ode.odeRateDependent with
% the same sys) are then imposed.
%
% With opts.repeat = false the sequence is run once from opts.y0, without
% a convergence check. Use this for sequences that do not repeat, e.g. a
% stochastic catalogue spanning many cycles (T = its total length).
%
% Each (super)cycle: integrate the creeping patches with ode45 from one
% event to the next and apply each event (unicycle.ode.applyEvent) at its
% time. The slip state is reset at the start of every interval and summed
% separately, so its total is the slip of that (super)cycle.
%
% Slip budget: the traction change on creeping patches over one cycle is
%   ΔT = Σ Δτ_k + K_cc Δs_inter + τ̇_load T
% so a periodic cycle (ΔT = 0) has a unique slip per cycle
%   Δs* = Σ slip_C,k - K_cc \ (Σ Δτ_k + τ̇_load T)
% With the pinned patches slipping Vpl T over the events of one cycle this
% is exactly Δs* = Vpl T on every creeping patch: the slip budget closes
% everywhere, wherever the events rupture. The cycle is converged when
%   max |Δs_c - Δs*| / (|Vpl| T)  < opts.tolSlip
%   max |Δzeta| between successive pre-cycle states < opts.tolZeta
%   max |Δrake| between successive pre-cycle states < opts.tolRake (deg)
% for opts.nConverged cycles in a row. The rake converges slowest, on the
% time scale |τ|/(k Vpl).
%
% INPUT:
%   sys   - structure from unicycle.ode.assembleODE
%   event - one event: structure with Trec, dtau (2Nc x 1) and optionally
%           slipC (Nc x 2) and slipP (Np x 2), e.g.
%           unicycle.ode.periodicEvent(sys,Trec,weightC)
%   seq   - or a sequence: structure with
%           T      - (super)cycle period (s)
%           events - struct array, one per event, with fields
%                    t     - event time, 0 <= t < T (s)
%                    slipP - (Np x 2) coseismic slip on pinned patches (m)
%                    slipC - (optional) (Nc x 2) coseismic slip on creeping
%                            patches (m)
%                    dtau  - (optional) (2Nc x 1) traction change on
%                            creeping patches (MPa) [K_cp slip_P + K_cc slip_C]
%                    e.g. from unicycle.ode.sequenceEvent
%   opts  - (optional) structure with fields (defaults in brackets)
%           tSince     - time since the start of the periodic (super)cycle
%                        (or of the single pass) at which to return the
%                        state (s), 0 <= tSince <= T; an event at t_k is
%                        included at t = t_k [T, i.e. just before the next
%                        (super)cycle]
%           y0         - initial state at the start of the first cycle
%                        [sys.y0, v = Vpl at the plate rake, close to the
%                        cycle-averaged state]
%           repeat     - repeat the sequence until periodic [true]; false
%                        runs it once from y0 to tSince
%           maxCycles  - [2000]
%           tolSlip    - [1e-3]
%           tolZeta    - [1e-3]
%           tolRake    - [1e-3] (deg)
%           nConverged - [3]
%           odeOpts    - [odeset('RelTol',1e-8,'AbsTol',1e-8)]
%           keepStates - also return the state at the start of every cycle [false]
%           dense      - intervals of the last pass with dense, log-spaced
%                        output: k = 1..M is the interval after event k (to
%                        the next event or tSince), 0 the interval before
%                        the first event; only with repeat = false (use
%                        unicycle.ode.runCycles for the periodic cycle) []
%           nOut       - log-spaced output times per dense interval [121]
%           tFirst     - first dense output time after an event (s) [1e-8 yr]
%           verbose    - print progress every this many cycles (repeat) or
%                        events (repeat = false) [50], 0 = off
%
% OUTPUT:
%   y    - (4Nc x 1) state [slip_s; slip_d; zeta; theta] a time tSince
%          after yPre, with slip set to zero
%   info - structure with
%          converged, nCycles, tSince (converged = true and nCycles = 1
%                      for repeat = false)
%          history     - per-cycle table: slipErr, dzeta, drakeDeg,
%                        minTplus, nReversed (traction reversed by an
%                        event), runtime
%          contraction - geometric decrease of slipErr per cycle over the
%                        last 10 cycles (NaN if fewer)
%          slipSince   - (Nc x 2) slip from yPre to tSince, including
%                        coseismic slip on creeping patches
%          yPre        - (4Nc x 1) state at the start of a (super)cycle,
%                        just before any event at t = 0: the periodic state
%                        if converged, y0 for repeat = false (slip zero)
%          slipTarget  - (Nc x 2) slip per cycle of the periodic solution
%          yStart      - (4Nc x nCycles+1) states at the start of each
%                        cycle, if keepStates
%          t           - (M x 1) event times, sorted; events are reordered
%                        by time in the outputs below
%          zetaPre, thetaPre   - (Nc x M) state just before and just after
%          zetaPost, thetaPost   each event of the last full cycle (or of
%                                the single pass; NaN for events after tSince)
%          slipPre     - (Nc x 2 x M) slip from the start of that cycle
%                        to just before each event
%          dense       - struct array, one per opts.dense interval, with
%                        fields interval, t (from the start of the cycle)
%                        and y (nt x 4Nc, slip from the start of the cycle)
%          seq         - the sequence (events sorted by time)
%          event       - the input event or sequence
%
% SEE ALSO: unicycle.ode.periodicEvent, unicycle.ode.sequenceEvent,
%           unicycle.ode.runCycles, unicycle.ode.applyEvent,
%           unicycle.ode.odeRateDependent, unicycle.ode.linvisSpinUp

if nargin < 3, opts = struct(); end
Nc = sys.Nc;
Np = 0;
if isfield(sys,'iP'), Np = numel(sys.iP); end
seq = asSequence(sys,event);
T = seq.T;
M = numel(seq.events);
def = struct('tSince',T,'y0',[],'repeat',true,'maxCycles',2000,'tolSlip',1e-3,'tolZeta',1e-3,...
    'tolRake',1e-3,'nConverged',3,'odeOpts',odeset('RelTol',1e-8,'AbsTol',1e-8),...
    'keepStates',false,'dense',[],'nOut',121,'tFirst',1e-8*3.1536e7,'verbose',50);
names = fieldnames(def);
for k = 1:numel(names)
    if ~isfield(opts,names{k}) || isempty(opts.(names{k}))
        opts.(names{k}) = def.(names{k});
    end
end
assert(opts.tSince >= 0 && opts.tSince <= T,'spinUp: tSince must be in [0, T]')
assert(~opts.repeat || isempty(opts.dense),'spinUp: dense output needs repeat = false (or use runCycles)')
assert(all(ismember(opts.dense,0:M)),'spinUp: dense must index intervals 0..M')
if isempty(opts.y0)
    opts.y0 = sys.y0;
end

% the locked region must release Vpl T over the events of one cycle
if Np > 0
    slipP = zeros(Np,2);
    for k = 1:M
        if ~isempty(seq.events(k).slipP), slipP = slipP + seq.events(k).slipP; end
    end
    budget = sys.Vpl(sys.iP,:)*T;
    if max(abs(slipP(:) - budget(:))) > 1e-9*max(abs(budget(:)))
        warning('unicycle:ode:spinUp:lockedBudget',...
            'spinUp: event slip on pinned patches differs from Vpl T; the slip budget of the locked region does not close')
    end
end

% slip per cycle of the periodic solution
dtauSum = sum([seq.events.dtau],2);
slipCSum = sum(cat(3,seq.events.slipC),3);
sTarget = slipCSum(:) - sys.Kcc \ (dtauSum + sys.taudotLoad*T);
sTarget = [sTarget(1:Nc),sTarget(Nc+1:end)];
VplMag = hypot(sys.Vpl(sys.iC,1),sys.Vpl(sys.iC,2))*T;
iZ = 2*Nc+1:3*Nc;
iT = 3*Nc+1:4*Nc;

y = opts.y0(:);
y(1:2*Nc) = 0;
verboseEvents = 0;
if ~opts.repeat
    opts.maxCycles = 1;
    verboseEvents = opts.verbose;
end
if opts.keepStates
    yStart = zeros(4*Nc,opts.maxCycles+1);
    yStart(:,1) = y;
end
H = zeros(opts.maxCycles,6);
nOk = 0;
converged = ~opts.repeat;
for k = 1:opts.maxCycles
    tCycle = tic;
    yMinus = y;
    tStop = T;
    if ~opts.repeat, tStop = opts.tSince; end
    [yEnd,ds,rec] = onePass(sys,seq,yMinus,tStop,opts,verboseEvents);

    slipErr = max(hypot(ds(:,1) - sTarget(:,1),ds(:,2) - sTarget(:,2))./VplMag);
    dzeta = max(abs(yEnd(iZ) - yMinus(iZ)));
    drake = max(abs(rad2deg(angle(exp(1i*(yEnd(iT) - yMinus(iT)))))));
    H(k,:) = [slipErr,dzeta,drake,rec.minTplus,rec.nReversed,toc(tCycle)];

    y = yEnd;
    if opts.keepStates
        yStart(:,k+1) = y;
    end
    if ~opts.repeat, break, end

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

slipSince = ds;
if opts.repeat
    % the periodic state at the start of a cycle, and the state a time
    % tSince into that cycle
    yPre = y;
    if opts.tSince < T
        [y,slipSince] = onePass(sys,seq,yPre,opts.tSince,opts,0);
    end
else
    yPre = opts.y0(:);
    yPre(1:2*Nc) = 0;
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
info.t = [seq.events.t]';
info.zetaPre = rec.zetaPre;
info.thetaPre = rec.thetaPre;
info.zetaPost = rec.zetaPost;
info.thetaPost = rec.thetaPost;
info.slipPre = rec.slipPre;
info.dense = rec.dense;
info.seq = seq;
info.event = event;

end

function seq = asSequence(sys,event)
% one event or a sequence -> sequence with sorted events and all fields
Nc = sys.Nc;
Np = 0;
if isfield(sys,'iP'), Np = numel(sys.iP); end
if isfield(event,'events')
    seq = event;
else
    assert(isscalar(event) && isfield(event,'Trec'),...
        'spinUp: give an event (Trec, dtau) or a sequence (T, events)')
    ev = struct('t',0,'dtau',event.dtau,'slipP',[],'slipC',[]);
    if isfield(event,'slipP'), ev.slipP = event.slipP; end
    if isfield(event,'slipC'), ev.slipC = event.slipC; end
    seq = struct('T',event.Trec,'events',ev);
end
T = seq.T;
assert(isscalar(T) && T > 0,'spinUp: seq.T must be a positive scalar')
ev = seq.events(:);
assert(~isempty(ev),'spinUp: seq.events is empty')
tk = [ev.t]';
assert(numel(tk) == numel(ev) && all(tk >= 0 & tk < T),'spinUp: event times must be in [0, T)')
[~,order] = sort(tk);
ev = ev(order);
out = repmat(struct('t',0,'dtau',zeros(2*Nc,1),'slipP',zeros(Np,2),'slipC',zeros(Nc,2)),numel(ev),1);
for k = 1:numel(ev)
    out(k).t = ev(k).t;
    if isfield(ev,'slipP') && ~isempty(ev(k).slipP), out(k).slipP = ev(k).slipP; end
    if isfield(ev,'slipC') && ~isempty(ev(k).slipC), out(k).slipC = ev(k).slipC; end
    if isfield(ev,'dtau') && ~isempty(ev(k).dtau)
        out(k).dtau = ev(k).dtau(:);
    else
        out(k).dtau = sys.Kcc*out(k).slipC(:);
        if Np > 0, out(k).dtau = out(k).dtau + sys.Kcp*out(k).slipP(:); end
    end
    assert(numel(out(k).dtau) == 2*Nc,'spinUp: dtau must be 2Nc x 1')
end
seq.events = out;
end

function [y,slip,rec] = onePass(sys,seq,y,tStop,opts,verbose)
% integrate from t = 0 (state y, slip zero) to tStop through the events
% with t_k <= tStop; slip is the total slip of the pass
Nc = sys.Nc;
iS = 1:2*Nc;
iZ = 2*Nc+1:3*Nc;
iT = 3*Nc+1:4*Nc;
rhs = @(t,y) unicycle.ode.odeRateDependent(t,y,sys);
ev = seq.events;
M = numel(ev);
tk = [ev.t]';
nDo = nnz(tk <= tStop);

rec.zetaPre = nan(Nc,M); rec.thetaPre = nan(Nc,M);
rec.zetaPost = nan(Nc,M); rec.thetaPost = nan(Nc,M);
rec.slipPre = nan(Nc,2,M);
rec.minTplus = Inf;
rec.nReversed = 0;
rec.dense = struct('interval',{},'t',{},'y',{});

slip = zeros(Nc,2);
y(iS) = 0;
tNow = 0;
for k = 1:nDo+1
    % creep from tNow to the next event (or tStop); interval k-1
    if k <= nDo, tNext = tk(k); else, tNext = tStop; end
    dt = tNext - tNow;
    isDense = any(opts.dense == k-1);
    if dt > 0
        if isDense
            tout = [0,logspace(log10(min(opts.tFirst,dt/2)),log10(dt),opts.nOut)];
        else
            tout = [0,dt/2,dt];
        end
        [tt,Y] = ode45(rhs,tout,y,opts.odeOpts);
        y = Y(end,:)';
        if isDense
            Y(:,iS) = Y(:,iS) + reshape(slip,1,[]);
            tt = tNow + tt;
            tt([1,end]) = [tNow,tNext];
            rec.dense(end+1) = struct('interval',k-1,'t',tt,'y',Y);
        end
    end
    slip = slip + reshape(y(iS),Nc,2);
    y(iS) = 0;
    y(iT) = angle(exp(1i*y(iT)));
    tNow = tNext;
    if k > nDo, break, end

    % event k
    rec.zetaPre(:,k) = y(iZ);
    rec.thetaPre(:,k) = y(iT);
    rec.slipPre(:,:,k) = slip;
    Tm = unicycle.ode.stateTraction(sys,y);
    y = unicycle.ode.applyEvent(sys,y,ev(k).dtau,ev(k).slipC);
    Tp = unicycle.ode.stateTraction(sys,y);
    rec.minTplus = min(rec.minTplus,min(hypot(Tp(:,1),Tp(:,2))));
    rec.nReversed = max(rec.nReversed,sum(sum(Tp.*Tm,2) < 0));
    slip = slip + ev(k).slipC;
    y(iS) = 0;
    rec.zetaPost(:,k) = y(iZ);
    rec.thetaPost(:,k) = y(iT);
    if verbose > 0 && (mod(k,verbose) == 0 || k == 1)
        fprintf('event %5d of %d at %10.1f yr\n',k,M,tk(k)/3.1536e7)
    end
end
end
