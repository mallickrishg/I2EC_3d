function [vPre,info] = linvisSpinUpNumerical(lin,event,opts)
% LINVISSPINUPNUMERICAL periodic cycle of the linear viscous model by
% time stepping: repeat event(s) + ode45 over the (super)cycle until
% periodic
%
%   [vPre,info] = unicycle.ode.linvisSpinUpNumerical(lin,event,opts)
%   [vPre,info] = unicycle.ode.linvisSpinUpNumerical(lin,seq,opts)
%
% Numerical counterpart of unicycle.ode.linvisSpinUp, used to verify the
% closed-form periodic solution.
% Each cycle applies, at each event time, the velocity jump
% (unicycle.ode.linvisEvent) and the coseismic slip on creeping patches,
% and integrates unicycle.ode.odeLinearViscous in between, over a period
% T (Trec for a single event). The cycle is converged when, for
% opts.nConverged cycles in a row,
%   max |Δs - Vpl T| / (|Vpl| T) < opts.tolSlip
%   max |v(0⁻, k+1) - v(0⁻, k)| / max |Vpl| < opts.tolV
%
% INPUT:
%   lin   - structure from unicycle.ode.linvisAssemble
%   event - structure with Trec, dtau and optionally slipC, e.g. from
%           unicycle.ode.periodicEvent(lin,Trec,weightC), the event at t = 0
%   seq   - or a sequence: structure with T and events (struct array with
%           t, dtau, optionally slipC), as for linvisSpinUp
%   opts  - (optional) structure with fields (defaults in brackets)
%           v0         - initial slip rate at t = 0⁻ [Vpl on creeping patches]
%           maxCycles  - [1000]
%           tolSlip    - [1e-8]
%           tolV       - [1e-8]
%           nConverged - [3]
%           odeOpts    - [odeset('RelTol',1e-10,'AbsTol',1e-14)]
%
% OUTPUT:
%   vPre - (2Nc x 1) periodic slip rate at t = 0⁻ (before the event)
%   info - structure with converged, nCycles, vPre and vPost (2Nc x M, just
%          before and after each event, sorted by time), slipCycle
%          (Nc x 2), history (slipErr and dv per cycle)
%
% SEE ALSO: unicycle.ode.linvisSpinUp, unicycle.ode.odeLinearViscous

if nargin < 3, opts = struct(); end
Nc = lin.Nc;
VplC = lin.Vpl(lin.iC,:);
def = struct('v0',reshape(VplC,[],1),'maxCycles',1000,'tolSlip',1e-8,'tolV',1e-8,...
    'nConverged',3,'odeOpts',odeset('RelTol',1e-10,'AbsTol',1e-14));
names = fieldnames(def);
for k = 1:numel(names)
    if ~isfield(opts,names{k}) || isempty(opts.(names{k}))
        opts.(names{k}) = def.(names{k});
    end
end
if isfield(event,'events')
    T = event.T;
    [~,order] = sort([event.events.t]);
    ev = event.events(order);
else
    T = event.Trec;
    ev = event;
    ev.t = 0;
end
M = numel(ev);
slipC = zeros(2*Nc,M);
for j = 1:M
    if isfield(ev,'slipC') && ~isempty(ev(j).slipC)
        slipC(:,j) = ev(j).slipC(:);
    end
end
VplMag = max(hypot(VplC(:,1),VplC(:,2)));
target = VplC*T;
rhs = @(t,y) unicycle.ode.odeLinearViscous(t,y,lin);

v = opts.v0(:);
H = zeros(opts.maxCycles,2);
nOk = 0;
converged = false;
vPreE = zeros(2*Nc,M);
vPostE = zeros(2*Nc,M);
for k = 1:opts.maxCycles
    y = [zeros(2*Nc,1);v];
    tPrev = 0;
    for j = 1:M
        y = step(rhs,y,tPrev,ev(j).t,opts.odeOpts);
        vPreE(:,j) = y(2*Nc+1:end);
        y(1:2*Nc) = y(1:2*Nc) + slipC(:,j);
        y(2*Nc+1:end) = unicycle.ode.linvisEvent(lin,y(2*Nc+1:end),ev(j).dtau);
        vPostE(:,j) = y(2*Nc+1:end);
        tPrev = ev(j).t;
    end
    yEnd = step(rhs,y,tPrev,T,opts.odeOpts);
    ds = [yEnd(1:Nc),yEnd(Nc+1:2*Nc)];
    vNew = yEnd(2*Nc+1:end);
    H(k,:) = [max(hypot(ds(:,1) - target(:,1),ds(:,2) - target(:,2)))/(VplMag*T),...
              max(abs(vNew - v))/VplMag];
    v = vNew;
    if H(k,1) < opts.tolSlip && H(k,2) < opts.tolV
        nOk = nOk + 1;
    else
        nOk = 0;
    end
    if nOk >= opts.nConverged
        converged = true;
        break
    end
end
if ~converged
    warning('unicycle:ode:linvisSpinUpNumerical:notConverged',...
        'linvisSpinUpNumerical: not periodic after %d cycles (slip err %.2e, dv %.2e)',k,H(k,1),H(k,2))
end

vPre = v;
info.converged = converged;
info.nCycles = k;
info.vPre = vPreE;
info.vPost = vPostE;
info.slipCycle = ds;
info.history = array2table(H(1:k,:),'VariableNames',{'slipErr','dv'});

end

function y = step(rhs,y,t0,t1,odeOpts)
% integrate from t0 to t1 (nothing to do between simultaneous events)
if t1 > t0
    [~,Y] = ode45(rhs,[t0,(t0 + t1)/2,t1],y,odeOpts);
    y = Y(end,:)';
end
end
