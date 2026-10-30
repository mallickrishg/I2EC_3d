function [vPre,info] = linvisSpinUpNumerical(lin,event,opts)
% LINVISSPINUPNUMERICAL periodic cycle of the linear viscous model by
% time stepping: repeat event + ode45 over Trec until periodic
%
%   [vPre,info] = unicycle.ode.linvisSpinUpNumerical(lin,event,opts)
%
% Numerical counterpart of unicycle.ode.linvisSpinUp, used to verify the
% closed-form periodic solution. Each cycle applies the velocity jump
% (unicycle.ode.linvisEvent) and the coseismic slip on creeping patches,
% then integrates unicycle.ode.odeLinearViscous for Trec. The cycle is
% converged when, for opts.nConverged cycles in a row,
%   max |Δs - Vpl Trec| / (|Vpl| Trec) < opts.tolSlip
%   max |v_pre(k+1) - v_pre(k)| / max |Vpl| < opts.tolV
%
% INPUT:
%   lin   - structure from unicycle.ode.linvisAssemble
%   event - structure with Trec, dtau and optionally slipC, e.g. from
%           unicycle.ode.periodicEvent(lin,Trec,weightC)
%   opts  - (optional) structure with fields (defaults in brackets)
%           v0         - initial pre-event slip rate [Vpl on creeping patches]
%           maxCycles  - [1000]
%           tolSlip    - [1e-8]
%           tolV       - [1e-8]
%           nConverged - [3]
%           odeOpts    - [odeset('RelTol',1e-10,'AbsTol',1e-14)]
%
% OUTPUT:
%   vPre - (2Nc x 1) periodic pre-event slip rate
%   info - structure with converged, nCycles, vPost, slipCycle (Nc x 2),
%          history (slipErr and dv per cycle)
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
Trec = event.Trec;
slipC = zeros(Nc,2);
if isfield(event,'slipC') && ~isempty(event.slipC)
    slipC = event.slipC;
end
VplMag = max(hypot(VplC(:,1),VplC(:,2)));
target = VplC*Trec;
rhs = @(t,y) unicycle.ode.odeLinearViscous(t,y,lin);

v = opts.v0(:);
H = zeros(opts.maxCycles,2);
nOk = 0;
converged = false;
for k = 1:opts.maxCycles
    vPost = unicycle.ode.linvisEvent(lin,v,event.dtau);
    [~,Y] = ode45(rhs,[0,Trec/2,Trec],[slipC(:);vPost],opts.odeOpts);
    yEnd = Y(end,:)';
    ds = [yEnd(1:Nc),yEnd(Nc+1:2*Nc)];
    vNew = yEnd(2*Nc+1:end);
    H(k,:) = [max(hypot(ds(:,1) - target(:,1),ds(:,2) - target(:,2)))/(VplMag*Trec),...
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
info.vPost = unicycle.ode.linvisEvent(lin,vPre,event.dtau);
info.slipCycle = ds;
info.history = array2table(H(1:k,:),'VariableNames',{'slipErr','dv'});

end
