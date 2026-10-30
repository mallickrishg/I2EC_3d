function sys = assembleODE(rcv,bigK,law)
% ASSEMBLEODE one-time setup of the interseismic ODE for rate-dependent
% friction: partitions the kernel into pinned and creeping patches and
% precomputes the loading rate and the steady-state initial condition
%
%   sys = unicycle.ode.assembleODE(rcv,bigK,law)
%
% INPUT:
%   rcv  - unicycle.geometry.triangleReceiver with the following properties
%          set by the user (N x 1 or scalar; units m, s, MPa)
%          a,b            - rate-and-state friction parameters
%          sigma          - effective normal stress
%          mu0            - reference friction coefficient (f0)
%          Vo             - reference slip rate (v0)
%          Vs             - shear-wave speed, radiation damping G/(2Vs)
%          Vpl,Vrake      - long-term slip rate and its rake (degrees)
%          pinnedPosition - indices of patches held at v = 0
%   bigK - [Kss,Kds;Ksd,Kdd] traction kernel (2N x 2N), see unicycle.greens.computeKernels
%   law  - (optional) steady-state friction law, 'asinh' (default) or 'log'
%          see unicycle.ode.frictionLaw
%
% OUTPUT: structure sys with
%   N,Nc,iC    - number of patches, number and index of creeping patches
%   Kcc        - (2Nc x 2Nc) kernel between creeping patches
%   taudotLoad - (2Nc x 1) stressing rate on creeping patches due to
%                loading at Vpl with pinned patches held at v = 0
%   Vpl        - (N x 2) long-term [strike,dip] slip rate
%   law        - friction law
%   ab,abSigma,rho,f0,v0,eta,sigma - friction properties on creeping
%                patches, rho = f0/(a-b)
%   y0         - steady-state initial condition at v = Vpl
%
% state vector is y = [slip_s; slip_d; zeta; theta] on creeping patches,
% with zeta = log(|v|/v0) and theta the rake from the strike direction
%
% SEE ALSO: unicycle.ode.odeRateDependent, unicycle.greens.computeKernels,
%           unicycle.ode.frictionLaw

if nargin < 3 || isempty(law), law = 'asinh'; end
assert(any(strcmp(law,{'asinh','log'})),'assembleODE: law must be ''asinh'' or ''log''')

N = rcv.N;
assert(isequal(size(bigK),[2*N,2*N]),'assembleODE: bigK must be 2N x 2N')

% read properties from the mesh, expanding scalars
p = struct();
names = {'a','b','sigma','mu0','Vo','Vs','Vpl','Vrake'};
for k = 1:numel(names)
    val = rcv.(names{k});
    assert(~isempty(val),'assembleODE: rcv.%s is not set',names{k})
    if isscalar(val)
        val = repmat(val,N,1);
    end
    assert(numel(val) == N,'assembleODE: rcv.%s must be N x 1 or scalar',names{k})
    p.(names{k}) = val(:);
end
if any(p.Vs == 3e3*3.1536e7)
    warning('unicycle:ode:assembleODE:defaultVs','assembleODE: rcv.Vs still holds the m/yr constructor default')
end
assert(any(p.Vpl > 0),'assembleODE: rcv.Vpl is zero everywhere')

isPinned = false(N,1);
isPinned(rcv.pinnedPosition) = true;
iC = find(~isPinned);
Nc = numel(iC);
assert(all(p.a(iC) > p.b(iC)),'assembleODE: creeping patches must be velocity strengthening (a > b)')

sys.N = N;
sys.iC = iC;
sys.Nc = Nc;
sys.law = law;

% long-term slip rate in strike and dip components
sys.Vpl = [p.Vpl.*cosd(p.Vrake),p.Vpl.*sind(p.Vrake)];

% rows/columns of creeping patches in the 2N-sized system
rowsC = [iC;N+iC];
sys.Kcc = bigK(rowsC,rowsC);

% loading on creeping patches, τ̇ = K (v - Vpl), with v = 0 on pinned
% patches:
%   τ̇_c = K_cc (v_c - Vpl_c) + K_cp (0 - Vpl_p)
%       = K_cc v_c - [K_cc, K_cp] Vpl
% the pinned columns therefore only enter through τ̇_load
sys.taudotLoad = -bigK(rowsC,:)*sys.Vpl(:);

% friction properties on creeping patches
sys.sigma = p.sigma(iC);
sys.ab = p.a(iC)-p.b(iC);
sys.abSigma = sys.ab.*sys.sigma;
sys.f0 = p.mu0(iC);
sys.rho = sys.f0./sys.ab;
sys.v0 = p.Vo(iC);
sys.eta = rcv.earthModel.G./(2*p.Vs(iC));

% steady-state initial condition (v = Vpl)
Vs = sys.Vpl(iC,1);
Vd = sys.Vpl(iC,2);
zeta0 = log(sqrt(Vs.^2+Vd.^2)./sys.v0);
theta0 = atan2(Vd,Vs);
sys.y0 = [zeros(2*Nc,1);zeta0;theta0];

end
