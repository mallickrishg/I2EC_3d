function lin = linvisAssemble(rcv,bigK,eta)
% LINVISASSEMBLE set up the linear viscous cycle ODE and its eigen
% decomposition, for the analytic solution
%
%   lin = unicycle.ode.linvisAssemble(rcv,bigK,eta)
%
% On creeping patches the friction is linear viscous and isotropic,
%   τ = τ0 + η v
% and with τ̇ = K_cc v + τ̇_load (pinned patches held at v = 0) the slip
% rate obeys the linear ODE
%   v̇ = A v + D⁻¹ τ̇_load,   A = D⁻¹ K_cc,   D = diag(η, η)
% τ0 drops out. With the steady state v* = -K_cc⁻¹ τ̇_load and u = v - v*,
%   u(t) = exp(A t) u(0) = V exp(Λ t) V⁻¹ u(0)
% using the eigendecomposition A = V Λ V⁻¹. K_cc from centroid collocation
% is not exactly symmetric (~2 % on the generic meshes), so Λ and V may be
% complex (conjugate pairs); results are real to round-off.
%
% INPUT:
%   rcv  - unicycle.geometry.triangleReceiver with Vpl, Vrake and
%          pinnedPosition set by the user (friction properties not needed)
%   bigK - [Kss,Kds;Ksd,Kdd] traction kernel (2N x 2N)
%   eta  - effective viscosity (MPa s/m): scalar, Nc x 1 (creeping
%          patches, in the order of lin.iC) or N x 1 (all patches)
%
% OUTPUT: structure lin with
%   N,Nc,iC,iP  - number of patches, creeping and pinned patch indices
%   Vpl         - (N x 2) long-term [strike,dip] slip rate
%   Kcc,Kcp     - kernel blocks, as in unicycle.ode.assembleODE (so
%                 unicycle.ode.periodicEvent and eventStressChange accept lin)
%   taudotLoad  - (2Nc x 1) stressing rate with pinned patches at v = 0
%   eta         - (Nc x 1) viscosity, d = [eta; eta] (2Nc x 1)
%   vStar       - (2Nc x 1) steady-state slip rate (independent of eta)
%   lambda      - (2Nc x 1) eigenvalues of A (1/s), all Re < 0
%   V, Vinv     - eigenvectors of A and their inverse
%   condV       - condition number of V
%
% SEE ALSO: unicycle.ode.linvisEvolve, unicycle.ode.linvisSpinUp

N = rcv.N;
assert(isequal(size(bigK),[2*N,2*N]),'linvisAssemble: bigK must be 2N x 2N')

Vpl = rcv.Vpl;
Vrake = rcv.Vrake;
assert(~isempty(Vpl) && ~isempty(Vrake),'linvisAssemble: rcv.Vpl and rcv.Vrake must be set')
if isscalar(Vpl), Vpl = repmat(Vpl,N,1); end
if isscalar(Vrake), Vrake = repmat(Vrake,N,1); end

isPinned = false(N,1);
isPinned(rcv.pinnedPosition) = true;
iC = find(~isPinned);
iP = find(isPinned);
Nc = numel(iC);

if isscalar(eta)
    eta = repmat(eta,Nc,1);
elseif numel(eta) == N
    eta = eta(iC);
end
eta = eta(:);
assert(numel(eta) == Nc,'linvisAssemble: eta must be scalar, Nc x 1 or N x 1')
assert(all(eta > 0 & isfinite(eta)),'linvisAssemble: eta must be positive and finite')

lin.N = N;
lin.Nc = Nc;
lin.iC = iC;
lin.iP = iP;
lin.Vpl = [Vpl(:).*cosd(Vrake(:)),Vpl(:).*sind(Vrake(:))];

rowsC = [iC;N+iC];
lin.Kcc = bigK(rowsC,rowsC);
lin.Kcp = bigK(rowsC,[iP;N+iP]);
lin.taudotLoad = -bigK(rowsC,:)*lin.Vpl(:);

lin.eta = eta;
lin.d = [eta;eta];
lin.vStar = -lin.Kcc \ lin.taudotLoad;

[V,L] = eig(lin.Kcc./lin.d);
lin.lambda = diag(L);
lin.V = V;
lin.Vinv = inv(V);
lin.condV = cond(V);
if any(real(lin.lambda) >= 0)
    warning('unicycle:ode:linvisAssemble:notDecaying',...
        'linvisAssemble: %d eigenvalues with Re >= 0; the solution does not relax',sum(real(lin.lambda) >= 0))
end
if lin.condV > 1e10
    warning('unicycle:ode:linvisAssemble:illConditioned',...
        'linvisAssemble: eigenvector matrix is ill conditioned (cond %.2e)',lin.condV)
end

end
