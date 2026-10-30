function mesh = prepareMesh(mesh,bigK)
% PREPAREMESH partitions the traction kernel into pinned and creeping
% patches and precomputes the quantities needed by odeRateDependent
%
%   mesh = prepareMesh(mesh,bigK)
%
% INPUT:
%   bigK - [Kss,Kds;Ksd,Kdd] traction kernel (2N x 2N), traction in
%          strike and dip direction due to strike and dip slip
%   mesh - structure with fields (N x 1 unless noted, scalars are expanded)
%          isPinned - logical, patches with imposed v = 0
%          Vpl      - (N x 2) long-term [strike,dip] slip rate
%          a,b      - rate-and-state friction parameters
%          sigma    - effective normal stress
%          f0       - reference friction coefficient
%          v0       - reference slip rate
%          eta      - (optional) radiation damping coefficient, default 0
%
% OUTPUT: mesh with added fields
%   N,Nc,iC - number of patches, number and index of creeping patches
%   Kcc     - (2Nc x 2Nc) kernel between creeping patches
%   taudotLoad - (2Nc x 1) stressing rate on creeping patches due to
%               loading at Vpl with pinned patches held at v = 0
%   y0      - steady-state initial condition at v = Vpl
%
% state vector is y = [slip_s; slip_d; zeta; theta] on creeping patches,
% with zeta = log(|v|/v0) and theta the rake angle from the strike direction
%
% SEE ALSO: odeRateDependent

N = size(mesh.Vpl,1);
mesh.N = N;

% expand scalar parameters
fields = {'a','b','sigma','f0','v0','eta'};
if ~isfield(mesh,'eta')
    mesh.eta = 0;
end
for k = 1:numel(fields)
    if isscalar(mesh.(fields{k}))
        mesh.(fields{k}) = repmat(mesh.(fields{k}),N,1);
    end
end
mesh.isPinned = logical(mesh.isPinned(:));

iC = find(~mesh.isPinned);
mesh.iC = iC;
mesh.Nc = numel(iC);
assert(all(mesh.a(iC) > mesh.b(iC)),'creeping patches must be velocity strengthening (a > b)')

% rows/columns of creeping patches in the 2N-sized system
rowsC = [iC;N+iC];
mesh.Kcc = bigK(rowsC,rowsC);

% τ̇ = K_cc (v_c - Vpl_c) + K_cp (0 - Vpl_p) = K_cc v_c - K(c,:) Vpl
mesh.taudotLoad = -bigK(rowsC,:)*mesh.Vpl(:);

% friction parameters on creeping patches
mesh.Ac = (mesh.a(iC)-mesh.b(iC)).*mesh.sigma(iC);
mesh.f0c = mesh.f0(iC);
mesh.v0c = mesh.v0(iC);
mesh.etac = mesh.eta(iC);
mesh.sigmac = mesh.sigma(iC);
mesh.abc = mesh.a(iC)-mesh.b(iC);

% steady-state initial condition (v = Vpl)
Vs = mesh.Vpl(iC,1);
Vd = mesh.Vpl(iC,2);
zeta0 = log(sqrt(Vs.^2+Vd.^2)./mesh.v0c);
theta0 = atan2(Vd,Vs);
mesh.y0 = [zeros(2*mesh.Nc,1);zeta0;theta0];

end
