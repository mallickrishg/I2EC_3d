function dydt = odeRateDependent(~,y,mesh)
% ODERATEDEPENDENT evolution of slip on creeping patches governed by
% steady-state rate-dependent friction, with slip velocity collinear with
% shear traction
%
%   dydt = odeRateDependent(t,y,mesh)
%
% friction:  |τ| = σ (f0 + (a-b) log(|v|/v0)),  v ∥ τ
% elasticity: τ̇ = K_cc v_c + τ̇_load   (see prepareMesh)
%
% with zeta = log(|v|/v0) and rake theta, project τ̇ onto the slip
% direction and its normal:
%   parallel      τ̇·v̂  = ((a-b)σ + η|v|) ζ̇     -- acceleration, f0 and v0 drop out
%   perpendicular τ̇·v̂⊥ = (|τ| + η|v|) θ̇         -- rotation, controlled by f0/(a-b)
%
% y = [slip_s; slip_d; zeta; theta] on creeping patches
%
% SEE ALSO: prepareMesh

Nc = mesh.Nc;
zeta = y(2*Nc+1:3*Nc);
theta = y(3*Nc+1:4*Nc);

% slip velocity
v = mesh.v0c.*exp(zeta);
c = cos(theta);
s = sin(theta);
vs = v.*c;
vd = v.*s;

% traction rate in strike and dip direction
taudot = mesh.Kcc*[vs;vd] + mesh.taudotLoad;
taudot_s = taudot(1:Nc);
taudot_d = taudot(Nc+1:end);

% project onto slip direction and its normal
taudot_par = c.*taudot_s + s.*taudot_d;
taudot_perp = -s.*taudot_s + c.*taudot_d;

% traction magnitude
tauMag = mesh.sigmac.*(mesh.f0c + mesh.abc.*zeta);

dzeta = taudot_par./(mesh.Ac + mesh.etac.*v);
dtheta = taudot_perp./(tauMag + mesh.etac.*v);

dydt = [vs;vd;dzeta;dtheta];

end
