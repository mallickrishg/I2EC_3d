function dydt = odeRateDependent(~,y,sys)
% ODERATEDEPENDENT evolution of slip on creeping patches governed by
% steady-state rate-dependent friction, with slip velocity collinear with
% shear traction
%
%   dydt = unicycle.ode.odeRateDependent(t,y,sys)
%
% friction:  |τ| = F(|v|) (see frictionLaw),  v ∥ τ
% elasticity: τ̇ = K_cc v_c + τ̇_load   (see assembleODE)
%
% with zeta = log(|v|/v0) and rake theta, project τ̇ onto the slip
% direction and its normal:
%   parallel      τ̇·v̂  = (dF/dζ + η|v|) ζ̇    -- acceleration, dF/dζ = (a-b)σ
%   perpendicular τ̇·v̂⊥ = (|τ| + η|v|) θ̇     -- rotation, controlled by f0/(a-b)
%
% y = [slip_s; slip_d; zeta; theta] on creeping patches
%
% SEE ALSO: unicycle.ode.assembleODE, unicycle.ode.frictionLaw

Nc = sys.Nc;
zeta = y(2*Nc+1:3*Nc);
theta = y(3*Nc+1:4*Nc);

% slip velocity
v = sys.v0.*exp(zeta);
c = cos(theta);
s = sin(theta);
vs = v.*c;
vd = v.*s;

% traction rate in strike and dip direction
taudot = sys.Kcc*[vs;vd] + sys.taudotLoad;
taudot_s = taudot(1:Nc);
taudot_d = taudot(Nc+1:end);

% project onto slip direction and its normal
taudot_par = c.*taudot_s + s.*taudot_d;
taudot_perp = -s.*taudot_s + c.*taudot_d;

[tauMag,dtaudzeta] = unicycle.ode.frictionLaw(sys,zeta);

dzeta = taudot_par./(dtaudzeta + sys.eta.*v);
dtheta = taudot_perp./(tauMag + sys.eta.*v);

dydt = [vs;vd;dzeta;dtheta];

end
