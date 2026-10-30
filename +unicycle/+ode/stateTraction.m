function T = stateTraction(sys,y)
% STATETRACTION shear traction vector on creeping patches implied by the
% ODE state
%
%   T = unicycle.ode.stateTraction(sys,y)
%
% The traction balances friction and radiation damping and is collinear
% with the slip rate:
%   |T| = F(v) + η v,   T = |T| [cos θ, sin θ]
%
% INPUT:
%   sys - structure from unicycle.ode.assembleODE
%   y   - state vector [slip_s; slip_d; zeta; theta] (4Nc x 1)
%
% OUTPUT:
%   T   - (Nc x 2) [strike,dip] shear traction (MPa)
%
% SEE ALSO: unicycle.ode.applyEvent, unicycle.ode.frictionLaw

Nc = sys.Nc;
zeta = y(2*Nc+1:3*Nc);
theta = y(3*Nc+1:4*Nc);
tauMag = unicycle.ode.frictionLaw(sys,zeta) + sys.eta.*sys.v0.*exp(zeta);
if any(tauMag < 0)
    error('unicycle:ode:stateTraction:negativeTraction',...
        ['stateTraction: %d patches have F(v) + eta v < 0 (log law below ',...
         'v0 exp(-f0/(a-b))); the traction direction is undefined, use the asinh law'],...
        sum(tauMag < 0))
end
T = tauMag.*[cos(theta),sin(theta)];

end
