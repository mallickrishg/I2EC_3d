function [y,v] = steadyState(sys)
% STEADYSTATE state vector of the interseismic steady state with pinned
% patches held at v = 0
%
%   [y,v] = unicycle.ode.steadyState(sys)
%
% The steady state has zero stressing rate on creeping patches,
%   K_cc v + τ̇_load = 0   =>   v = K_cc \ (-τ̇_load)
% and is independent of the friction law. Starting from it instead of
% v = Vpl at the plate rake (sys.y0) avoids the slow rake relaxation,
% whose time scale is |τ|/(k Vpl), i.e. ~1e4 yr for f0 ~ 0.6.
%
% INPUT:
%   sys - structure from unicycle.ode.assembleODE
%
% OUTPUT:
%   y   - state vector [zeros; zeros; zeta; theta] (4Nc x 1)
%   v   - (Nc x 2) [strike,dip] steady-state slip rate (m/s)
%
% SEE ALSO: unicycle.ode.assembleODE, unicycle.ode.odeRateDependent

Nc = sys.Nc;
vv = sys.Kcc \ (-sys.taudotLoad);
v = [vv(1:Nc),vv(Nc+1:end)];
vMag = hypot(v(:,1),v(:,2));
assert(all(vMag > 0),'steadyState: %d creeping patches have zero steady-state slip rate',sum(vMag == 0))

y = [zeros(2*Nc,1);log(vMag./sys.v0);atan2(v(:,2),v(:,1))];

end
