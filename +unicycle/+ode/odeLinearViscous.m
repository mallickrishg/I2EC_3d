function dydt = odeLinearViscous(~,y,lin)
% ODELINEARVISCOUS right-hand side of the cycle ODE for linear viscous
% friction, for numerical integration (e.g. ode45)
%
%   dydt = unicycle.ode.odeLinearViscous(t,y,lin)
%
% With τ = τ0 + η v on creeping patches and τ̇ = K_cc v + τ̇_load,
%   η v̇ = K_cc v + τ̇_load,   ṡ = v
% The state is y = [slip_s; slip_d; v_s; v_d] on creeping patches
% (4Nc x 1). This is the numerical counterpart of the closed-form solution
% in unicycle.ode.linvisEvolve.
%
% INPUT:
%   lin - structure from unicycle.ode.linvisAssemble (uses Kcc,
%         taudotLoad, d and Nc)
%
% SEE ALSO: unicycle.ode.linvisEvolve, unicycle.ode.linvisSpinUpNumerical

n = 2*lin.Nc;
v = y(n+1:2*n);
dydt = [v;(lin.Kcc*v + lin.taudotLoad)./lin.d];

end
