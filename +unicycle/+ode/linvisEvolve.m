function [v,s] = linvisEvolve(lin,v0,t)
% LINVISEVOLVE analytic slip rate and slip of the linear viscous cycle ODE
%
%   [v,s] = unicycle.ode.linvisEvolve(lin,v0,t)
%
% With u = v - v* and the modal amplitudes c = V⁻¹ u(0),
%   v(t) = v* + V (exp(λ t) .* c)
%   s(t) = v* t + V ((exp(λ t) - 1)./λ .* c)
% s is the slip accumulated since t = 0; e^z - 1 is evaluated with
% unicycle.ode.linvisExpm1 to avoid cancellation for small |λ t|. No time
% stepping is involved; the cost is two (2Nc x 2Nc)(2Nc x nt) products.
%
% INPUT:
%   lin - structure from unicycle.ode.linvisAssemble
%   v0  - (2Nc x 1) [strike; dip] slip rate at t = 0 on creeping patches (m/s)
%   t   - (nt x 1) times (s), any order
%
% OUTPUT:
%   v   - (2Nc x nt) slip rate (m/s)
%   s   - (2Nc x nt) slip since t = 0 (m)
%
% SEE ALSO: unicycle.ode.linvisAssemble, unicycle.ode.linvisSpinUp

assert(numel(v0) == 2*lin.Nc,'linvisEvolve: v0 must be 2Nc x 1')
t = t(:)';
c = lin.Vinv*(v0(:) - lin.vStar);
z = lin.lambda.*t;                       % (2Nc x nt)
v = lin.vStar + real(lin.V*(exp(z).*c));
if nargout > 1
    s = lin.vStar.*t + real(lin.V*(unicycle.ode.linvisExpm1(z)./lin.lambda.*c));
end

end
