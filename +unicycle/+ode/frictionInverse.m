function zeta = frictionInverse(sys,tauMag,useDamping)
% FRICTIONINVERSE zeta = log(|v|/v0) for a given shear traction magnitude,
% inverse of unicycle.ode.frictionLaw, optionally with radiation damping
%
%   zeta = unicycle.ode.frictionInverse(sys,tauMag,useDamping)
%
% without damping (useDamping = false, or sys.eta = 0), solves
%   |τ| = F(v)
% explicitly:
%   'log'   zeta = (|τ|/σ - f0)/(a-b)
%   'asinh' v = 2 v0 exp(-f0/(a-b)) sinh(|τ|/((a-b)σ))
%
% with damping (default), the traction is shared between friction and
% radiation damping, and v solves the implicit equation
%   g(zeta) = F(zeta) + η v0 exp(zeta) - |τ| = 0
% g is increasing and convex in zeta for both laws, so the root is unique
% and Newton's method started from a point with g >= 0 converges
% monotonically from above. Starting points: the undamped solution
% (g = η v >= 0) or the damping cap v = |τ|/η when F >= 0 there.
%
% SEE ALSO: unicycle.ode.frictionLaw

if nargin < 3 || isempty(useDamping), useDamping = true; end

% undamped inverse
switch sys.law
    case 'log'
        zeta = (tauMag./sys.sigma - sys.f0)./sys.ab;
    case 'asinh'
        % log(2 sinh(s)) = s + log(1 - exp(-2s)), s > 0
        s = tauMag./sys.abSigma;
        zeta = s + log(-expm1(-2*s)) - sys.rho;
    otherwise
        error('frictionInverse: unknown law %s',sys.law)
end

if ~useDamping || ~isfield(sys,'eta') || all(sys.eta == 0)
    return
end
eta = sys.eta.*ones(size(zeta));
v0 = sys.v0.*ones(size(zeta));

% start from the damping cap where it is closer and still has g >= 0
damped = eta > 0;
zetaCap = zeta;
zetaCap(damped) = log(tauMag(damped)./(eta(damped).*v0(damped)));
gCap = g(sys,zetaCap,tauMag,eta);
useCap = damped & zetaCap < zeta & gCap >= 0;
zeta(useCap) = zetaCap(useCap);

% Newton iterations
for iter = 1:100
    [gz,dg] = g(sys,zeta,tauMag,eta);
    step = gz./dg;
    zeta = zeta - step;
    if all(abs(step) < 1e-13*max(1,abs(zeta)))
        return
    end
end
error('frictionInverse: Newton iterations did not converge')

end

function [gz,dg] = g(sys,zeta,tauMag,eta)
% residual of traction balance and its derivative with respect to zeta
[F,dF] = unicycle.ode.frictionLaw(sys,zeta);
etav = eta.*sys.v0.*exp(zeta);
gz = F + etav - tauMag;
dg = dF + etav;
end
