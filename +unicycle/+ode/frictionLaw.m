function [tauMag,dtaudzeta] = frictionLaw(sys,zeta)
% FRICTIONLAW magnitude of shear traction at steady-state sliding and its
% derivative with respect to zeta = log(|v|/v0)
%
%   [tauMag,dtaudzeta] = unicycle.ode.frictionLaw(sys,zeta)
%
% sys.law selects
%   'log'   |τ| = σ (f0 + (a-b) log(v/v0))
%   'asinh' |τ| = (a-b)σ asinh( v/(2 v0) exp(f0/(a-b)) )
%
% the two laws are identical when v/(2 v0) exp(f0/(a-b)) >> 1, i.e.,
% unless v is smaller than about v0 exp(-f0/(a-b)). The asinh form stays
% positive and well defined as v -> 0.
%
% SEE ALSO: unicycle.ode.frictionInverse, unicycle.ode.odeRateDependent

switch sys.law
    case 'log'
        tauMag = sys.sigma.*(sys.f0 + sys.ab.*zeta);
        dtaudzeta = sys.abSigma.*ones(size(zeta));
    case 'asinh'
        % x = exp(u)/2 with u = zeta + rho, computed without overflow
        u = zeta + sys.rho;
        tauMag = sys.abSigma.*asinhHalfExp(u);
        % d asinh(x)/du = x/sqrt(1+x^2) = 1/sqrt(1+4 exp(-2u))
        dtaudzeta = sys.abSigma./sqrt(1 + 4*exp(-2*u));
    otherwise
        error('frictionLaw: unknown law %s',sys.law)
end

end

function y = asinhHalfExp(u)
% asinh(exp(u)/2) = u + log(1/2 + sqrt(1/4 + exp(-2u))) for u > 0
y = asinh(exp(u)/2);
pos = u > 0;
y(pos) = u(pos) + log(0.5 + sqrt(0.25 + exp(-2*u(pos))));
end
