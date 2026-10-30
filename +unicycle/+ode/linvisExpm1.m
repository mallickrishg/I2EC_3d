function y = linvisExpm1(z)
% LINVISEXPM1 exp(z) - 1 for real or complex z, without cancellation
%
%   y = unicycle.ode.linvisExpm1(z)
%
% For small |z| (large viscosity, short times) exp(z) - 1 loses digits to
% cancellation; there it is evaluated as 2 exp(z/2) sinh(z/2). For large
% |z| (e.g. Re z -> -Inf for small viscosity) the direct form is exact and
% avoids 0*Inf in the sinh form.
%
% SEE ALSO: unicycle.ode.linvisEvolve, unicycle.ode.linvisSpinUp

y = exp(z) - 1;
small = abs(z) < 1;
y(small) = 2*exp(z(small)/2).*sinh(z(small)/2);

end
