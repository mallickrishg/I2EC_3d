function sys = oneSliderSystem(k,Vpl,rakePl,law,f0,eta)
% ONESLIDERSYSTEM a single creeping patch loaded by an isotropic spring,
% τ̇ = -k (v - Vpl), in the format returned by unicycle.ode.assembleODE
%
%   sys = oneSliderSystem(k,Vpl,rakePl,law,f0,eta)
%
% k [MPa/m], Vpl [m/s], rakePl [rad]

if nargin < 4, law = 'log'; end
if nargin < 5, f0 = 0.6; end
if nargin < 6, eta = 0; end

sys.N = 1;
sys.Nc = 1;
sys.iC = 1;
sys.law = law;
sys.Vpl = Vpl*[cos(rakePl),sin(rakePl)];
sys.Kcc = -k*eye(2);
sys.taudotLoad = k*sys.Vpl(:);
sys.sigma = 50;
sys.ab = 0.005;
sys.abSigma = sys.ab*sys.sigma;
sys.f0 = f0;
sys.rho = f0/sys.ab;
sys.v0 = 1e-6;
sys.eta = eta;
sys.y0 = [0;0;log(Vpl/sys.v0);rakePl];

end
