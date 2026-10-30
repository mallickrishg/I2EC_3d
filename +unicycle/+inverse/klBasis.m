function basis = klBasis(rcv,iC,opts)
% KLBASIS Karhunen-Loève basis of a Matérn covariance on fault patches,
% for low-rank fields such as log effective viscosity
%
%   basis = unicycle.inverse.klBasis(rcv,iC,opts)
%
% A field m on the patches is written as
%   m = m̄ + Phi c,   Phi = phi diag(sqrt(lambda)),   c ~ N(0, I)
% where (lambda, phi) are eigenpairs of the covariance operator
%   C_ij = sigma² Matérn_nu(r_ij / ell)
% between patch centroids. The eigenproblem is area weighted,
%   C W phi = lambda phi,   phi' W phi = I,   W = diag(area)
% (the discrete version of the integral operator), so the modes and
% eigenvalues do not depend on the mesh resolution. The mean level m̄ is
% a separate unknown: by default the covariance is centred on the
% area-weighted mean (removeMean), so all modes have zero area-weighted
% mean and only describe the shape of the field. The modes are nested:
% using the first r of them for r = 0, 1, 2, ... gives nested models.
%
% INPUT:
%   rcv  - unicycle.geometry.triangleReceiver (uses xc and area)
%   iC   - (optional) indices of the patches to build the basis on, e.g.
%          lin.iC or sys.iC [all patches]
%   opts - (optional) structure with fields (defaults in brackets)
%          ell        - correlation length (m) [L/2, with L the largest
%                       principal extent of the centroids]
%          nu         - Matérn smoothness: 0.5 (exponential), 1.5, 2.5,
%                       Inf (Gaussian) or any positive value [0.5]
%                       (general values use the Bessel form in log space)
%          sigma      - prior standard deviation of the field [1]
%          removeMean - centre the covariance on the area-weighted mean [true]
%          rMax       - number of modes returned [20]
%
% OUTPUT: structure basis with
%   Phi     - (n x rMax) scaled modes, m = m̄ + Phi(:,1:r) c(1:r)
%   phi     - (n x rMax) modes, orthonormal in the area-weighted inner product
%   lambda  - (rMax x 1) eigenvalues, decreasing
%   varFraction - (rMax x 1) cumulative fraction of the prior variance
%             captured by the first r modes
%   ell, nu, sigma, removeMean, L, iC, xc, area
%
% SEE ALSO: unicycle.ode.linvisAssemble

if nargin < 2 || isempty(iC), iC = (1:rcv.N)'; end
if nargin < 3, opts = struct(); end
iC = iC(:);
xc = rcv.xc(iC,:);
area = rcv.area(iC);
area = area(:);
n = numel(iC);

% largest principal extent of the centroids
xm = xc - mean(xc,1);
[~,~,Vp] = svd(xm,'econ');
proj = xm*Vp(:,1);
L = max(proj) - min(proj);

def = struct('ell',L/2,'nu',0.5,'sigma',1,'removeMean',true,'rMax',20);
names = fieldnames(def);
for k = 1:numel(names)
    if ~isfield(opts,names{k}) || isempty(opts.(names{k}))
        opts.(names{k}) = def.(names{k});
    end
end
assert(opts.ell > 0 && opts.nu > 0 && opts.sigma > 0,'klBasis: ell, nu and sigma must be positive')
rMax = min(opts.rMax,n - opts.removeMean);

% covariance between centroids
r = sqrt(max(sum(xc.^2,2) + sum(xc.^2,2)' - 2*(xc*xc'),0));
C = opts.sigma^2*matern(r/opts.ell,opts.nu);
C = (C + C')/2;

% area-weighted centring: C <- P C P', P = I - 1 w', w = area/sum(area)
if opts.removeMean
    w = area/sum(area);
    P = eye(n) - ones(n,1)*w';
    C = P*C*P';
end

% symmetric form of C W phi = lambda phi: (W^½ C W^½) psi = lambda psi
sw = sqrt(area);
B = sw.*C.*sw';
B = (B + B')/2;
[Psi,Lam] = eig(B);
[lam,order] = sort(diag(Lam),'descend');
Psi = Psi(:,order);
lam = max(lam,0);
phi = Psi./sw;

% fix the sign of each mode (largest-magnitude entry positive) so the
% basis is reproducible
[~,imax] = max(abs(phi),[],1);
sgn = sign(phi(sub2ind(size(phi),imax,1:n)));
phi = phi.*sgn;

basis.phi = phi(:,1:rMax);
basis.lambda = lam(1:rMax);
basis.Phi = basis.phi.*sqrt(basis.lambda)';
basis.varFraction = cumsum(basis.lambda)/sum(lam);
basis.ell = opts.ell;
basis.nu = opts.nu;
basis.sigma = opts.sigma;
basis.removeMean = opts.removeMean;
basis.L = L;
basis.iC = iC;
basis.xc = xc;
basis.area = area;

end

function c = matern(x,nu)
% Matérn correlation as a function of x = r/ell
switch nu
    case 0.5
        c = exp(-x);
    case 1.5
        s = sqrt(3)*x;
        c = (1 + s).*exp(-s);
    case 2.5
        s = sqrt(5)*x;
        c = (1 + s + s.^2/3).*exp(-s);
    case Inf
        c = exp(-x.^2/2);
    otherwise
        % log form with the scaled Bessel function to avoid overflow of
        % gamma(nu) and besselk for large nu; where besselk still overflows
        % (very small s) use the small-distance expansion
        s = sqrt(2*nu)*x;
        logc = (1 - nu)*log(2) - gammaln(nu) + nu*log(s) + log(besselk(nu,s,1)) - s;
        c = exp(logc);
        bad = ~isfinite(logc);
        if nu > 1
            c(bad) = 1 - s(bad).^2/(4*(nu - 1));
        end
        c(x == 0) = 1;
end
end
