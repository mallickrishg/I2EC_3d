function [vSurf,uSurf] = surfaceDeformation(G,sys,t,Y,opts)
% SURFACEDEFORMATION surface velocity and displacement time series at
% observation points from an ODE solution
%
%   [vSurf,uSurf] = unicycle.ode.surfaceDeformation(G,sys,t,Y,opts)
%
% Surface velocity is the kernel times the on-fault slip rate,
%   v_surf(t) = G_C v_C(t)          (pinned patches have v = 0)
% and surface displacement is its time integral. Because the ODE state
% already holds the slip (the time integral of the slip rate, integrated
% by ode45 with the rest of the state), the displacement is evaluated
% exactly from it rather than by quadrature:
%   u_surf(t) = G_C slip_C(t) + G_P slip_P
% where slip_C is the slip state and slip_P the slip of the pinned patches
% (e.g. the coseismic slip of the event that starts the time series).
%
% With opts.removeLongTerm, the deformation due to steady slip at Vpl on
% the whole mesh is removed,
%   v_surf(t) - G Vpl,     u_surf(t) - G Vpl t
% which leaves the signal relative to the long-term (plate) motion, i.e.
% the effect of the slip deficit (G (v - Vpl), the back-slip view).
%
% INPUT:
%   G    - (3D x 2N) kernel from unicycle.greens.computeSurfaceKernels
%   sys  - structure from unicycle.ode.assembleODE
%   t    - (nt x 1) times (s)
%   Y    - (nt x 4Nc) states [slip_s, slip_d, zeta, theta] at t, e.g. from
%          ode45 or unicycle.ode.runCycles (out.t, out.y)
%   opts - (optional) structure with fields (defaults in brackets)
%          removeLongTerm - remove G Vpl from the velocity and G Vpl t
%                           from the displacement [false]
%          pinnedSlip     - (Np x 2) slip of the pinned patches, held
%                           fixed over the time series [zeros]. For a cycle
%                           from runCycles pass out.event.slipP so that the
%                           displacement includes the coseismic step.
%          reference      - 'zero' (displacement relative to zero slip
%                           state and zero pinnedSlip; for runCycles output
%                           that is just before the event) or 'first'
%                           (relative to the first time) ['zero']
%
% OUTPUT:
%   vSurf - (nt x D x 3) [east, north, up] surface velocity (m/s)
%   uSurf - (nt x D x 3) [east, north, up] surface displacement (m)
%
% SEE ALSO: unicycle.greens.computeSurfaceKernels, unicycle.ode.runCycles

if nargin < 5, opts = struct(); end
if ~isfield(opts,'removeLongTerm') || isempty(opts.removeLongTerm), opts.removeLongTerm = false; end
if ~isfield(opts,'pinnedSlip') || isempty(opts.pinnedSlip), opts.pinnedSlip = zeros(numel(sys.iP),2); end
if ~isfield(opts,'reference') || isempty(opts.reference), opts.reference = 'zero'; end
assert(any(strcmp(opts.reference,{'zero','first'})),'surfaceDeformation: reference must be ''zero'' or ''first''')

N = sys.N;
Nc = sys.Nc;
iC = sys.iC;
iP = sys.iP;
assert(size(G,2) == 2*N,'surfaceDeformation: G must be 3D x 2N')
assert(mod(size(G,1),3) == 0,'surfaceDeformation: G must have 3 rows per station')
assert(size(Y,2) == 4*Nc && size(Y,1) == numel(t),'surfaceDeformation: Y must be nt x 4Nc')
assert(isequal(size(opts.pinnedSlip),[numel(iP),2]),'surfaceDeformation: pinnedSlip must be Np x 2')
D = size(G,1)/3;
nt = numel(t);
t = t(:);

GC = G(:,[iC;N+iC]);
GP = G(:,[iP;N+iP]);

% on-fault slip rate and slip on creeping patches (2Nc x nt)
vMag = sys.v0.*exp(Y(:,2*Nc+1:3*Nc)');
theta = Y(:,3*Nc+1:end)';
vC = [vMag.*cos(theta);vMag.*sin(theta)];
sC = Y(:,1:2*Nc)';

v = GC*vC;
u = GC*sC + GP*opts.pinnedSlip(:);
if opts.removeLongTerm
    vLong = G*sys.Vpl(:);
    v = v - vLong;
    u = u - vLong*t';
end
if strcmp(opts.reference,'first')
    u = u - u(:,1);
end

% (3D x nt), rows [e1 n1 u1 e2 ...] -> (nt x D x 3)
vSurf = permute(reshape(v,3,D,nt),[3,2,1]);
uSurf = permute(reshape(u,3,D,nt),[3,2,1]);

end
