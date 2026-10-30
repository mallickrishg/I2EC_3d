function G = computeSurfaceKernels(rcv,xObs)
% COMPUTESURFACEKERNELS displacement at observation points due to unit
% strike slip and dip slip on every patch of the mesh
%
%   G = unicycle.greens.computeSurfaceKernels(rcv,xObs)
%
% INPUT:
%   rcv  - unicycle.geometry.triangleReceiver (or any receiver with
%          displacementKernels)
%   xObs - (D x 3) observation points [east, north, up] (m), the same
%          coordinates as rcv.xc; up = 0 at the surface
%
% OUTPUT:
%   G    - (3D x 2N) [G_strike, G_dip], rows interlaced by station as
%          [e1; n1; u1; e2; n2; u2; ...] (m of displacement per m of slip)
%
% Surface displacement for slip s (N x 2) is G*s(:), and surface velocity
% for slip rate v (N x 2) is G*v(:). See unicycle.ode.surfaceDeformation.
%
% SEE ALSO: unicycle.ode.surfaceDeformation, unicycle.greens.computeKernels

assert(size(xObs,2) == 3,'computeSurfaceKernels: xObs must be D x 3 [east, north, up]')
[Gs,Gd] = rcv.displacementKernels(xObs,3);
G = [Gs,Gd];

end
