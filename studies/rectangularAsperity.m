function pinned = rectangularAsperity(rcv,xlim,wlim,dip)
% RECTANGULARASPERITY indices of patches whose centroid lies in a rectangle
% in fault-plane coordinates
%
%   pinned = rectangularAsperity(rcv,xlim,wlim,dip)
%
% xlim = [xmin,xmax] along strike (east, m), wlim = [wmin,wmax] down-dip
% distance from the trench (m), dip in degrees. Down-dip distance is
% recovered from the centroid north coordinate as w = north/cos(dip), which
% holds for the planar make_fault_mesh geometry.
%
% Put the rectangle edges on mesh grid lines, so that every grid cell is
% either fully in or fully out and no lone triangles are left unlocked.

xc = rcv.xc;
w = xc(:,2)./cosd(dip);
pinned = find(xc(:,1) > xlim(1) & xc(:,1) < xlim(2) & ...
              w > wlim(1) & w < wlim(2));
end
