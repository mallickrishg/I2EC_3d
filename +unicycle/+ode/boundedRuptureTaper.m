function [weightC,A,dtauPar] = boundedRuptureTaper(sys,Trec,shapeC,dtauBounds,direction,weight)
% BOUNDEDRUPTURETAPER largest coseismic slip on creeping patches, with a
% given taper shape, that keeps the traction change on every creeping patch
% within bounds
%
%   [weightC,A,dtauPar] = unicycle.ode.boundedRuptureTaper(sys,Trec,shapeC,dtauBounds,direction,weight)
%
% The slip on creeping patches inside the rupture is weightC = A*shapeC (in
% units of Vpl Trec, see unicycle.ode.eventStressChange). The traction
% change is linear in A,
%   Δτ(A) = Δτ_P + A Δτ_shape
% and its component along a unit direction d on each creeping patch must
% satisfy
%   dtauBounds(1) <= Δτ(A)·d <= dtauBounds(2)
% A is the largest value in [0,1] for which this holds on every creeping
% patch. Use the pre-event traction direction for d (from
% unicycle.ode.stateTraction), so the lower bound limits how much the
% traction magnitude can drop.
%
% INPUT:
%   sys        - structure from unicycle.ode.assembleODE
%   Trec       - recurrence interval (s)
%   shapeC     - Nc x 1 taper shape in [0,1], zero outside the rupture
%   dtauBounds - [min,max] traction change along d (MPa), min <= 0 <= max
%   direction  - Nc x 2 [strike,dip] direction (normalised here)
%   weight     - (optional) slip weight on pinned patches (default 1)
%
% OUTPUT:
%   weightC - Nc x 1 weight to pass to unicycle.ode.eventStressChange
%   A       - taper amplitude in [0,1]
%   dtauPar - Nc x 1 traction change along d with this taper (MPa)
%
% SEE ALSO: unicycle.ode.eventStressChange, unicycle.ode.stateTraction

if nargin < 6 || isempty(weight), weight = 1; end
Nc = sys.Nc;
assert(numel(shapeC) == Nc,'boundedRuptureTaper: shapeC must be Nc x 1')
assert(all(shapeC >= 0 & shapeC <= 1),'boundedRuptureTaper: shapeC must be in [0,1]')
assert(isequal(size(direction),[Nc,2]),'boundedRuptureTaper: direction must be Nc x 2')
d = direction./hypot(direction(:,1),direction(:,2));

project = @(dtau) d(:,1).*dtau(1:Nc) + d(:,2).*dtau(Nc+1:end);
b = project(unicycle.ode.eventStressChange(sys,Trec,weight,0));
dtauShape = unicycle.ode.eventStressChange(sys,Trec,0,shapeC(:));
c = project(dtauShape);
assert(all(b >= dtauBounds(1) & b <= dtauBounds(2)),...
    ['boundedRuptureTaper: the pinned-patch slip alone gives a traction change ',...
     'in [%.3g, %.3g] MPa, outside the bounds'],min(b),max(b))

% each patch bounds A from above where c pushes Δτ towards a bound
A = 1;
up = c > 0;
A = min([A;(dtauBounds(2) - b(up))./c(up)]);
down = c < 0;
A = min([A;(dtauBounds(1) - b(down))./c(down)]);
A = max(A,0);

weightC = A*shapeC(:);
dtauPar = b + A*c;

end
