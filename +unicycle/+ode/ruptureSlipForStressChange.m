function [weightC,dtauPar] = ruptureSlipForStressChange(sys,Trec,inRuptureC,dtauTarget,direction,weight)
% RUPTURESLIPFORSTRESSCHANGE coseismic slip on the creeping patches inside
% a rupture that gives them a prescribed traction change
%
%   [weightC,dtauPar] = unicycle.ode.ruptureSlipForStressChange(sys,Trec,inRuptureC,dtauTarget,direction,weight)
%
% Creeping patches inside the rupture (set E) slip along the plate rake by
% weightC Vpl Trec (see unicycle.ode.eventStressChange). Their weights are
% chosen so that the traction change along a unit direction d is the
% target on every patch of E:
%   (Δτ_P + K_cc slip_E(w))·d = dtauTarget   on E
% which is a square linear system for the weights. Use the pre-event
% traction direction for d (from unicycle.ode.stateTraction), so that
% dtauTarget < 0 is a stress drop of that size.
%
% INPUT:
%   sys        - structure from unicycle.ode.assembleODE
%   Trec       - recurrence interval (s)
%   inRuptureC - Nc x 1 logical, creeping patches inside the rupture
%   dtauTarget - scalar or nE x 1 traction change along d on those patches (MPa)
%   direction  - Nc x 2 [strike,dip] direction (normalised here)
%   weight     - (optional) slip weight on pinned patches (default 1)
%
% OUTPUT:
%   weightC - Nc x 1 weight to pass to unicycle.ode.eventStressChange,
%             zero outside the rupture
%   dtauPar - Nc x 1 traction change along d on all creeping patches (MPa)
%
% SEE ALSO: unicycle.ode.eventStressChange, unicycle.ode.stateTraction

if nargin < 6 || isempty(weight), weight = 1; end
Nc = sys.Nc;
inRuptureC = logical(inRuptureC(:));
assert(numel(inRuptureC) == Nc,'ruptureSlipForStressChange: inRuptureC must be Nc x 1')
assert(isequal(size(direction),[Nc,2]),'ruptureSlipForStressChange: direction must be Nc x 2')
d = direction./hypot(direction(:,1),direction(:,2));
iE = find(inRuptureC);
nE = numel(iE);

% traction change along d from the pinned patches
dtauP = unicycle.ode.eventStressChange(sys,Trec,weight,0);
b = d(:,1).*dtauP(1:Nc) + d(:,2).*dtauP(Nc+1:end);

% traction change along d on all creeping patches per unit weight on E
u = sys.Vpl(sys.iC(iE),:)*Trec;               % slip per unit weight, nE x 2
C = d(:,1).*(sys.Kcc(1:Nc,iE).*u(:,1)' + sys.Kcc(1:Nc,Nc+iE).*u(:,2)') + ...
    d(:,2).*(sys.Kcc(Nc+1:end,iE).*u(:,1)' + sys.Kcc(Nc+1:end,Nc+iE).*u(:,2)');

wE = C(iE,:) \ (dtauTarget(:).*ones(nE,1) - b(iE));
weightC = zeros(Nc,1);
weightC(iE) = wE;
dtauPar = b + C*wE;

end
