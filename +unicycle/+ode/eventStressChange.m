function [dtau,slipP,slipC] = eventStressChange(sys,Trec,weight,weightC)
% EVENTSTRESSCHANGE coseismic slip of an imposed earthquake and the
% traction change it imposes on the creeping patches
%
%   [dtau,slipP,slipC] = unicycle.ode.eventStressChange(sys,Trec,weight,weightC)
%
% Pinned patches are held at v = 0 between events, so over a recurrence
% interval Trec they accumulate a slip deficit Vpl Trec along the plate
% rake. The imposed earthquake releases it:
%   slip_P = weight .* Vpl_P Trec
% The rupture may also extend onto creeping patches (weightC > 0), which
% then slip coseismically too:
%   slip_C = weightC .* Vpl_C Trec
% and the traction change on creeping patches is
%   Δτ_c = K_cp slip_P + K_cc slip_C
% Creeping patches inside the rupture get a stress drop from their own
% slip (K_cc has a negative diagonal).
%
% INPUT:
%   sys     - structure from unicycle.ode.assembleODE
%   Trec    - recurrence interval (s)
%   weight  - (optional) Np x 1 or scalar fraction of the deficit released
%             on each pinned patch, e.g. an edge taper (default 1). With
%             weight < 1 the deficit is not fully released and the
%             long-term slip rate of the asperity is weight*Vpl
%   weightC - (optional) Nc x 1 or scalar coseismic slip on creeping
%             patches in units of Vpl Trec (default 0)
%
% OUTPUT:
%   dtau  - (2Nc x 1) [strike; dip] traction change on creeping patches (MPa)
%   slipP - (Np x 2) [strike,dip] coseismic slip on pinned patches (m)
%   slipC - (Nc x 2) [strike,dip] coseismic slip on creeping patches (m),
%           pass it to unicycle.ode.applyEvent
%
% SEE ALSO: unicycle.ode.applyEvent, unicycle.ode.assembleODE

if nargin < 3 || isempty(weight), weight = 1; end
if nargin < 4 || isempty(weightC), weightC = 0; end
Np = numel(sys.iP);
assert(Np > 0,'eventStressChange: there are no pinned patches')
assert(isscalar(weight) || numel(weight) == Np,'eventStressChange: weight must be scalar or Np x 1')
assert(isscalar(weightC) || numel(weightC) == sys.Nc,'eventStressChange: weightC must be scalar or Nc x 1')

slipP = weight(:).*sys.Vpl(sys.iP,:)*Trec;
slipC = weightC(:).*sys.Vpl(sys.iC,:)*Trec;
dtau = sys.Kcp*slipP(:);
if any(slipC(:))
    dtau = dtau + sys.Kcc*slipC(:);
end

end
