function y = applyEvent(sys,y,dtau,slipC)
% APPLYEVENT coseismic jump of the ODE state for an imposed traction change
%
%   y = unicycle.ode.applyEvent(sys,y,dtau,slipC)
%
% The jump is exact in traction space. With T⁻ the traction before the
% event (see unicycle.ode.stateTraction),
%   T⁺ = T⁻ + Δτ
%   θ⁺ = atan2(T⁺_d, T⁺_s)
%   F(v⁺) + η v⁺ = |T⁺|
% so the postseismic velocity is limited by radiation damping (solved with
% unicycle.ode.frictionInverse). Slip on creeping patches jumps by the
% imposed coseismic slip slipC (zero unless the rupture extends onto
% creeping patches).
%
% A large stress drop can reverse the traction (T⁺·T⁻ < 0). The asinh law
% handles this (v ∥ T⁺, backward slip); the log law has no meaning for
% small |T⁺|.
%
% INPUT:
%   sys  - structure from unicycle.ode.assembleODE
%   y    - state vector before the event [slip_s; slip_d; zeta; theta]
%   dtau - (2Nc x 1) [strike; dip] traction change on creeping patches,
%          e.g. from unicycle.ode.eventStressChange
%   slipC - (optional) Nc x 2 [strike,dip] coseismic slip on creeping
%          patches (default 0)
%
% OUTPUT:
%   y    - state vector after the event
%
% SEE ALSO: unicycle.ode.eventStressChange, unicycle.ode.stateTraction

Nc = sys.Nc;
assert(numel(dtau) == 2*Nc,'applyEvent: dtau must be 2Nc x 1')
if nargin > 3 && ~isempty(slipC)
    assert(isequal(size(slipC),[Nc,2]),'applyEvent: slipC must be Nc x 2')
    y(1:2*Nc) = y(1:2*Nc) + slipC(:);
end

T = unicycle.ode.stateTraction(sys,y) + [dtau(1:Nc),dtau(Nc+1:end)];
tauMag = hypot(T(:,1),T(:,2));
% a traction that cancels exactly would give v = 0 and zeta = -Inf
tauMag = max(tauMag,eps*sys.abSigma);

y(2*Nc+1:3*Nc) = unicycle.ode.frictionInverse(sys,tauMag);
y(3*Nc+1:4*Nc) = atan2(T(:,2),T(:,1));

end
