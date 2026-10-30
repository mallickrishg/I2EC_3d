function v = linvisEvent(lin,v,dtau)
% LINVISEVENT slip-rate jump of the linear viscous model for an imposed
% traction change
%
%   v = unicycle.ode.linvisEvent(lin,v,dtau)
%
% With τ = τ0 + η v the response to a traction step Δτ is instantaneous,
%   v⁺ = v⁻ + Δτ / η
% (no radiation-damping cap is needed, the law is linear). Slip on
% creeping patches inside the rupture, if any, jumps separately by slip_C.
%
% INPUT:
%   lin  - structure from unicycle.ode.linvisAssemble
%   v    - (2Nc x 1) slip rate before the event (m/s)
%   dtau - (2Nc x 1) traction change, e.g. from unicycle.ode.eventStressChange
%
% SEE ALSO: unicycle.ode.linvisEvolve, unicycle.ode.linvisSpinUp

assert(numel(dtau) == 2*lin.Nc && numel(v) == 2*lin.Nc,'linvisEvent: v and dtau must be 2Nc x 1')
v = v(:) + dtau(:)./lin.d;

end
