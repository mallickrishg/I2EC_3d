function [v,info] = linvisSpinUp(lin,event,tSince)
% LINVISSPINUP exact periodic cycle of the linear viscous ODE for a
% characteristic event repeating every Trec
%
%   [v,info] = unicycle.ode.linvisSpinUp(lin,event,tSince)
%
% An event with traction change Δτ makes the slip rate jump by
% b = D⁻¹ Δτ (and the slip on creeping rupture patches by slip_C).
% Between events u = v - v* decays as exp(A t). With E = exp(A Trec) and
% p the pre-event u, periodicity p = E (p + b) gives
%   p = (I - E)⁻¹ E b,   i.e. p_i = b_i e^{λ_i Trec} / (1 - e^{λ_i Trec})
% in modal coordinates, and the post-event state is v* + (I - E)⁻¹ b. Only
% I - E must be invertible (all Re λ < 0). As η -> ∞ (λ -> 0) the solution
% tends to steady creep at Vpl; as η -> 0 it relaxes instantly to v*.
% This is the limit of repeating the event indefinitely (the counterpart of
% unicycle.ode.spinUp), obtained without time stepping. The slip per cycle,
%   Δs = v* Trec - A⁻¹ b + slip_C = -K_cc⁻¹ (Δτ + τ̇_load Trec) + slip_C
% does not depend on η; it equals Vpl Trec when the event releases the
% slip deficit of the locked patches (unicycle.ode.periodicEvent).
%
% INPUT:
%   lin    - structure from unicycle.ode.linvisAssemble
%   event  - structure with Trec, dtau (2Nc x 1) and optionally slipC
%            (Nc x 2), e.g. unicycle.ode.periodicEvent(lin,Trec,weightC)
%   tSince - (optional) time since the last event at which to return the
%            slip rate, 0 <= tSince <= Trec [Trec, just before the next event]
%
% OUTPUT:
%   v    - (2Nc x 1) slip rate at tSince after the last event (m/s)
%   info - structure with
%          vPre, vPost - (2Nc x 1) periodic slip rate just before and just
%                        after the event
%          slipCycle   - (Nc x 2) slip per cycle, coseismic + interseismic
%          slipSince   - (Nc x 2) slip from the event to tSince, including
%                        coseismic slip on creeping patches
%          spectralRadius - max |exp(λ Trec)|, the per-cycle contraction
%                        of the cycle map (< 1)
%
% SEE ALSO: unicycle.ode.linvisAssemble, unicycle.ode.linvisEvolve,
%           unicycle.ode.periodicEvent

Nc = lin.Nc;
Trec = event.Trec;
if nargin < 3 || isempty(tSince), tSince = Trec; end
assert(tSince >= 0 && tSince <= Trec,'linvisSpinUp: tSince must be in [0, Trec]')
slipC = zeros(Nc,2);
if isfield(event,'slipC') && ~isempty(event.slipC)
    slipC = event.slipC;
end

% modal coordinates of the velocity jump and the periodic pre-event state.
% e^{λT} - 1 is evaluated without cancellation: for large η (λT -> 0) both
% b and 1 - e^{λT} vanish like 1/η and their ratio tends to the finite
% limit -K_cc⁻¹Δτ/T (steady creep at Vpl)
bm = lin.Vinv*(event.dtau(:)./lin.d);
z = lin.lambda*Trec;
em1 = unicycle.ode.linvisExpm1(z);
eT = 1 + em1;
pm = -eT./em1.*bm;

info.vPre = lin.vStar + real(lin.V*pm);
info.vPost = info.vPre + event.dtau(:)./lin.d;
[v,s] = unicycle.ode.linvisEvolve(lin,info.vPost,tSince);
[~,sT] = unicycle.ode.linvisEvolve(lin,info.vPost,Trec);
info.slipCycle = [sT(1:Nc),sT(Nc+1:end)] + slipC;
info.slipSince = [s(1:Nc),s(Nc+1:end)] + slipC;
info.spectralRadius = max(abs(eT));

end
