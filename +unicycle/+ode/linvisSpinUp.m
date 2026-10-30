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
%            (Nc x 2) and slipP (Np x 2; if given, sum must be Vpl Trec on
%            pinned patches), e.g. unicycle.ode.periodicEvent(lin,Trec,weightC)
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
% SEE ALSO: unicycle.ode.linvisSpinUpSequence (sequences of events),
%           unicycle.ode.linvisAssemble, unicycle.ode.linvisEvolve,
%           unicycle.ode.periodicEvent

% the single-event case of unicycle.ode.linvisSpinUpSequence; the locked
% slip budget is checked when event.slipP is given
assert(isscalar(event),'linvisSpinUp: one event (use linvisSpinUpSequence for a sequence)')
Trec = event.Trec;
if nargin < 3 || isempty(tSince), tSince = Trec; end
assert(tSince >= 0 && tSince <= Trec,'linvisSpinUp: tSince must be in [0, Trec]')
ev = struct('t',0,'dtau',event.dtau,'slipP',[],'slipC',[]);
hasP = isfield(event,'slipP') && ~isempty(event.slipP);
if hasP, ev.slipP = event.slipP; end
if isfield(event,'slipC'), ev.slipC = event.slipC; end
[v,infoS] = unicycle.ode.linvisSpinUpSequence(lin,struct('T',Trec,'events',ev),tSince,...
    struct('checkBudget',hasP));
info.vPre = infoS.vPre;
info.vPost = infoS.vPost;
info.slipCycle = infoS.slipCycle;
info.slipSince = reshape(infoS.slip,[],2);
info.spectralRadius = infoS.spectralRadius;

end
