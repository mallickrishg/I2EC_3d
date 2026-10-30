function [v,info] = linvisSpinUpSequence(lin,seq,tOut,opts)
% LINVISSPINUPSEQUENCE exact periodic cycle of the linear viscous ODE for a
% sequence of events that repeats every supercycle T
%
%   [v,info] = unicycle.ode.linvisSpinUpSequence(lin,seq,tOut,opts)
%
% The slip rate obeys the inhomogeneous linear ODE
%   v̇ = A v + D⁻¹ τ̇_load,   A = D⁻¹ K_cc
% and event k at time t_k (0 <= t_k < T) makes it jump by b_k = D⁻¹ Δτ_k.
% The constant forcing is absorbed by the steady state v* = -K_cc⁻¹ τ̇_load:
% u = v - v* obeys u̇ = A u, so u(t) = E(t) u(0) with E(t) = exp(A t), and
% the jumps of u equal those of v. With p = u(0⁻), one supercycle gives
%   u(T⁻) = E(T) p + Σ_k E(T - t_k) b_k
% and periodicity u(T⁻) = p gives
%   p = (I - E(T))⁻¹ Σ_k E(T - t_k) b_k
% in modal coordinates p_i = Σ_k b_ik e^{λ_i (T - t_k)} / (1 - e^{λ_i T}).
% Equivalently u(t) = Σ_k G_T(t - t_k) b_k with the periodic Green's
% function G_T(s) = E(s mod T) (I - E(T))⁻¹: the periodic cycle of a
% sequence is the sum of single-event periodic cycles of period T, shifted
% to the event times. Slip since t = 0⁻ is
%   s(t) = v* t + ∫₀ᵗ u + Σ_{t_k <= t} slip_C,k
% For a single event at t = 0 this is unicycle.ode.linvisSpinUp.
%
% Slip budget. Over a supercycle the slip on creeping patches is
%   Δs_C = -K_cc⁻¹ (Σ_k Δτ_k + τ̇_load T) + Σ_k slip_C,k
% independent of η. With Δτ_k = K_cp slip_P,k + K_cc slip_C,k it equals
%   Vpl_C T - K_cc⁻¹ K_cp (Σ_k slip_P,k - Vpl_P T)
% so it is Vpl_C T on every creeping patch exactly when the events release
% the slip deficit of the locked patches, Σ_k slip_P,k = Vpl_P T, wherever
% else they rupture. That condition is checked (an error if it fails, see
% opts); the creeping budget is reported in info.budgetErrC.
%
% INPUT:
%   lin  - structure from unicycle.ode.linvisAssemble
%   seq  - structure with
%          T      - supercycle period (s)
%          events - struct array, one per event, with fields
%                   t     - event time, 0 <= t < T (s)
%                   slipP - (Np x 2) coseismic slip on pinned patches (m)
%                   slipC - (optional) (Nc x 2) coseismic slip on creeping
%                           patches (m)
%                   dtau  - (optional) (2Nc x 1) traction change on creeping
%                           patches (MPa) [K_cp slip_P + K_cc slip_C]
%                   e.g. from unicycle.ode.sequenceEvent
%   tOut - (optional) times in [0, T] at which to return the slip rate [T].
%          An event at t_k is included at t = t_k (post-event state); t = T
%          is the state just before the first event of the next supercycle
%   opts - (optional) structure with
%          checkBudget - error if the locked slip budget does not close [true]
%          tolBudget   - tolerance on max |Σ slip_P - Vpl_P T| / (|Vpl| T) [1e-6]
%
% OUTPUT:
%   v    - (2Nc x nt) slip rate at tOut (m/s)
%   info - structure with
%          t           - (M x 1) event times, sorted; events are reordered
%                        by time in the outputs below (order = sort index)
%          vPre, vPost - (2Nc x M) periodic slip rate just before and just
%                        after each event
%          slip        - (2Nc x nt) slip since t = 0⁻ at tOut, including
%                        coseismic slip on creeping patches
%          slipCycle   - (Nc x 2) slip per supercycle on creeping patches
%          slipCycleP  - (Np x 2) coseismic slip per supercycle on pinned
%                        patches
%          budgetErrP  - max |Σ slip_P - Vpl_P T| / (|Vpl| T) (NaN if slipP
%                        is not given and the check is off)
%          budgetErrC  - max |Δs_C - Vpl_C T| / (|Vpl| T)
%          spectralRadius - max |exp(λ T)|, the per-supercycle contraction
%
% SEE ALSO: unicycle.ode.linvisSpinUp, unicycle.ode.sequenceEvent,
%           unicycle.ode.linvisEvolve, unicycle.ode.linvisSpinUpNumerical

if nargin < 4, opts = struct(); end
def = struct('checkBudget',true,'tolBudget',1e-6);
names = fieldnames(def);
for k = 1:numel(names)
    if ~isfield(opts,names{k}) || isempty(opts.(names{k}))
        opts.(names{k}) = def.(names{k});
    end
end
Nc = lin.Nc;
Np = numel(lin.iP);
n = 2*Nc;
T = seq.T;
assert(isscalar(T) && T > 0,'linvisSpinUpSequence: seq.T must be a positive scalar')
ev = seq.events(:);
M = numel(ev);
assert(M >= 1,'linvisSpinUpSequence: seq.events is empty')
tk = [ev.t]';
assert(numel(tk) == M && all(tk >= 0 & tk < T),'linvisSpinUpSequence: event times must be in [0, T)')
[tk,order] = sort(tk);
ev = ev(order);
if nargin < 3 || isempty(tOut), tOut = T; end
tOut = tOut(:)';
assert(all(tOut >= 0 & tOut <= T),'linvisSpinUpSequence: tOut must be in [0, T]')

% coseismic slip and traction change of each event
dtau = zeros(n,M);
slipC = zeros(n,M);
slipP = zeros(2*Np,M);
hasP = true;
for k = 1:M
    if isfield(ev,'slipC') && ~isempty(ev(k).slipC)
        assert(numel(ev(k).slipC) == n,'linvisSpinUpSequence: slipC must be Nc x 2')
        slipC(:,k) = ev(k).slipC(:);
    end
    if isfield(ev,'slipP') && ~isempty(ev(k).slipP)
        assert(numel(ev(k).slipP) == 2*Np,'linvisSpinUpSequence: slipP must be Np x 2')
        slipP(:,k) = ev(k).slipP(:);
    elseif Np > 0
        hasP = false;
    end
    if isfield(ev,'dtau') && ~isempty(ev(k).dtau)
        assert(numel(ev(k).dtau) == n,'linvisSpinUpSequence: dtau must be 2Nc x 1')
        dtau(:,k) = ev(k).dtau(:);
    else
        assert(hasP,'linvisSpinUpSequence: event %d has neither dtau nor slipP',order(k))
        dtau(:,k) = lin.Kcp*slipP(:,k) + lin.Kcc*slipC(:,k);
    end
end

% slip budget of the locked patches
VplMag = max(hypot(lin.Vpl(:,1),lin.Vpl(:,2)));
info.slipCycleP = reshape(sum(slipP,2),Np,2);
info.budgetErrP = NaN;
if hasP
    mis = info.slipCycleP - lin.Vpl(lin.iP,:)*T;
    info.budgetErrP = max([0;hypot(mis(:,1),mis(:,2))])/(VplMag*T);
end
if opts.checkBudget
    assert(hasP,'linvisSpinUpSequence: the budget check needs slipP for every event (or opts.checkBudget = false)')
    if info.budgetErrP > opts.tolBudget
        error('unicycle:ode:linvisSpinUpSequence:budget',...
            'linvisSpinUpSequence: locked slip per supercycle differs from Vpl T by up to %.2e |Vpl| T (tolerance %.1e)',...
            info.budgetErrP,opts.tolBudget)
    end
end

% periodic pre-cycle state p = u(0⁻) in modal coordinates. e^{λT} - 1 is
% evaluated without cancellation: for large η both b and e^{λT} - 1 vanish
% like 1/η and p stays finite
lam = lin.lambda;
b = lin.Vinv*(dtau./lin.d);
em1 = unicycle.ode.linvisExpm1(lam*T);
p = -sum(exp(lam.*(T - tk')).*b,2)./em1;

[u,su,on] = modalState(lam,p,b,tk,tOut,true);
v = lin.vStar + real(lin.V*u);
info.t = tk;
info.order = order;
info.slip = lin.vStar.*tOut + real(lin.V*su) + slipC*on;
uPre = modalState(lam,p,b,tk,tk',false);
info.vPre = lin.vStar + real(lin.V*uPre);
info.vPost = info.vPre + dtau./lin.d;
[~,suT] = modalState(lam,p,b,tk,T,true);
sT = lin.vStar*T + real(lin.V*suT) + sum(slipC,2);
info.slipCycle = reshape(sT,Nc,2);
mis = info.slipCycle - lin.Vpl(lin.iC,:)*T;
info.budgetErrC = max(hypot(mis(:,1),mis(:,2)))/(VplMag*T);
info.spectralRadius = max(abs(1 + em1));

end

function [u,su,on] = modalState(lam,p,b,tk,t,inclusive)
% modal u(t) and ∫₀ᵗ u from u(0⁻) = p, adding the jumps b_k of events with
% t_k < t (t_k <= t if inclusive). on (M x nt) marks the included events
z = lam.*t;
u = exp(z).*p;
su = unicycle.ode.linvisExpm1(z)./lam.*p;
on = false(numel(tk),numel(t));
for k = 1:numel(tk)
    dt = t - tk(k);
    on(k,:) = dt > 0 | (inclusive & dt == 0);
    dt(~on(k,:)) = 0;
    z = lam.*dt;
    u = u + on(k,:).*exp(z).*b(:,k);
    su = su + on(k,:).*unicycle.ode.linvisExpm1(z)./lam.*b(:,k);
end
end
