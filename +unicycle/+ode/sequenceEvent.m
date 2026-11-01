function event = sequenceEvent(sys,t,slip)
% SEQUENCEEVENT one imposed earthquake of an event sequence, from its slip
% on any patches (locked, creeping or both)
%
%   event = unicycle.ode.sequenceEvent(sys,t,slip)
%
% The slip is split into the pinned and creeping parts and the traction
% change on creeping patches is
%   Δτ_c = K_cp slip_P + K_cc slip_C
% as in unicycle.ode.eventStressChange. Collect the events in a struct
% array, seq = struct('T',T,'events',events), for unicycle.ode.spinUp,
% unicycle.ode.runCycles or unicycle.ode.linvisSpinUp. There need not be
% any pinned patches.
%
% INPUT:
%   sys  - structure from unicycle.ode.assembleODE or linvisAssemble
%   t    - time of the event within the supercycle (s)
%   slip - (N x 2) [strike,dip] coseismic slip on all patches (m), zero
%          outside the rupture
%
% OUTPUT: structure event with
%   t           - event time (s)
%   dtau        - (2Nc x 1) traction change on creeping patches (MPa)
%   slipP,slipC - (Np x 2), (Nc x 2) coseismic slip (m)
%
% SEE ALSO: unicycle.ode.spinUp, unicycle.ode.runCycles,
%           unicycle.ode.linvisSpinUp, unicycle.ode.eventStressChange

N = numel(sys.iC) + numel(sys.iP);
assert(isequal(size(slip),[N,2]),'sequenceEvent: slip must be N x 2')
slipP = slip(sys.iP,:);
slipC = slip(sys.iC,:);
dtau = sys.Kcp*slipP(:) + sys.Kcc*slipC(:);
event = struct('t',t,'dtau',dtau,'slipP',slipP,'slipC',slipC);

end
