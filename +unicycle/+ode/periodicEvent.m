function event = periodicEvent(sys,Trec,weightC)
% PERIODICEVENT imposed earthquake that closes the slip budget of the
% locked region every recurrence interval
%
%   event = unicycle.ode.periodicEvent(sys,Trec,weightC)
%
% The pinned patches slip exactly Vpl Trec along the plate rake in every
% event, so their long-term slip rate is Vpl. Creeping patches inside the
% rupture may slip too (weightC, in units of Vpl Trec, e.g. from
% unicycle.ode.ruptureSlipForStressChange); in a periodic cycle their
% coseismic plus interseismic slip still adds up to Vpl Trec.
%
% INPUT:
%   sys     - structure from unicycle.ode.assembleODE
%   Trec    - recurrence interval (s)
%   weightC - (optional) Nc x 1 coseismic slip on creeping patches in units
%             of Vpl Trec (default 0)
%
% OUTPUT: structure event with
%   Trec        - recurrence interval (s)
%   dtau        - (2Nc x 1) traction change on creeping patches (MPa)
%   slipP,slipC - (Np x 2), (Nc x 2) coseismic slip (m)
%
% SEE ALSO: unicycle.ode.runCycles, unicycle.ode.eventStressChange

if nargin < 3 || isempty(weightC), weightC = 0; end
[dtau,slipP,slipC] = unicycle.ode.eventStressChange(sys,Trec,1,weightC);
event = struct('Trec',Trec,'dtau',dtau,'slipP',slipP,'slipC',slipC);

end
