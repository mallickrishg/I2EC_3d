function tests = testRunCycles
% imposed earthquake cycles: unicycle.ode.periodicEvent,
% unicycle.ode.spinUp and unicycle.ode.runCycles, for one event, a
% periodic sequence of events and one pass through a sequence
tests = functiontests(localfunctions);
end

function setupOnce(testCase)
addpath(fileparts(fileparts(mfilename('fullpath'))));
testCase.TestData.yr = 3.1536e7;
testCase.TestData.dir = tempname;
mkdir(testCase.TestData.dir);
testCase.TestData.base = makeTestMesh(testCase.TestData.dir,'plane');
rcv = unicycle.geometry.triangleReceiver(testCase.TestData.base,unicycle.greens.nikkhoo15(30e3,0.25));
testCase.TestData.bigK = unicycle.greens.computeKernels(rcv,'plane',fullfile(testCase.TestData.dir,'kernels'));
end

function teardownOnce(testCase)
rmdir(testCase.TestData.dir,'s');
end

function [sys,rcv] = meshSystem(testCase,f0,pinned)
% pinned = false: no locked zone, every patch creeps [true]
if nargin < 3, pinned = true; end
rcv = unicycle.geometry.triangleReceiver(testCase.TestData.base,unicycle.greens.nikkhoo15(30e3,0.25));
rcv = setTestProperties(rcv);
rcv.mu0 = f0;
rcv.Vo = rcv.Vpl;
if pinned
    rcv.pinnedPosition = find(hypot(rcv.xc(:,1),rcv.xc(:,2)-14e3) < 11e3);
end
sys = unicycle.ode.assembleODE(rcv,testCase.TestData.bigK,'asinh');
end

function seq = squareSequence(sys,rcv,tYr,slipMag,Tyr)
% uniform slip along the plate rake on squares of patches (no locked zone)
yr = 3.1536e7;
centres = [-15e3,10e3;5e3,15e3;15e3,8e3];
for k = 1:numel(tYr)
    c = centres(mod(k-1,3)+1,:);
    in = abs(rcv.xc(:,1) - c(1)) < 12e3 & abs(rcv.xc(:,2) - c(2)) < 12e3;
    slip = zeros(rcv.N,2);
    slip(in,:) = slipMag(k)*repmat([cosd(rcv.Vrake),sind(rcv.Vrake)],nnz(in),1);
    events(k) = unicycle.ode.sequenceEvent(sys,tYr(k)*yr,slip); %#ok<AGROW>
end
seq = struct('T',Tyr*yr,'events',events);
end

function testPeriodicEventClosesLockedBudget(testCase)
% with Vpl Trec released on the pinned patches, the periodic slip per
% cycle is Vpl Trec on every creeping patch, with or without coseismic
% slip on creeping patches
sys = meshSystem(testCase,0.6);
Nc = sys.Nc;
Trec = 100*testCase.TestData.yr;
for wC = {0,linspace(0,0.5,Nc)'}
    event = unicycle.ode.periodicEvent(sys,Trec,wC{1});
    verifyEqual(testCase,event.slipP,sys.Vpl(sys.iP,:)*Trec,'RelTol',1e-12)
    s = event.slipC(:) - sys.Kcc \ (event.dtau + sys.taudotLoad*Trec);
    VplC = sys.Vpl(sys.iC,:)*Trec;
    verifyEqual(testCase,s,VplC(:),'AbsTol',1e-9*max(abs(VplC(:))))
end
end

function testOneSliderCycleConverges(testCase)
% spring slider with a stress step Δτ every cycle: the periodic slip per
% cycle is Vpl Trec + Δτ/k along the loading direction
yr = testCase.TestData.yr;
k = 0.1; Vpl = 1e-9; rake = deg2rad(30); dtau = 2;
sys = oneSliderSystem(k,Vpl,rake,'asinh',0.6,30e3/(2*3e3));
sys.v0 = Vpl;
sys.y0(3) = 0;
Trec = 100*yr;
event = struct('Trec',Trec,'dtau',dtau*[cos(rake);sin(rake)],'slipP',[],'slipC',[]);
out = unicycle.ode.runCycles(sys,event,struct('maxCycles',200,'verbose',0,'y0',sys.y0));
verifyTrue(testCase,out.converged)
verifyEqual(testCase,out.slipTarget,(Vpl*Trec + dtau/k)*[cos(rake),sin(rake)],'RelTol',1e-12)
verifyEqual(testCase,out.slipCycle,out.slipTarget,'AbsTol',2e-3*Vpl*Trec)
verifyEqual(testCase,out.y(:,4),rake*ones(numel(out.t),1),'AbsTol',1e-8)
end

function testMeshCyclesCloseSlipBudget(testCase)
% small mesh, asperity releasing Vpl Trec every 50 yr, f0 = 0.1: the
% spin-up converges to a cycle in which every creeping patch slips Vpl Trec
yr = testCase.TestData.yr;
sys = meshSystem(testCase,0.1);
Trec = 50*yr;
event = unicycle.ode.periodicEvent(sys,Trec);
out = unicycle.ode.runCycles(sys,event,struct('maxCycles',500,'verbose',0));
verifyTrue(testCase,out.converged)
VplC = sys.Vpl(sys.iC,:)*Trec;
verifyEqual(testCase,out.slipTarget,VplC,'AbsTol',1e-9*max(abs(VplC(:))))
verifyEqual(testCase,out.slipCycle,VplC,'AbsTol',2e-3*max(abs(VplC(:))))
verifyEqual(testCase,out.history.nReversed,zeros(out.nCycles,1))
% the stored cycle ends where it started
Nc = sys.Nc;
yEnd = out.y(end,:)';
y0 = out.yStart(:,end);
verifyLessThan(testCase,max(abs(yEnd(2*Nc+1:3*Nc) - y0(2*Nc+1:3*Nc))),2e-3)
verifyLessThan(testCase,max(abs(angle(exp(1i*(yEnd(3*Nc+1:end) - y0(3*Nc+1:end)))))),deg2rad(2e-3))
end

function testLockedBudgetWarning(testCase)
sys = meshSystem(testCase,0.6);
event = unicycle.ode.periodicEvent(sys,100*testCase.TestData.yr);
event.slipP = 0.9*event.slipP;
verifyWarning(testCase,@() unicycle.ode.runCycles(sys,event,struct('maxCycles',1,'verbose',0)),...
    'unicycle:ode:spinUp:lockedBudget')
end

function testSpinUpMatchesRunCycles(testCase)
% spinUp returns the periodic pre-event state, or the state a time tSince
% into the periodic cycle, consistent with the cycle stored by runCycles
yr = testCase.TestData.yr;
sys = meshSystem(testCase,0.1);
Nc = sys.Nc;
Trec = 50*yr;
event = unicycle.ode.periodicEvent(sys,Trec);
opts = struct('maxCycles',500,'verbose',0);
out = unicycle.ode.runCycles(sys,event,opts);
[y,info] = unicycle.ode.spinUp(sys,event,opts);
verifyTrue(testCase,info.converged)
verifyEqual(testCase,info.nCycles,out.nCycles)
verifyEqual(testCase,y,out.yStart(:,end),'AbsTol',1e-12)
verifyEqual(testCase,y(1:2*Nc),zeros(2*Nc,1))
% part way through the cycle
k = find(out.t > 3*yr,1);
opts.tSince = out.t(k);
[y,info] = unicycle.ode.spinUp(sys,event,opts);
verifyEqual(testCase,y(2*Nc+1:end),out.y(k,2*Nc+1:end)','AbsTol',1e-6)
verifyEqual(testCase,info.slipSince(:),out.y(k,1:2*Nc)','AbsTol',1e-6*max(abs(out.y(k,1:2*Nc))))
verifyEqual(testCase,y(1:2*Nc),zeros(2*Nc,1))
% right after the event
opts.tSince = 0;
y = unicycle.ode.spinUp(sys,event,opts);
verifyEqual(testCase,y(2*Nc+1:end),out.y(1,2*Nc+1:end)','AbsTol',1e-12)
end

function testSpinUpIsIndependentOfStart(testCase)
% the cycle map is contractive: starting from v = Vpl at the plate rake or
% from the locked-asperity steady state gives the same primed state
sys = meshSystem(testCase,0.1);
Nc = sys.Nc;
event = unicycle.ode.periodicEvent(sys,50*testCase.TestData.yr);
opts = struct('maxCycles',1000,'verbose',0,'tolSlip',1e-7,'tolZeta',1e-7,'tolRake',1e-6);
[y1,info1] = unicycle.ode.spinUp(sys,event,opts);
opts.y0 = unicycle.ode.steadyState(sys);
[y2,info2] = unicycle.ode.spinUp(sys,event,opts);
verifyTrue(testCase,info1.converged && info2.converged)
verifyLessThan(testCase,[info1.contraction,info2.contraction],1)
verifyEqual(testCase,y1(2*Nc+1:3*Nc),y2(2*Nc+1:3*Nc),'AbsTol',1e-5)
verifyEqual(testCase,angle(exp(1i*(y1(3*Nc+1:end) - y2(3*Nc+1:end)))),zeros(Nc,1),'AbsTol',1e-5)
end

function testSpinUpRejectsPhaseOutsideCycle(testCase)
sys = meshSystem(testCase,0.6);
event = unicycle.ode.periodicEvent(sys,100*testCase.TestData.yr);
verifyError(testCase,@() unicycle.ode.spinUp(sys,event,struct('tSince',2*event.Trec,'verbose',0)),?MException)
end

function testSequenceOfOneEventMatchesEvent(testCase)
% a periodic event is a sequence of one event at t = 0 with T = Trec; dtau
% is computed from the slip when not given
yr = testCase.TestData.yr;
sys = meshSystem(testCase,0.1);
Trec = 50*yr;
event = unicycle.ode.periodicEvent(sys,Trec);
opts = struct('maxCycles',500,'verbose',0,'tSince',20*yr);
[y1,info1] = unicycle.ode.spinUp(sys,event,opts);
seq = struct('T',Trec,'events',struct('t',0,'slipP',event.slipP,'slipC',event.slipC));
[y2,info2] = unicycle.ode.spinUp(sys,seq,opts);
verifyEqual(testCase,info2.seq.events.dtau,event.dtau,'AbsTol',1e-12*max(abs(event.dtau)))
verifyEqual(testCase,info2.nCycles,info1.nCycles)
verifyEqual(testCase,y2,y1,'AbsTol',1e-9)
verifyEqual(testCase,info2.slipSince,info1.slipSince,'AbsTol',1e-9)
end

function testSinglePassMatchesRepeatedEvent(testCase)
% one pass (repeat = false) through the characteristic event repeated
% every Trec gives the same pre-event states as the spin-up loop
yr = testCase.TestData.yr;
sys = meshSystem(testCase,0.1);
Nc = sys.Nc;
Trec = 50*yr;
event = unicycle.ode.periodicEvent(sys,Trec);
nCyc = 4;
w = warning('off','unicycle:ode:spinUp:notConverged');
c = onCleanup(@() warning(w));
[~,info] = unicycle.ode.spinUp(sys,event,struct('maxCycles',nCyc,'tolSlip',0,...
    'keepStates',true,'verbose',0));
ev = repmat(struct('t',0,'dtau',event.dtau,'slipP',event.slipP,'slipC',event.slipC),nCyc,1);
for k = 1:nCyc, ev(k).t = (k-1)*Trec; end
seq = struct('T',nCyc*Trec,'events',ev);
w2 = warning('off','unicycle:ode:spinUp:lockedBudget');
c2 = onCleanup(@() warning(w2));
[y,pass] = unicycle.ode.spinUp(sys,seq,struct('repeat',false,'verbose',0));
verifyTrue(testCase,pass.converged)
verifyEqual(testCase,pass.nCycles,1)
verifyEqual(testCase,pass.zetaPre,info.yStart(2*Nc+1:3*Nc,1:nCyc),'AbsTol',1e-6)
dTheta = angle(exp(1i*(pass.thetaPre - info.yStart(3*Nc+1:end,1:nCyc))));
verifyLessThan(testCase,max(abs(dTheta(:))),1e-6)
verifyEqual(testCase,y(2*Nc+1:3*Nc),info.yStart(2*Nc+1:3*Nc,end),'AbsTol',1e-6)
end

function testSequenceSupercycleClosesBudget(testCase)
% a periodic sequence whose events release Vpl T on the locked patches
% (two halves at different times, one in two parts, plus an event on
% creeping patches only) converges to a supercycle in which every creeping
% patch slips Vpl T; dropping an event warns
yr = testCase.TestData.yr;
[sys,rcv] = meshSystem(testCase,0.1);
Nc = sys.Nc;
T = 60*yr;
VplP = sys.Vpl(sys.iP,:)*T;
west = rcv.xc(sys.iP,1) < 0;
deep = rcv.xc(sys.iC,2) > 24e3;
mk = @(t,sP,sC) struct('t',t,'slipP',sP,'slipC',sC);
ev(1) = mk(0,VplP.*west,zeros(Nc,2));
ev(2) = mk(35*yr,0.4*VplP.*~west,0.3*sys.Vpl(sys.iC,:)*T.*deep);
ev(3) = mk(20*yr,0.6*VplP.*~west,zeros(Nc,2));
seq = struct('T',T,'events',ev);
out = unicycle.ode.runCycles(sys,seq,struct('maxCycles',500,'verbose',0,'nOut',41));
verifyTrue(testCase,out.converged)
VplC = sys.Vpl(sys.iC,:)*T;
verifyEqual(testCase,out.slipTarget,VplC,'AbsTol',1e-9*max(abs(VplC(:))))
verifyEqual(testCase,out.slipCycle,VplC,'AbsTol',2e-3*max(abs(VplC(:))))
verifyEqual(testCase,out.tEvents,[0;20;35]*yr)
% stored supercycle: ordered in time, events appear at their times, ends
% where it started
verifyTrue(testCase,all(diff(out.t) >= 0))
verifyEqual(testCase,nnz(out.t == 20*yr),2)
verifyEqual(testCase,out.t(end),T)
verifyEqual(testCase,out.y(end,1:2*Nc)',out.slipCycle(:),'AbsTol',1e-12)
verifyLessThan(testCase,max(abs(out.y(end,2*Nc+1:3*Nc) - out.yStart(2*Nc+1:3*Nc,end)')),2e-3)
% pre-event state of the event at 20 yr = dense state at 20 yr
i20 = find(out.t == 20*yr,1);
verifyEqual(testCase,out.zetaPre(:,2),out.y(i20,2*Nc+1:3*Nc)','AbsTol',1e-12)
seq.events = ev([1,3]);
verifyWarning(testCase,@() unicycle.ode.spinUp(sys,seq,struct('maxCycles',1,'verbose',0)),...
    'unicycle:ode:spinUp:lockedBudget')
end

function testSinglePassWithoutLockedZone(testCase)
% no locked zone, one pass through a non-repeating sequence: the jumps are
% applyEvent of the pre-event states, slip adds up, dense output joins the
% event states, and stopping at tSince truncates the pass consistently
yr = testCase.TestData.yr;
[sys,rcv] = meshSystem(testCase,0.3,false);
verifyEmpty(testCase,sys.iP)
Nc = sys.Nc;
seq = squareSequence(sys,rcv,[0,20,30,150,170,180],[2,3,2.5,3,2,2.5],300);
[y,info] = unicycle.ode.spinUp(sys,seq,struct('repeat',false,'dense',[0,1,2,6],'nOut',41,'verbose',0));
M = numel(seq.events);
for k = 1:M
    yk = [zeros(2*Nc,1);info.zetaPre(:,k);info.thetaPre(:,k)];
    yP = unicycle.ode.applyEvent(sys,yk,seq.events(k).dtau,seq.events(k).slipC);
    verifyEqual(testCase,info.zetaPost(:,k),yP(2*Nc+1:3*Nc),'AbsTol',1e-12)
end
% creep between events is forward along the plate rake on average
slipC = cat(3,seq.events.slipC);
dA = info.slipPre(:,:,2:end) - info.slipPre(:,:,1:end-1) - slipC(:,:,1:end-1);
along = dA(:,1,:)*cosd(rcv.Vrake) + dA(:,2,:)*sind(rcv.Vrake);
verifyGreaterThan(testCase,sum(along(:)),0)
% dense output: interval 0 is empty (event at t = 0); others join the states
verifyEqual(testCase,[info.dense.interval],[1,2,6])
for d = info.dense(1:2)
    k = d.interval;
    verifyEqual(testCase,d.t([1,end]),[seq.events(k).t;seq.events(k+1).t],'RelTol',1e-12)
    verifyEqual(testCase,reshape(d.y(1,1:2*Nc),Nc,2),info.slipPre(:,:,k) + seq.events(k).slipC,'AbsTol',1e-12)
    verifyEqual(testCase,reshape(d.y(end,1:2*Nc),Nc,2),info.slipPre(:,:,k+1),'AbsTol',1e-9)
    verifyEqual(testCase,d.y(end,2*Nc+1:3*Nc)',info.zetaPre(:,k+1),'AbsTol',1e-12)
end
d = info.dense(3);
verifyEqual(testCase,d.t(end),seq.T,'RelTol',1e-12)
verifyEqual(testCase,reshape(d.y(end,1:2*Nc),Nc,2),info.slipSince,'AbsTol',1e-9)
verifyEqual(testCase,d.y(end,2*Nc+1:end)',y(2*Nc+1:end),'AbsTol',1e-12)
% stop at 160 yr: events after it are not applied (dense output changes
% the first ode45 step, hence 1e-8)
[y2,info2] = unicycle.ode.spinUp(sys,seq,struct('repeat',false,'tSince',160*yr,'verbose',0));
verifyEqual(testCase,info2.zetaPre(:,1:4),info.zetaPre(:,1:4),'AbsTol',1e-8)
verifyTrue(testCase,all(isnan(info2.zetaPre(:,5:6)),'all'))
[y3,info3] = unicycle.ode.spinUp(sys,seq,struct('repeat',false,'tSince',170*yr,'verbose',0));
verifyEqual(testCase,info3.zetaPre(:,5),info.zetaPre(:,5),'AbsTol',1e-8)
verifyEqual(testCase,y3(2*Nc+1:3*Nc),info.zetaPost(:,5),'AbsTol',1e-8)
verifyEqual(testCase,info3.slipSince,info.slipPre(:,:,5) + seq.events(5).slipC,'AbsTol',1e-8)
verifyEqual(testCase,y2(1:2*Nc),zeros(2*Nc,1))
end
