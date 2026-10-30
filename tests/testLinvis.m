function tests = testLinvis
% analytic solution of the linear viscous cycle ODE:
% unicycle.ode.linvisAssemble, linvisEvolve, linvisEvent, linvisSpinUp
tests = functiontests(localfunctions);
end

function setupOnce(testCase)
addpath(fileparts(fileparts(mfilename('fullpath'))));
yr = 3.1536e7;
testCase.TestData.yr = yr;
testCase.TestData.dir = tempname;
mkdir(testCase.TestData.dir);
testCase.TestData.base = makeTestMesh(testCase.TestData.dir,'plane');
rcv = unicycle.geometry.triangleReceiver(testCase.TestData.base,unicycle.greens.nikkhoo15(30e3,0.25));
testCase.TestData.bigK = unicycle.greens.computeKernels(rcv,'plane',fullfile(testCase.TestData.dir,'kernels'));
% tangent viscosity of RSF at plate rate, (a-b) sigma / Vpl
testCase.TestData.eta0 = 0.25/(0.05/yr);
end

function teardownOnce(testCase)
rmdir(testCase.TestData.dir,'s');
end

function rcv = meshReceiver(testCase,pinned)
rcv = unicycle.geometry.triangleReceiver(testCase.TestData.base,unicycle.greens.nikkhoo15(30e3,0.25));
rcv = setTestProperties(rcv);
if pinned
    rcv.pinnedPosition = find(hypot(rcv.xc(:,1),rcv.xc(:,2)-14e3) < 11e3);
end
end

function eta = heterogeneousEta(testCase,rcv,lin)
% viscosity increasing with depth, varying by a factor of ~10
z = -rcv.xc(lin.iC,3);
eta = testCase.TestData.eta0*10.^(z/max(z));
end

function testAssemblyMatchesODEPartition(testCase)
rcv = meshReceiver(testCase,true);
bigK = testCase.TestData.bigK;
lin = unicycle.ode.linvisAssemble(rcv,bigK,testCase.TestData.eta0);
sys = unicycle.ode.assembleODE(rcv,bigK);
verifyEqual(testCase,lin.iC,sys.iC)
verifyEqual(testCase,lin.iP,sys.iP)
verifyEqual(testCase,lin.Kcc,sys.Kcc)
verifyEqual(testCase,lin.Kcp,sys.Kcp)
verifyEqual(testCase,lin.taudotLoad,sys.taudotLoad,'RelTol',1e-12)
[~,vSS] = unicycle.ode.steadyState(sys);
verifyEqual(testCase,lin.vStar,vSS(:),'RelTol',1e-10)
% A = D⁻¹K decays and its eigendecomposition reproduces it
verifyLessThan(testCase,real(lin.lambda),0)
A = lin.Kcc./lin.d;
verifyEqual(testCase,real(lin.V*diag(lin.lambda)*lin.Vinv),A,'AbsTol',1e-10*max(abs(A(:))))
% eta given on all patches or on creeping patches only
etaAll = (1:rcv.N)'*testCase.TestData.eta0;
lin2 = unicycle.ode.linvisAssemble(rcv,bigK,etaAll);
verifyEqual(testCase,lin2.eta,etaAll(lin.iC))
verifyError(testCase,@() unicycle.ode.linvisAssemble(rcv,bigK,-1),?MException)
end

function testEvolveMatchesOde45(testCase)
% the analytic slip rate and slip match a numerical solution of
% v̇ = D⁻¹(K v + τ̇_load), ṡ = v (unicycle.ode.odeLinearViscous), for
% heterogeneous eta
yr = testCase.TestData.yr;
rcv = meshReceiver(testCase,true);
lin = unicycle.ode.linvisAssemble(rcv,testCase.TestData.bigK,1);
lin = unicycle.ode.linvisAssemble(rcv,testCase.TestData.bigK,heterogeneousEta(testCase,rcv,lin));
n = 2*lin.Nc;
v0 = lin.vStar + 3*rcv.Vpl*cos((1:n)');
t = [0,logspace(-3,3,30)]'*yr;
[v,s] = unicycle.ode.linvisEvolve(lin,v0,t);
rhs = @(t,y) unicycle.ode.odeLinearViscous(t,y,lin);
[~,Y] = ode45(rhs,t,[zeros(n,1);v0],odeset('RelTol',1e-11,'AbsTol',1e-14));
verifyEqual(testCase,v',Y(:,n+1:end),'AbsTol',1e-7*rcv.Vpl)
verifyEqual(testCase,s',Y(:,1:n),'AbsTol',1e-7*rcv.Vpl*t(end))
% t = 0 returns the initial condition exactly
verifyEqual(testCase,v(:,1),v0,'AbsTol',1e-12*max(abs(v0)))
verifyEqual(testCase,s(:,1),zeros(n,1))
end

function testSteadyStateIsFixedPoint(testCase)
rcv = meshReceiver(testCase,true);
lin = unicycle.ode.linvisAssemble(rcv,testCase.TestData.bigK,testCase.TestData.eta0);
t = [0,1,10,1000]*testCase.TestData.yr;
[v,s] = unicycle.ode.linvisEvolve(lin,lin.vStar,t);
verifyEqual(testCase,v,repmat(lin.vStar,1,numel(t)),'AbsTol',1e-10*rcv.Vpl)
verifyEqual(testCase,s,lin.vStar.*t,'AbsTol',1e-10*rcv.Vpl*t(end))
end

function testViscosityScalesTime(testCase)
% uniform eta: multiplying eta by c is the same as dividing time by c
yr = testCase.TestData.yr;
rcv = meshReceiver(testCase,true);
lin1 = unicycle.ode.linvisAssemble(rcv,testCase.TestData.bigK,testCase.TestData.eta0);
lin3 = unicycle.ode.linvisAssemble(rcv,testCase.TestData.bigK,3*testCase.TestData.eta0);
v0 = lin1.vStar + rcv.Vpl;
t = [0.5,5,50]*yr;
verifyEqual(testCase,unicycle.ode.linvisEvolve(lin3,v0,3*t),unicycle.ode.linvisEvolve(lin1,v0,t),...
    'AbsTol',1e-10*rcv.Vpl)
end

function testEventJump(testCase)
rcv = meshReceiver(testCase,true);
lin = unicycle.ode.linvisAssemble(rcv,testCase.TestData.bigK,testCase.TestData.eta0);
dtau = linspace(-1,2,2*lin.Nc)';
v = unicycle.ode.linvisEvent(lin,lin.vStar,dtau);
verifyEqual(testCase,v,lin.vStar + dtau./[lin.eta;lin.eta],'RelTol',1e-14)
end

function testSpinUpIsPeriodic(testCase)
% the closed-form pre-event state returns to itself after one event and
% one interseismic period, and equals the limit of iterating the cycle map
yr = testCase.TestData.yr;
rcv = meshReceiver(testCase,true);
lin = unicycle.ode.linvisAssemble(rcv,testCase.TestData.bigK,1);
lin = unicycle.ode.linvisAssemble(rcv,testCase.TestData.bigK,heterogeneousEta(testCase,rcv,lin));
Trec = 50*yr;
event = unicycle.ode.periodicEvent(lin,Trec);
[vPre,info] = unicycle.ode.linvisSpinUp(lin,event);
verifyEqual(testCase,vPre,info.vPre,'AbsTol',1e-12*rcv.Vpl)
vEnd = unicycle.ode.linvisEvolve(lin,unicycle.ode.linvisEvent(lin,vPre,event.dtau),Trec);
verifyEqual(testCase,vEnd,vPre,'AbsTol',1e-9*rcv.Vpl)
verifyLessThan(testCase,info.spectralRadius,1)
% iterate the cycle map from v = Vpl
v = reshape(lin.Vpl(lin.iC,:),[],1);
for k = 1:200
    v = unicycle.ode.linvisEvolve(lin,unicycle.ode.linvisEvent(lin,v,event.dtau),Trec);
end
verifyEqual(testCase,v,vPre,'AbsTol',1e-6*rcv.Vpl)
% state part way through the cycle
tS = 7*yr;
[vS,infoS] = unicycle.ode.linvisSpinUp(lin,event,tS);
[vRef,sRef] = unicycle.ode.linvisEvolve(lin,info.vPost,tS);
verifyEqual(testCase,vS,vRef,'AbsTol',1e-12*rcv.Vpl)
verifyEqual(testCase,infoS.slipSince(:),sRef,'AbsTol',1e-12)
end

function testSlipBudgetClosesIndependentOfViscosity(testCase)
% slip per cycle = Vpl Trec on every creeping patch, for any eta, with or
% without coseismic slip on creeping patches
yr = testCase.TestData.yr;
rcv = meshReceiver(testCase,true);
lin1 = unicycle.ode.linvisAssemble(rcv,testCase.TestData.bigK,testCase.TestData.eta0);
lin2 = unicycle.ode.linvisAssemble(rcv,testCase.TestData.bigK,heterogeneousEta(testCase,rcv,lin1));
Trec = 50*yr;
VplC = lin1.Vpl(lin1.iC,:)*Trec;
for wC = {0,linspace(0,0.5,lin1.Nc)'}
    for lin = {lin1,lin2}
        event = unicycle.ode.periodicEvent(lin{1},Trec,wC{1});
        [~,info] = unicycle.ode.linvisSpinUp(lin{1},event);
        verifyEqual(testCase,info.slipCycle,VplC,'AbsTol',1e-9*max(abs(VplC(:))))
    end
end
end

function testSurfaceDeformationAcceptsSlipHistory(testCase)
% the linear solution feeds unicycle.ode.surfaceDeformation directly, and
% the relative displacement returns to zero after a periodic cycle
yr = testCase.TestData.yr;
rcv = meshReceiver(testCase,true);
lin = unicycle.ode.linvisAssemble(rcv,testCase.TestData.bigK,testCase.TestData.eta0);
Trec = 50*yr;
event = unicycle.ode.periodicEvent(lin,Trec);
[~,info] = unicycle.ode.linvisSpinUp(lin,event);
t = linspace(0,Trec,21)';
[v,s] = unicycle.ode.linvisEvolve(lin,info.vPost,t);
G = unicycle.greens.computeSurfaceKernels(rcv,[0,14e3,0;0,40e3,0]);
o = struct('removeLongTerm',true,'pinnedSlip',event.slipP);
[vS,uS] = unicycle.ode.surfaceDeformation(G,lin,t,struct('v',v,'s',s),o);
verifyEqual(testCase,size(vS),[numel(t),2,3])
verifyEqual(testCase,uS(end,:,:),zeros(1,2,3),'AbsTol',1e-6*max(abs(uS(:))))
end

function testViscosityLimits(testCase)
% eta -> infinity: steady creep at Vpl through the whole cycle (the events
% are averaged out); eta -> 0: instant relaxation to the steady state v*
yr = testCase.TestData.yr;
rcv = meshReceiver(testCase,true);
Trec = 50*yr;
VplC = [];
for etaScale = [1e8,1e12]
    lin = unicycle.ode.linvisAssemble(rcv,testCase.TestData.bigK,etaScale*testCase.TestData.eta0);
    VplC = reshape(lin.Vpl(lin.iC,:),[],1);
    event = unicycle.ode.periodicEvent(lin,Trec);
    [v,info] = unicycle.ode.linvisSpinUp(lin,event,Trec/3);
    verifyEqual(testCase,info.vPost,VplC,'AbsTol',1e-6*rcv.Vpl)
    verifyEqual(testCase,v,VplC,'AbsTol',1e-6*rcv.Vpl)
    verifyEqual(testCase,info.slipCycle(:),VplC*Trec,'AbsTol',1e-6*rcv.Vpl*Trec)
end
lin = unicycle.ode.linvisAssemble(rcv,testCase.TestData.bigK,1e-6*testCase.TestData.eta0);
event = unicycle.ode.periodicEvent(lin,Trec);
[v,info] = unicycle.ode.linvisSpinUp(lin,event,Trec/3);
verifyEqual(testCase,v,lin.vStar,'AbsTol',1e-9*rcv.Vpl)
verifyEqual(testCase,info.vPre,lin.vStar,'AbsTol',1e-9*rcv.Vpl)
verifyEqual(testCase,info.slipCycle(:),VplC*Trec,'AbsTol',1e-6*rcv.Vpl*Trec)
end

function testAnalyticSpinUpMatchesNumericalSpinUp(testCase)
% the closed-form periodic cycle equals the cycle obtained by time stepping
% (event + ode45 over Trec, repeated until periodic), for heterogeneous
% eta, with and without coseismic slip on creeping patches, and for a
% slowly relaxing (high-eta) case that needs many cycles
yr = testCase.TestData.yr;
rcv = meshReceiver(testCase,true);
lin0 = unicycle.ode.linvisAssemble(rcv,testCase.TestData.bigK,1);
etaH = heterogeneousEta(testCase,rcv,lin0);
Trec = 50*yr;
Vpl = rcv.Vpl;
cases = {etaH,0; etaH,linspace(0,0.4,lin0.Nc)'; 30*etaH,0};
for k = 1:size(cases,1)
    lin = unicycle.ode.linvisAssemble(rcv,testCase.TestData.bigK,cases{k,1});
    event = unicycle.ode.periodicEvent(lin,Trec,cases{k,2});
    [~,infoA] = unicycle.ode.linvisSpinUp(lin,event);
    [vPreN,infoN] = unicycle.ode.linvisSpinUpNumerical(lin,event);
    verifyTrue(testCase,infoN.converged)
    verifyEqual(testCase,vPreN,infoA.vPre,'AbsTol',1e-7*Vpl)
    verifyEqual(testCase,infoN.vPost,infoA.vPost,'AbsTol',1e-7*Vpl)
    verifyEqual(testCase,infoN.slipCycle,infoA.slipCycle,'AbsTol',1e-7*Vpl*Trec)
    % the whole periodic cycle, not just its end points
    t = [0,logspace(-3,log10(Trec/yr),40)]'*yr;
    n = 2*lin.Nc;
    [vA,sA] = unicycle.ode.linvisEvolve(lin,infoA.vPost,t);
    [~,Y] = ode45(@(t,y) unicycle.ode.odeLinearViscous(t,y,lin),t,[zeros(n,1);infoN.vPost],...
        odeset('RelTol',1e-10,'AbsTol',1e-14));
    verifyEqual(testCase,Y(:,n+1:end)',vA,'AbsTol',1e-7*Vpl)
    verifyEqual(testCase,Y(:,1:n)',sA,'AbsTol',1e-7*Vpl*Trec)
end
end

function testExtendedRuptureAnalyticMatchesNumerical(testCase)
% rupture 20 % larger in area than the locked zone: uniform slip Vpl Trec
% on the locked patches, half of that on the creeping patches inside the
% extension (the nearest creeping patches, until their area reaches 20 %
% of the locked area). The closed-form periodic cycle must match the
% numerically spun-up one, and the slip budget still closes
yr = testCase.TestData.yr;
rcv = meshReceiver(testCase,true);
lin0 = unicycle.ode.linvisAssemble(rcv,testCase.TestData.bigK,1);
Trec = 50*yr;
Vpl = rcv.Vpl;
iC = lin0.iC;
iP = lin0.iP;
% creeping patches nearest the locked zone
dmin = min(sqrt(sum((permute(rcv.xc(iC,:),[1,3,2]) - permute(rcv.xc(iP,:),[3,1,2])).^2,3)),[],2);
[~,order] = sort(dmin);
areaC = rcv.area(iC);
nExt = find(cumsum(areaC(order)) >= 0.2*sum(rcv.area(iP)),1);
weightC = zeros(lin0.Nc,1);
weightC(order(1:nExt)) = 0.5;
verifyEqual(testCase,sum(areaC(weightC > 0))/sum(rcv.area(iP)),0.2,'AbsTol',0.1)
for etaCase = {heterogeneousEta(testCase,rcv,lin0),10*heterogeneousEta(testCase,rcv,lin0)}
    lin = unicycle.ode.linvisAssemble(rcv,testCase.TestData.bigK,etaCase{1});
    event = unicycle.ode.periodicEvent(lin,Trec,weightC);
    verifyEqual(testCase,hypot(event.slipP(:,1),event.slipP(:,2)),Vpl*Trec*ones(numel(iP),1),'RelTol',1e-12)
    verifyEqual(testCase,hypot(event.slipC(weightC > 0,1),event.slipC(weightC > 0,2)),...
        0.5*Vpl*Trec*ones(nExt,1),'RelTol',1e-12)
    % slip on the extension lowers its traction along the slip direction
    % relative to a rupture confined to the locked zone (with half the slip,
    % loading by the adjacent locked-zone slip can still win on some patches)
    event0 = unicycle.ode.periodicEvent(lin,Trec);
    ext = find(weightC > 0);
    slipDir = lin.Vpl(iC(ext),:)./Vpl;
    par = @(dtau) sum([dtau(ext),dtau(lin.Nc + ext)].*slipDir,2);
    verifyLessThan(testCase,par(event.dtau) - par(event0.dtau),0)

    [~,infoA] = unicycle.ode.linvisSpinUp(lin,event);
    [vPreN,infoN] = unicycle.ode.linvisSpinUpNumerical(lin,event);
    verifyTrue(testCase,infoN.converged)
    verifyEqual(testCase,vPreN,infoA.vPre,'AbsTol',1e-7*Vpl)
    verifyEqual(testCase,infoN.vPost,infoA.vPost,'AbsTol',1e-7*Vpl)
    verifyEqual(testCase,infoN.slipCycle,infoA.slipCycle,'AbsTol',1e-7*Vpl*Trec)
    verifyEqual(testCase,infoA.slipCycle,lin.Vpl(iC,:)*Trec,'AbsTol',1e-9*Vpl*Trec)
    t = [0,logspace(-3,log10(Trec/yr),40)]'*yr;
    n = 2*lin.Nc;
    [vA,sA] = unicycle.ode.linvisEvolve(lin,infoA.vPost,t);
    [~,Y] = ode45(@(t,y) unicycle.ode.odeLinearViscous(t,y,lin),t,[event.slipC(:);infoN.vPost],...
        odeset('RelTol',1e-10,'AbsTol',1e-14));
    verifyEqual(testCase,Y(:,n+1:end)',vA,'AbsTol',1e-7*Vpl)
    verifyEqual(testCase,Y(:,1:n)',sA + event.slipC(:),'AbsTol',1e-7*Vpl*Trec)
end
end

function seq = mixedSequence(rcv,lin,T)
% supercycle of four events: the west half of the locked zone with Vpl T,
% extending onto the nearest creeping patches (half of that); the east half
% in two events (0.6 and 0.4 Vpl T), the second extending onto creeping
% patches too; and an event confined to deep creeping patches. Every
% locked patch slips Vpl T per supercycle
yr = 3.1536e7;
iC = lin.iC;
iP = lin.iP;
west = rcv.xc(:,1) < 0;
dmin = min(sqrt(sum((permute(rcv.xc(iC,:),[1,3,2]) - permute(rcv.xc(iP,:),[3,1,2])).^2,3)),[],2);
near = false(rcv.N,1);
near(iC(dmin < prctile(dmin,10))) = true;
deep = false(rcv.N,1);
z = -rcv.xc(:,3);
deep(iC(z(iC) > prctile(z(iC),85))) = true;
isP = false(rcv.N,1);
isP(iP) = true;
w = {double(isP & west) + 0.5*(near & west),...
     0.6*(isP & ~west),...
     0.3*deep,...
     0.4*(isP & ~west) + 0.5*(near & ~west)};
tk = [0,35,50,80]*yr;
for k = 1:4
    events(k) = unicycle.ode.sequenceEvent(lin,tk(k),w{k}.*lin.Vpl*T); %#ok<AGROW>
end
% given out of order, the solver sorts them
seq = struct('T',T,'events',events([3,1,4,2]));
end

function testSequenceAnalyticMatchesNumerical(testCase)
% a mixed sequence (partial ruptures of the locked zone, ruptures extending
% onto creeping patches, a purely creeping event): the closed-form
% periodic supercycle matches the numerically spun-up one, and the slip
% budget closes on every creeping patch because it closes on the locked ones
yr = testCase.TestData.yr;
rcv = meshReceiver(testCase,true);
lin0 = unicycle.ode.linvisAssemble(rcv,testCase.TestData.bigK,1);
T = 120*yr;
Vpl = rcv.Vpl;
for etaCase = {heterogeneousEta(testCase,rcv,lin0),10*heterogeneousEta(testCase,rcv,lin0)}
    lin = unicycle.ode.linvisAssemble(rcv,testCase.TestData.bigK,etaCase{1});
    seq = mixedSequence(rcv,lin,T);
    [~,infoA] = unicycle.ode.linvisSpinUpSequence(lin,seq);
    verifyEqual(testCase,infoA.t',[0,35,50,80]*yr)
    verifyLessThan(testCase,infoA.budgetErrP,1e-12)
    verifyLessThan(testCase,infoA.budgetErrC,1e-9)
    verifyEqual(testCase,infoA.slipCycle,lin.Vpl(lin.iC,:)*T,'AbsTol',1e-9*Vpl*T)
    [vPreN,infoN] = unicycle.ode.linvisSpinUpNumerical(lin,seq);
    verifyTrue(testCase,infoN.converged)
    verifyEqual(testCase,vPreN,infoA.vPre(:,1),'AbsTol',1e-7*Vpl)
    verifyEqual(testCase,infoN.vPre,infoA.vPre,'AbsTol',1e-7*Vpl)
    verifyEqual(testCase,infoN.vPost,infoA.vPost,'AbsTol',1e-7*Vpl)
    verifyEqual(testCase,infoN.slipCycle,infoA.slipCycle,'AbsTol',1e-7*Vpl*T)
    % output times: post-event state at an event time, the state before the
    % next supercycle at T, and evolution between events
    tOut = [0,20,35,60,T/yr]*yr;
    [v,info] = unicycle.ode.linvisSpinUpSequence(lin,seq,tOut);
    verifyEqual(testCase,v(:,[1,3]),infoA.vPost(:,[1,2]),'AbsTol',1e-12*Vpl)
    verifyEqual(testCase,v(:,5),infoA.vPre(:,1),'AbsTol',1e-12*Vpl)
    [vE,sE] = unicycle.ode.linvisEvolve(lin,infoA.vPost(:,3),10*yr);
    verifyEqual(testCase,v(:,4),vE,'AbsTol',1e-10*Vpl)
    slipTo50 = info.slip(:,4) - sE;      % slip from 0⁻ to just after the event at 50 yr
    [~,s50] = unicycle.ode.linvisSpinUpSequence(lin,seq,50*yr);
    verifyEqual(testCase,slipTo50,s50.slip,'AbsTol',1e-9*Vpl*T)
    verifyEqual(testCase,info.slip(:,5),infoA.slipCycle(:),'AbsTol',1e-12*Vpl*T)
end
end

function testSequenceIsSumOfShiftedSingleEventCycles(testCase)
% periodic Green's function: u(t) of the sequence is the sum over events of
% the single-event periodic cycles of period T, each shifted to its event
% time; a single event at t_k is the linvisSpinUp cycle shifted by t_k
yr = testCase.TestData.yr;
rcv = meshReceiver(testCase,true);
lin = unicycle.ode.linvisAssemble(rcv,testCase.TestData.bigK,1);
lin = unicycle.ode.linvisAssemble(rcv,testCase.TestData.bigK,heterogeneousEta(testCase,rcv,lin));
T = 120*yr;
seq = mixedSequence(rcv,lin,T);
tOut = linspace(0,T,13);
vSeq = unicycle.ode.linvisSpinUpSequence(lin,seq,tOut);
uSum = zeros(size(vSeq));
noCheck = struct('checkBudget',false);
for k = 1:numel(seq.events)
    vk = unicycle.ode.linvisSpinUpSequence(lin,struct('T',T,'events',seq.events(k)),tOut,noCheck);
    uSum = uSum + vk - lin.vStar;
end
verifyEqual(testCase,vSeq - lin.vStar,uSum,'AbsTol',1e-10*rcv.Vpl)
% one event at t_k vs the characteristic cycle shifted by t_k
ev = seq.events([seq.events.t] == 35*yr);
ev1 = struct('Trec',T,'dtau',ev.dtau,'slipC',ev.slipC);
for tS = [0,10,40,90,120]*yr
    vShift = unicycle.ode.linvisSpinUpSequence(lin,struct('T',T,'events',ev),tS,noCheck);
    tSince = tS - ev.t + T*(tS < ev.t);   % tS = T is just before the event at 0 of the next supercycle
    vRef = unicycle.ode.linvisSpinUp(lin,ev1,tSince);
    verifyEqual(testCase,vShift,vRef,'AbsTol',1e-12*rcv.Vpl)
end
end

function testSequenceSlipBudgetCheck(testCase)
% the locked slip budget must close per supercycle; when it does not, the
% creeping slip per supercycle is off by K_cc⁻¹ K_cp (Σ slip_P - Vpl_P T)
yr = testCase.TestData.yr;
rcv = meshReceiver(testCase,true);
lin = unicycle.ode.linvisAssemble(rcv,testCase.TestData.bigK,testCase.TestData.eta0);
T = 120*yr;
seq = mixedSequence(rcv,lin,T);
bad = seq;
bad.events = seq.events([seq.events.t] ~= 0);           % drop the west-half rupture
verifyError(testCase,@() unicycle.ode.linvisSpinUpSequence(lin,bad),...
    'unicycle:ode:linvisSpinUpSequence:budget')
[~,info] = unicycle.ode.linvisSpinUpSequence(lin,bad,[],struct('checkBudget',false));
verifyEqual(testCase,info.budgetErrP,1,'RelTol',1e-12)  % the west half never slips
mis = info.slipCycleP - lin.Vpl(lin.iP,:)*T;
expected = lin.Vpl(lin.iC,:)*T - reshape(lin.Kcc\(lin.Kcp*mis(:)),[],2);
verifyEqual(testCase,info.slipCycle,expected,'AbsTol',1e-9*rcv.Vpl*T)
verifyGreaterThan(testCase,info.budgetErrC,1e-3)
% a missing slipP cannot be checked
noP = rmfield(seq,'events');
noP.events = rmfield(seq.events,'slipP');
verifyError(testCase,@() unicycle.ode.linvisSpinUpSequence(lin,noP),?MException)
verifyWarningFree(testCase,@() unicycle.ode.linvisSpinUpSequence(lin,noP,[],struct('checkBudget',false)))
% large eta: steady creep through the supercycle at the rate that, with
% the coseismic slip on creeping patches, adds up to Vpl T
linH = unicycle.ode.linvisAssemble(rcv,testCase.TestData.bigK,1e10*testCase.TestData.eta0);
seqH = mixedSequence(rcv,linH,T);
v = unicycle.ode.linvisSpinUpSequence(lin,seq,[0,40,T/yr]*yr);
vH = unicycle.ode.linvisSpinUpSequence(linH,seqH,[0,40,T/yr]*yr);
vSteady = reshape(lin.Vpl(lin.iC,:) - sum(cat(3,seqH.events.slipC),3)/T,[],1);
verifyEqual(testCase,vH,repmat(vSteady,1,3),'AbsTol',1e-6*rcv.Vpl)
verifyGreaterThan(testCase,max(abs(v(:) - vH(:))),0.1*rcv.Vpl)
end
