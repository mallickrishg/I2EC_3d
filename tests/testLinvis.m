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
% v̇ = D⁻¹(K v + τ̇_load), ṡ = v, for heterogeneous eta
yr = testCase.TestData.yr;
rcv = meshReceiver(testCase,true);
lin = unicycle.ode.linvisAssemble(rcv,testCase.TestData.bigK,1);
lin = unicycle.ode.linvisAssemble(rcv,testCase.TestData.bigK,heterogeneousEta(testCase,rcv,lin));
n = 2*lin.Nc;
v0 = lin.vStar + 3*rcv.Vpl*cos((1:n)');
t = [0,logspace(-3,3,30)]'*yr;
[v,s] = unicycle.ode.linvisEvolve(lin,v0,t);
rhs = @(t,y) [(lin.Kcc*y(1:n) + lin.taudotLoad)./lin.d;y(1:n)];
[~,Y] = ode45(rhs,t,[v0;zeros(n,1)],odeset('RelTol',1e-11,'AbsTol',1e-14));
verifyEqual(testCase,v',Y(:,1:n),'AbsTol',1e-7*rcv.Vpl)
verifyEqual(testCase,s',Y(:,n+1:end),'AbsTol',1e-7*rcv.Vpl*t(end))
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
