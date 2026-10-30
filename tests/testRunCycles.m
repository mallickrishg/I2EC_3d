function tests = testRunCycles
% periodic imposed earthquake cycles: unicycle.ode.periodicEvent and
% unicycle.ode.runCycles
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

function sys = meshSystem(testCase,f0)
rcv = unicycle.geometry.triangleReceiver(testCase.TestData.base,unicycle.greens.nikkhoo15(30e3,0.25));
rcv = setTestProperties(rcv);
rcv.mu0 = f0;
rcv.Vo = rcv.Vpl;
rcv.pinnedPosition = find(hypot(rcv.xc(:,1),rcv.xc(:,2)-14e3) < 11e3);
sys = unicycle.ode.assembleODE(rcv,testCase.TestData.bigK,'asinh');
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
    'unicycle:ode:runCycles:lockedBudget')
end
