function tests = testImposedEvent
% imposed earthquakes: unicycle.ode.eventStressChange, stateTraction,
% applyEvent and steadyState
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

function setup(testCase)
rcv = unicycle.geometry.triangleReceiver(testCase.TestData.base,unicycle.greens.nikkhoo15(30e3,0.25));
rcv = setTestProperties(rcv);
% a small asperity in the middle of the mesh
rcv.pinnedPosition = find(hypot(rcv.xc(:,1),rcv.xc(:,2)-14e3) < 11e3);
testCase.TestData.rcv = rcv;
end

function testKcpPartition(testCase)
rcv = testCase.TestData.rcv;
bigK = testCase.TestData.bigK;
sys = unicycle.ode.assembleODE(rcv,bigK);
N = rcv.N;
iP = sort(rcv.pinnedPosition(:));
verifyEqual(testCase,sys.iP,iP)
verifyEqual(testCase,sys.Kcp,bigK([sys.iC;N+sys.iC],[iP;N+iP]))
end

function testStressChangeMatchesKernel(testCase)
% Δτ_c = K_cp slip_P with slip_P = Vpl Trec along the plate rake
rcv = testCase.TestData.rcv;
bigK = testCase.TestData.bigK;
sys = unicycle.ode.assembleODE(rcv,bigK);
N = rcv.N;
Trec = 100*testCase.TestData.yr;
[dtau,slipP] = unicycle.ode.eventStressChange(sys,Trec);
Np = numel(sys.iP);
verifyEqual(testCase,slipP,repmat(rcv.Vpl*Trec*[cosd(70),sind(70)],Np,1),'RelTol',1e-12)
% independent evaluation from the full kernel with zero slip elsewhere
slipAll = zeros(N,2);
slipAll(sys.iP,:) = slipP;
dtauAll = bigK*slipAll(:);
verifyEqual(testCase,dtau,dtauAll([sys.iC;N+sys.iC]),'AbsTol',1e-12*max(abs(dtau)))
% linear in the weight
w = linspace(0.2,1,Np)';
[dtauW,slipW] = unicycle.ode.eventStressChange(sys,Trec,w);
verifyEqual(testCase,slipW,w.*slipP,'RelTol',1e-12)
verifyEqual(testCase,dtauW,sys.Kcp*slipW(:),'AbsTol',1e-12*max(abs(dtau)))
end

function testEventLoadsNearestNeighbours(testCase)
% creeping patches next to a slipping asperity are loaded in the slip
% direction
rcv = testCase.TestData.rcv;
sys = unicycle.ode.assembleODE(rcv,testCase.TestData.bigK);
Nc = sys.Nc;
dtau = unicycle.ode.eventStressChange(sys,100*testCase.TestData.yr);
xP = rcv.xc(sys.iP,:);
xC = rcv.xc(sys.iC,:);
dmin = min(sqrt(sum((permute(xC,[1,3,2]) - permute(xP,[3,1,2])).^2,3)),[],2);
near = dmin < 7e3;
verifyGreaterThan(testCase,sum(near),0)
slipDir = [cosd(70);sind(70)];
dtauPar = [dtau(1:Nc),dtau(Nc+1:end)]*slipDir;
verifyGreaterThan(testCase,dtauPar(near),0)
end

function testZeroStressChangeLeavesStateUnchanged(testCase)
% stateTraction and applyEvent are exact inverses of each other
for law = {'log','asinh'}
    sys = oneSliderSystem(0.1,1e-9,deg2rad(30),law{1},0.6,5);
    for zeta = [log(1e-12/sys.v0),sys.y0(3),log(1e-2/sys.v0)]
        y = [1;-2;zeta;deg2rad(-40)];
        verifyEqual(testCase,unicycle.ode.applyEvent(sys,y,[0;0]),y,'AbsTol',1e-10)
    end
end
end

function testTractionBalanceAfterEvent(testCase)
% F(v⁺) + η v⁺ = |T⁻ + Δτ| and v⁺ ∥ T⁻ + Δτ, for steps that rotate,
% amplify and reduce the traction
eta = 30e3/(2*3e3);
for law = {'log','asinh'}
    sys = oneSliderSystem(0.1,1e-9,deg2rad(30),law{1},0.6,eta);
    y = sys.y0;
    Tm = unicycle.ode.stateTraction(sys,y);
    for dtau = {[3;0],[0;3],[-2;1],[10;10],[0.01;-0.02]}
        yp = unicycle.ode.applyEvent(sys,y,dtau{1});
        Tp = Tm + dtau{1}';
        vp = sys.v0*exp(yp(3));
        verifyEqual(testCase,unicycle.ode.frictionLaw(sys,yp(3)) + eta*vp,norm(Tp),'RelTol',1e-12)
        verifyEqual(testCase,yp(4),atan2(Tp(2),Tp(1)),'AbsTol',1e-12)
        verifyEqual(testCase,yp(1:2),y(1:2))
        verifyEqual(testCase,unicycle.ode.stateTraction(sys,yp),Tp,'RelTol',1e-12)
    end
end
end

function testRadiationDampingLimitsCoseismicVelocity(testCase)
% a large stress step gives v⁺ < |T⁺|/η, far below the undamped value
eta = 30e3/(2*3e3);
sys = oneSliderSystem(0.1,1e-9,0,'asinh',0.6,eta);
yp = unicycle.ode.applyEvent(sys,sys.y0,[15;0]);
Tp = unicycle.ode.stateTraction(sys,sys.y0) + [15,0];
vp = sys.v0*exp(yp(3));
vUndamped = sys.v0*exp(unicycle.ode.frictionInverse(sys,norm(Tp),false));
verifyLessThan(testCase,vp,norm(Tp)/eta)
verifyLessThan(testCase,vp,1e-3*vUndamped)
end

function testSpringSliderAfterslipReleasesStressStep(testCase)
% after a stress step Δτ aligned with sliding, afterslip relaxes the
% spring back to steady state: slip - Vpl t -> Δτ/k
yr = testCase.TestData.yr;
eta = 30e3/(2*3e3);
k = 0.1; Vpl = 1e-9; dtau = 3;
opts = odeset('RelTol',1e-10,'AbsTol',1e-10);
for law = {'log','asinh'}
    sys = oneSliderSystem(k,Vpl,deg2rad(30),law{1},0.6,eta);
    y0 = unicycle.ode.applyEvent(sys,sys.y0,dtau*[cosd(30);sind(30)]);
    % relaxation time (a-b)σ/(k Vpl) ~ 80 yr
    tEnd = 3000*yr;
    [~,y] = ode45(@(t,y) unicycle.ode.odeRateDependent(t,y,sys),[0,tEnd],y0,opts);
    afterslip = hypot(y(end,1),y(end,2)) - Vpl*tEnd;
    verifyEqual(testCase,afterslip,dtau/k,'RelTol',1e-4)
    verifyEqual(testCase,y(end,4),deg2rad(30),'AbsTol',1e-8)
end
end

function testSteadyStateIsFixedPoint(testCase)
rcv = testCase.TestData.rcv;
sys = unicycle.ode.assembleODE(rcv,testCase.TestData.bigK);
Nc = sys.Nc;
[y,v] = unicycle.ode.steadyState(sys);
verifyEqual(testCase,v(:),sys.Kcc \ (-sys.taudotLoad),'RelTol',1e-12)
dydt = unicycle.ode.odeRateDependent(0,y,sys);
verifyEqual(testCase,dydt(1:2*Nc),v(:),'RelTol',1e-12)
% natural rate scale k Vpl/((a-b)σ) with k ~ G/dx
rate = 30e3/10e3*rcv.Vpl/0.25;
verifyLessThan(testCase,max(abs(dydt(2*Nc+1:end))),1e-9*rate)
end

function testNegativeLogTractionErrors(testCase)
% the log law has F < 0 below v0 exp(-f0/(a-b)); the traction vector is
% then undefined and the event must not be applied silently
sys = oneSliderSystem(0.1,1e-9,0,'log',0.01,0);
y = [0;0;log(1e-12/sys.v0);0];
verifyError(testCase,@() unicycle.ode.applyEvent(sys,y,[1;0]),...
    'unicycle:ode:stateTraction:negativeTraction')
end

function testRuptureOntoCreepingPatches(testCase)
% coseismic slip on creeping patches adds K_cc slip_C to Δτ (a stress drop
% on the slipping creeping patches) and jumps their slip state
rcv = testCase.TestData.rcv;
bigK = testCase.TestData.bigK;
sys = unicycle.ode.assembleODE(rcv,bigK);
N = rcv.N; Nc = sys.Nc;
Trec = 100*testCase.TestData.yr;
% rupture also covers the creeping patches nearest the asperity
xP = rcv.xc(sys.iP,:);
xC = rcv.xc(sys.iC,:);
dmin = min(sqrt(sum((permute(xC,[1,3,2]) - permute(xP,[3,1,2])).^2,3)),[],2);
weightC = double(dmin < 7e3);
[dtau,slipP,slipC] = unicycle.ode.eventStressChange(sys,Trec,1,weightC);
verifyEqual(testCase,slipC,weightC.*sys.Vpl(sys.iC,:)*Trec,'RelTol',1e-12)
slipAll = zeros(N,2);
slipAll(sys.iP,:) = slipP;
slipAll(sys.iC,:) = slipC;
dtauAll = bigK*slipAll(:);
verifyEqual(testCase,dtau,dtauAll([sys.iC;N+sys.iC]),'AbsTol',1e-12*max(abs(dtau)))
% the rupturing creeping patches have a stress drop along the slip direction
slipDir = [cosd(70);sind(70)];
dtauPar = [dtau(1:Nc),dtau(Nc+1:end)]*slipDir;
verifyLessThan(testCase,dtauPar(weightC > 0),0)
% applyEvent adds the coseismic slip to the state
y0 = unicycle.ode.steadyState(sys);
y1 = unicycle.ode.applyEvent(sys,y0,dtau,slipC);
verifyEqual(testCase,y1(1:2*Nc),y0(1:2*Nc) + slipC(:))
end

function testTractionReversal(testCase)
% a stress drop larger than the traction reverses it; with the asinh law
% the patch slips backward with v ∥ T⁺
eta = 30e3/(2*3e3);
sys = oneSliderSystem(0.1,1e-9,deg2rad(30),'asinh',0.1,eta);
Tm = unicycle.ode.stateTraction(sys,sys.y0);
dtau = -2*Tm' + [0;1e-3];
yp = unicycle.ode.applyEvent(sys,sys.y0,dtau);
Tp = Tm + dtau';
vp = sys.v0*exp(yp(3))*[cos(yp(4)),sin(yp(4))];
verifyLessThan(testCase,vp*[cosd(30);sind(30)],0)
verifyEqual(testCase,vp/norm(vp),Tp/norm(Tp),'AbsTol',1e-12)
verifyEqual(testCase,unicycle.ode.frictionLaw(sys,yp(3)) + eta*norm(vp),norm(Tp),'RelTol',1e-12)
end

function testRuptureSlipForStressChange(testCase)
% slip on the creeping rupture patches gives exactly the target traction
% change along the pre-event traction on each of them
rcv = testCase.TestData.rcv;
sys = unicycle.ode.assembleODE(rcv,testCase.TestData.bigK);
Nc = sys.Nc;
Trec = 100*testCase.TestData.yr;
xP = rcv.xc(sys.iP,:);
xC = rcv.xc(sys.iC,:);
dmin = min(sqrt(sum((permute(xC,[1,3,2]) - permute(xP,[3,1,2])).^2,3)),[],2);
inE = dmin < 7e3;
y = unicycle.ode.steadyState(sys);
d = unicycle.ode.stateTraction(sys,y);
target = -0.5;
[wC,dtauPar] = unicycle.ode.ruptureSlipForStressChange(sys,Trec,inE,target,d);
verifyEqual(testCase,wC(~inE),zeros(sum(~inE),1))
verifyEqual(testCase,dtauPar(inE),target*ones(sum(inE),1),'AbsTol',1e-9)
% consistent with eventStressChange
dtau = unicycle.ode.eventStressChange(sys,Trec,1,wC);
dn = d./hypot(d(:,1),d(:,2));
verifyEqual(testCase,dn(:,1).*dtau(1:Nc) + dn(:,2).*dtau(Nc+1:end),dtauPar,'AbsTol',1e-9)
end
