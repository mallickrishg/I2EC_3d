function tests = testSurfaceDeformation
% surface kernels and surface velocity/displacement time series:
% unicycle.greens.computeSurfaceKernels, unicycle.ode.surfaceDeformation
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
% stations [east, north, up]: hanging wall above the asperity, near the
% trench, far landward, and off to the side
testCase.TestData.xObs = [0,14e3,0;0,3e3,0;0,60e3,0;25e3,14e3,0];
testCase.TestData.G = unicycle.greens.computeSurfaceKernels(rcv,testCase.TestData.xObs);
end

function teardownOnce(testCase)
rmdir(testCase.TestData.dir,'s');
end

function [sys,rcv] = meshSystem(testCase,pinned)
rcv = unicycle.geometry.triangleReceiver(testCase.TestData.base,unicycle.greens.nikkhoo15(30e3,0.25));
rcv = setTestProperties(rcv);
rcv.Vo = rcv.Vpl;
if pinned
    rcv.pinnedPosition = find(hypot(rcv.xc(:,1),rcv.xc(:,2)-14e3) < 11e3);
end
sys = unicycle.ode.assembleODE(rcv,testCase.TestData.bigK,'asinh');
end

function testKernelLayoutAndSlipSigns(testCase)
% rows [e1 n1 u1 e2 ...]; uniform thrust moves the hanging wall up-dip
% (toward the trench, -north) and up; positive strike slip moves it along
% the strike vector (west for this mesh)
G = testCase.TestData.G;
[~,rcv] = meshSystem(testCase,false);
N = rcv.N;
D = size(testCase.TestData.xObs,1);
verifyEqual(testCase,size(G),[3*D,2*N])
uThrust = reshape(G*[zeros(N,1);ones(N,1)],3,D);
verifyLessThan(testCase,uThrust(2,1),0)
verifyGreaterThan(testCase,uThrust(3,1),0)
uStrike = reshape(G*[ones(N,1);zeros(N,1)],3,D);
verifyEqual(testCase,rcv.sv(1,:),[-1,0,0],'AbsTol',1e-6)
verifyLessThan(testCase,uStrike(1,1),0)
% consistent with calling displacementKernels directly
[Gs,Gd] = rcv.displacementKernels(testCase.TestData.xObs,3);
verifyEqual(testCase,G,[Gs,Gd])
end

function testSteadySlidingIsRemovedByLongTerm(testCase)
% with no locked patches and v = Vpl everywhere, the deformation equals
% the long-term term exactly
yr = testCase.TestData.yr;
sys = meshSystem(testCase,false);
Nc = sys.Nc;
t = linspace(0,100*yr,11)';
Y = [t*sys.Vpl(:,1)',t*sys.Vpl(:,2)',repmat(sys.y0(2*Nc+1:end)',numel(t),1)];
G = testCase.TestData.G;
[v,u] = unicycle.ode.surfaceDeformation(G,sys,t,Y);
[vr,ur] = unicycle.ode.surfaceDeformation(G,sys,t,Y,struct('removeLongTerm',true));
scaleV = max(abs(v(:)));
scaleU = max(abs(u(:)));
verifyGreaterThan(testCase,scaleV,0)
verifyEqual(testCase,vr,zeros(size(vr)),'AbsTol',1e-12*scaleV)
verifyEqual(testCase,ur,zeros(size(ur)),'AbsTol',1e-12*scaleU)
% long-term velocity is G Vpl at every time
vLong = reshape(G*sys.Vpl(:),3,[])';
verifyEqual(testCase,squeeze(v(end,:,:)),vLong,'RelTol',1e-12)
end

function testInterseismicBackslipIsLandward(testCase)
% relative to the long-term motion, a station above a locked patch moves
% landward (+north, away from the trench), as for back slip
sys = meshSystem(testCase,true);
Nc = sys.Nc;
y = unicycle.ode.steadyState(sys);
[vr,~] = unicycle.ode.surfaceDeformation(testCase.TestData.G,sys,0,y',struct('removeLongTerm',true));
verifyGreaterThan(testCase,vr(1,1,2),0)
verifyEqual(testCase,size(vr),[1,size(testCase.TestData.xObs,1),3])
verifyEqual(testCase,numel(y),4*Nc)
end

function testDisplacementIsTimeIntegralOfVelocity(testCase)
% postseismic transient after an imposed event: the displacement from the
% slip state matches the trapezoidal integral of the surface velocity, with
% and without the long-term term
yr = testCase.TestData.yr;
sys = meshSystem(testCase,true);
event = unicycle.ode.periodicEvent(sys,50*yr);
y0 = unicycle.ode.applyEvent(sys,unicycle.ode.steadyState(sys),event.dtau,event.slipC);
% the coseismic slip rate is ~1e8 Vpl and decays within seconds, so the
% quadrature grid starts at 1e-12 yr
t = [0,logspace(-12,1,6000)*yr]';
[t,Y] = ode45(@(t,y) unicycle.ode.odeRateDependent(t,y,sys),t,y0,odeset('RelTol',1e-10,'AbsTol',1e-10));
G = testCase.TestData.G;
for removeLongTerm = [false,true]
    o = struct('removeLongTerm',removeLongTerm,'reference','first');
    [v,u] = unicycle.ode.surfaceDeformation(G,sys,t,Y,o);
    uInt = cumtrapz(t,v,1);
    scale = max(abs(u(end,:,:)),[],'all');
    verifyEqual(testCase,u(end,:,:),uInt(end,:,:),'AbsTol',1e-3*scale)
    verifyEqual(testCase,u(1,:,:),zeros(size(u(1,:,:))))
end
end

function testCoseismicStepFromPinnedSlip(testCase)
% for a cycle that starts at the event, reference 'zero' includes the
% coseismic displacement G*slip of pinned and creeping rupture patches
yr = testCase.TestData.yr;
sys = meshSystem(testCase,true);
N = sys.N; Nc = sys.Nc;
weightC = zeros(Nc,1);
weightC(1:3) = 0.5;
event = unicycle.ode.periodicEvent(sys,50*yr,weightC);
yPlus = unicycle.ode.applyEvent(sys,sys.y0 .* [zeros(2*Nc,1);ones(2*Nc,1)],event.dtau,event.slipC);
G = testCase.TestData.G;
[~,u] = unicycle.ode.surfaceDeformation(G,sys,0,yPlus',struct('pinnedSlip',event.slipP));
slipAll = zeros(N,2);
slipAll(sys.iP,:) = event.slipP;
slipAll(sys.iC,:) = event.slipC;
verifyEqual(testCase,squeeze(u(1,:,:)),reshape(G*slipAll(:),3,[])','AbsTol',1e-12*max(abs(G*slipAll(:))))
end
