function tests = testAssembleODE
% setup of the interseismic system by unicycle.ode.assembleODE
tests = functiontests(localfunctions);
end

function setupOnce(testCase)
addpath(fileparts(fileparts(mfilename('fullpath'))));
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
testCase.TestData.rcv = setTestProperties(rcv);
end

function testScalarExpansionAndPartition(testCase)
rcv = testCase.TestData.rcv;
rcv.pinnedPosition = [3,7,8];
sys = unicycle.ode.assembleODE(rcv,testCase.TestData.bigK);
Nc = rcv.N - 3;
verifyEqual(testCase,sys.Nc,Nc)
verifyEqual(testCase,sys.iC,setdiff(1:rcv.N,[3,7,8])')
verifyEqual(testCase,size(sys.Kcc),[2*Nc,2*Nc])
verifyEqual(testCase,size(sys.taudotLoad),[2*Nc,1])
verifyEqual(testCase,sys.abSigma,0.005*50*ones(Nc,1),'AbsTol',1e-15)
verifyEqual(testCase,sys.rho,0.6/0.005*ones(Nc,1),'RelTol',1e-12)
verifyEqual(testCase,sys.eta,30e3/(2*3e3)*ones(Nc,1))
verifyEqual(testCase,sys.law,'asinh')
end

function testPlateRateComponents(testCase)
rcv = testCase.TestData.rcv;
sys = unicycle.ode.assembleODE(rcv,testCase.TestData.bigK);
verifyEqual(testCase,sys.Vpl,repmat(rcv.Vpl*[cosd(70),sind(70)],rcv.N,1),'RelTol',1e-12)
end

function testLoadingWithPinnedPatchesAtZero(testCase)
% τ̇_c from the RHS must equal [K (v - Vpl)]_c with v = 0 on pinned patches
rcv = testCase.TestData.rcv;
bigK = testCase.TestData.bigK;
N = rcv.N;
rcv.pinnedPosition = 10:16;
sys = unicycle.ode.assembleODE(rcv,bigK);
rowsC = [sys.iC;N+sys.iC];
vc = randn(2*sys.Nc,1)*1e-9;
v = zeros(2*N,1);
v(rowsC) = vc;
full = bigK*(v - sys.Vpl(:));
verifyEqual(testCase,sys.Kcc*vc + sys.taudotLoad,full(rowsC),'RelTol',1e-12)
end

function testNoPinningIsSteadyState(testCase)
% with nothing pinned, v = Vpl everywhere gives zero stressing
rcv = testCase.TestData.rcv;
sys = unicycle.ode.assembleODE(rcv,testCase.TestData.bigK);
vc = reshape(sys.Vpl(sys.iC,:),[],1);
verifyEqual(testCase,sys.Kcc*vc + sys.taudotLoad,zeros(2*sys.Nc,1),'AbsTol',1e-12*norm(sys.taudotLoad))
end

function testInitialCondition(testCase)
rcv = testCase.TestData.rcv;
sys = unicycle.ode.assembleODE(rcv,testCase.TestData.bigK);
Nc = sys.Nc;
verifyEqual(testCase,sys.y0(1:2*Nc),zeros(2*Nc,1))
verifyEqual(testCase,sys.v0.*exp(sys.y0(2*Nc+1:3*Nc)),rcv.Vpl*ones(Nc,1),'RelTol',1e-12)
verifyEqual(testCase,sys.y0(3*Nc+1:end),deg2rad(70)*ones(Nc,1),'AbsTol',1e-12)
end

function testValidation(testCase)
bigK = testCase.TestData.bigK;
rcv = testCase.TestData.rcv;
verifyError(testCase,@() unicycle.ode.assembleODE(rcv,bigK,'linear'),?MException)
verifyError(testCase,@() unicycle.ode.assembleODE(rcv,bigK(1:end-2,1:end-2)),?MException)
rcv.b = 0.02;
verifyError(testCase,@() unicycle.ode.assembleODE(rcv,bigK),?MException)
% velocity weakening is allowed on pinned patches only
rcv.b = 0.01*ones(rcv.N,1);
rcv.b(5) = 0.02;
verifyError(testCase,@() unicycle.ode.assembleODE(rcv,bigK),?MException)
rcv.pinnedPosition = 5;
verifyWarningFree(testCase,@() unicycle.ode.assembleODE(rcv,bigK))
rcv.Vpl = 0;
verifyError(testCase,@() unicycle.ode.assembleODE(rcv,bigK),?MException)
rcv.Vpl = 1e-9;
rcv.a = [0.015;0.015];
verifyError(testCase,@() unicycle.ode.assembleODE(rcv,bigK),?MException)
end

function testWarnsOnDefaultShearWaveSpeed(testCase)
rcv = unicycle.geometry.triangleReceiver(testCase.TestData.base,unicycle.greens.nikkhoo15(30e3,0.25));
rcv.Vpl = 1e-9;
verifyWarning(testCase,@() unicycle.ode.assembleODE(rcv,testCase.TestData.bigK),'unicycle:ode:assembleODE:defaultVs')
end
