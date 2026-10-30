function tests = testComputeKernels
% caching and geometry checks of unicycle.greens.computeKernels
tests = functiontests(localfunctions);
end

function setupOnce(testCase)
addpath(fileparts(fileparts(mfilename('fullpath'))));
testCase.TestData.dir = tempname;
mkdir(testCase.TestData.dir);
testCase.TestData.base = makeTestMesh(testCase.TestData.dir,'plane');
end

function teardownOnce(testCase)
rmdir(testCase.TestData.dir,'s');
end

function setup(testCase)
% fresh cache folder and mesh for each test
testCase.TestData.kdir = tempname(testCase.TestData.dir);
testCase.TestData.rcv = unicycle.geometry.triangleReceiver(testCase.TestData.base,...
    unicycle.greens.nikkhoo15(30e3,0.25));
end

function testBlockOrdering(testCase)
rcv = testCase.TestData.rcv;
N = rcv.N;
[bigK,bigKn] = unicycle.greens.computeKernels(rcv,'plane',testCase.TestData.kdir);
[Kss,Kds,Ksd,Kdd,Ksn,Kdn] = rcv.tractionKernels(rcv);
verifyEqual(testCase,size(bigK),[2*N,2*N])
% [τs;τd] = [Kss,Kds;Ksd,Kdd]*[s_s;s_d]
verifyEqual(testCase,bigK,[Kss,Kds;Ksd,Kdd])
verifyEqual(testCase,bigKn,[Ksn,Kdn])
end

function testReloadIsIdentical(testCase)
rcv = testCase.TestData.rcv;
[K1,~,f1] = unicycle.greens.computeKernels(rcv,'plane',testCase.TestData.kdir);
[K2,~,f2] = unicycle.greens.computeKernels(rcv,'plane',testCase.TestData.kdir,rcv.N);
verifyEqual(testCase,f2,f1)
verifyEqual(testCase,K2,K1)
end

function testElasticModelChangesFile(testCase)
rcv = testCase.TestData.rcv;
[~,~,f1] = unicycle.greens.computeKernels(rcv,'plane',testCase.TestData.kdir);
rcv2 = unicycle.geometry.triangleReceiver(testCase.TestData.base,unicycle.greens.nikkhoo15(40e3,0.25));
[K2,~,f2] = unicycle.greens.computeKernels(rcv2,'plane',testCase.TestData.kdir);
rcv3 = unicycle.geometry.triangleReceiver(testCase.TestData.base,unicycle.greens.nikkhoo15(30e3,0.3));
[~,~,f3] = unicycle.greens.computeKernels(rcv3,'plane',testCase.TestData.kdir);
verifyNotEqual(testCase,f2,f1)
verifyNotEqual(testCase,f3,f1)
% kernels scale with rigidity
[K1] = unicycle.greens.computeKernels(rcv,'plane',testCase.TestData.kdir);
verifyEqual(testCase,K2,K1*4/3,'RelTol',1e-10)
end

function testGeometryChangesFile(testCase)
rcv = testCase.TestData.rcv;
[~,~,f1] = unicycle.greens.computeKernels(rcv,'plane',testCase.TestData.kdir);
rcv.x(1,3) = rcv.x(1,3) - 1;
[~,~,f2] = unicycle.greens.computeKernels(rcv,'plane',testCase.TestData.kdir);
verifyNotEqual(testCase,f2,f1)
end

function testTamperedCacheIsDetected(testCase)
rcv = testCase.TestData.rcv;
[~,~,f1] = unicycle.greens.computeKernels(rcv,'plane',testCase.TestData.kdir);
m = matfile(f1,'Writable',true);
Kdd = m.Kdd;
Kdd(:,1) = 2*Kdd(:,1);
m.Kdd = Kdd;
% check every column so the tampered one is hit
verifyError(testCase,@() unicycle.greens.computeKernels(rcv,'plane',testCase.TestData.kdir,rcv.N),?MException)
end

function testForeignCacheIsDetected(testCase)
rcv = testCase.TestData.rcv;
[~,~,f1] = unicycle.greens.computeKernels(rcv,'plane',testCase.TestData.kdir);
rcv2 = unicycle.geometry.triangleReceiver(testCase.TestData.base,unicycle.greens.nikkhoo15(40e3,0.25));
[~,~,f2] = unicycle.greens.computeKernels(rcv2,'plane',testCase.TestData.kdir);
% another mesh's kernels saved under this mesh's file name
copyfile(f2,f1);
verifyError(testCase,@() unicycle.greens.computeKernels(rcv,'plane',testCase.TestData.kdir,0),?MException)
end
