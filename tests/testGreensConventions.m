function tests = testGreensConventions
% sign and axis conventions of triangleReceiver and the nikkhoo15 kernels
tests = functiontests(localfunctions);
end

function setupOnce(testCase)
addpath(fileparts(fileparts(mfilename('fullpath'))));
testCase.TestData.dir = tempname;
mkdir(testCase.TestData.dir);
base = makeTestMesh(testCase.TestData.dir,'plane');
testCase.TestData.rcv = unicycle.geometry.triangleReceiver(base,unicycle.greens.nikkhoo15(30e3,0.25));
end

function teardownOnce(testCase)
rmdir(testCase.TestData.dir,'s');
end

function testAxesAfterReading(testCase)
% .ned is [north, east, depth]; rcv.x is [east, north, up]
rcv = testCase.TestData.rcv;
ned = load([fullfile(testCase.TestData.dir,'plane') '.ned']);
verifyEqual(testCase,rcv.x,[ned(:,3),ned(:,2),-ned(:,4)])
% trench along east, dipping toward north
verifyEqual(testCase,[min(rcv.x(:,1)),max(rcv.x(:,1))],[-30e3,30e3],'AbsTol',1e-6)
verifyGreaterThanOrEqual(testCase,rcv.x(:,2),0)
verifyLessThanOrEqual(testCase,rcv.x(:,3),0)
end

function testUnitVectors(testCase)
rcv = testCase.TestData.rcv;
% upward normal, dip of 20 degrees (vertices are written to the mm, so
% the dip is exact to about 1e-5 degrees)
verifyGreaterThan(testCase,rcv.nv(:,3),0)
verifyEqual(testCase,acosd(rcv.nv(:,3)),20*ones(rcv.N,1),'AbsTol',1e-5)
% dip vector points up-dip (toward the trench, -north, +up)
verifyLessThan(testCase,rcv.dv(:,2),0)
verifyGreaterThan(testCase,rcv.dv(:,3),0)
% right-handed: fault dips to the right of strike
verifyEqual(testCase,cross(rcv.sv,rcv.dv,2),rcv.nv,'AbsTol',1e-12)
end

function testPositiveDipSlipIsThrust(testCase)
rcv = testCase.TestData.rcv;
% surface points above the hanging wall and the footwall
obs = [0,10e3,0;0,-10e3,0];
[~,Gd] = rcv.displacementKernels(obs,3);
u = reshape(Gd*ones(rcv.N,1),3,[]);
% hanging wall moves up-dip (toward -north) and up
verifyLessThan(testCase,u(2,1),0)
verifyGreaterThan(testCase,u(3,1),0)
% footwall moves the opposite way horizontally
verifyGreaterThan(testCase,u(2,2),0)
end

function testPositiveStrikeSlipIsLeftLateral(testCase)
rcv = testCase.TestData.rcv;
obs = [0,10e3,0;0,-10e3,0];
[Gs,~] = rcv.displacementKernels(obs,3);
u = reshape(Gs*ones(rcv.N,1),3,[]);
% strike vector is -east here; hanging wall (north) moves west and
% footwall moves east: left-lateral
verifyEqual(testCase,rcv.sv(1,:),[-1,0,0],'AbsTol',1e-12)
verifyLessThan(testCase,u(1,1),0)
verifyGreaterThan(testCase,u(1,2),0)
end

function testSelfStressReducesTraction(testCase)
% slip on a patch reduces shear traction on itself in the slip direction
rcv = testCase.TestData.rcv;
[Kss,~,~,Kdd] = rcv.tractionKernels(rcv);
verifyLessThan(testCase,diag(Kss),0)
verifyLessThan(testCase,diag(Kdd),0)
% uniform thrust on the whole fault drops dip traction on average
verifyLessThan(testCase,mean(Kdd*ones(rcv.N,1)),0)
end
