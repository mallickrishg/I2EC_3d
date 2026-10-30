function tests = testKLBasis
% Karhunen-Loève basis for low-rank fields: unicycle.inverse.klBasis
tests = functiontests(localfunctions);
end

function setupOnce(testCase)
addpath(fileparts(fileparts(mfilename('fullpath'))));
testCase.TestData.dir = tempname;
mkdir(testCase.TestData.dir);
eM = unicycle.greens.nikkhoo15(30e3,0.25);
% 100 km x 50 km plane at two resolutions
testCase.TestData.coarse = unicycle.geometry.triangleReceiver(...
    makeTestMesh(testCase.TestData.dir,'coarse',100e3,50e3,20,10e3),eM);
testCase.TestData.fine = unicycle.geometry.triangleReceiver(...
    makeTestMesh(testCase.TestData.dir,'fine',100e3,50e3,20,5e3),eM);
end

function teardownOnce(testCase)
rmdir(testCase.TestData.dir,'s');
end

function testOrthonormalAndReconstructsCovariance(testCase)
% phi' W phi = I, C W phi = lambda phi, and with all modes the covariance
% is reproduced exactly
rcv = testCase.TestData.coarse;
n = rcv.N;
ell = 30e3;
basis = unicycle.inverse.klBasis(rcv,[],struct('ell',ell,'removeMean',false,'rMax',n,'sigma',2));
W = diag(basis.area);
verifyEqual(testCase,basis.phi'*W*basis.phi,eye(n),'AbsTol',1e-10)
xc = rcv.xc;
r = sqrt(max(sum(xc.^2,2) + sum(xc.^2,2)' - 2*(xc*xc'),0));
C = 4*exp(-r/ell);
verifyEqual(testCase,C*W*basis.phi,basis.phi.*basis.lambda','AbsTol',1e-8*max(basis.lambda))
verifyEqual(testCase,basis.Phi*basis.Phi',C,'AbsTol',1e-10*max(C(:)))
verifyEqual(testCase,basis.varFraction(end),1,'AbsTol',1e-12)
verifyLessThanOrEqual(testCase,diff(basis.lambda),1e-12*basis.lambda(1))
end

function testRemoveMeanGivesZeroMeanModes(testCase)
% with removeMean (default) every mode has zero area-weighted mean, so the
% mean level is a separate parameter
rcv = testCase.TestData.coarse;
basis = unicycle.inverse.klBasis(rcv);
verifyTrue(testCase,basis.removeMean)
verifyEqual(testCase,basis.area'*basis.phi,zeros(1,size(basis.phi,2)),'AbsTol',1e-10*sqrt(sum(basis.area)))
verifyEqual(testCase,basis.varFraction(end) <= 1,true)
end

function testDefaultLengthIsHalfTheMeshExtent(testCase)
% L is the largest principal extent of the centroids: for a 100 km x 50 km
% plane that is the along-strike extent of the centroids
rcv = testCase.TestData.coarse;
basis = unicycle.inverse.klBasis(rcv);
xc = rcv.xc;
verifyEqual(testCase,basis.L,max(xc(:,1)) - min(xc(:,1)),'RelTol',1e-10)
verifyEqual(testCase,basis.ell,basis.L/2)
end

function testLongCorrelationGivesSmoothLeadingModes(testCase)
% for ell >> L (here 100 L) the leading mode without centring is nearly
% constant and holds all but ~L/ell of the variance, and
% with centring it is nearly a linear gradient along the longest axis
rcv = testCase.TestData.fine;
xc = rcv.xc;
b1 = unicycle.inverse.klBasis(rcv,[],struct('ell',1e7,'removeMean',false,'rMax',3));
verifyLessThan(testCase,std(b1.phi(:,1))/abs(mean(b1.phi(:,1))),1e-3)
verifyGreaterThan(testCase,b1.varFraction(1),0.99)
b2 = unicycle.inverse.klBasis(rcv,[],struct('ell',1e7,'rMax',3));
cc = corrcoef(b2.phi(:,1),xc(:,1));
verifyGreaterThan(testCase,abs(cc(1,2)),0.99)
end

function testMaternKernels(testCase)
% closed-form nu = 0.5, 1.5, 2.5, Inf agree with the general Bessel form
rcv = testCase.TestData.coarse;
ell = 25e3;
for nu = [0.5,1.5,2.5]
    b = unicycle.inverse.klBasis(rcv,[],struct('ell',ell,'nu',nu,'removeMean',false,'rMax',rcv.N));
    bg = unicycle.inverse.klBasis(rcv,[],struct('ell',ell,'nu',nu + 1e-9,'removeMean',false,'rMax',rcv.N));
    verifyEqual(testCase,bg.Phi*bg.Phi',b.Phi*b.Phi','AbsTol',1e-6)
end
bInf = unicycle.inverse.klBasis(rcv,[],struct('ell',ell,'nu',Inf,'removeMean',false,'rMax',rcv.N));
% large nu tends to the Gaussian kernel, and the general form stays finite
bBig = unicycle.inverse.klBasis(rcv,[],struct('ell',ell,'nu',200,'removeMean',false,'rMax',rcv.N));
verifyTrue(testCase,all(isfinite(bBig.Phi(:))))
verifyEqual(testCase,bBig.Phi*bBig.Phi',bInf.Phi*bInf.Phi','AbsTol',1e-2)
end

function testModesAreMeshIndependent(testCase)
% area weighting makes the leading eigenvalues and modes (nearly)
% independent of the mesh resolution
ell = 40e3;
bc = unicycle.inverse.klBasis(testCase.TestData.coarse,[],struct('ell',ell,'rMax',3));
bf = unicycle.inverse.klBasis(testCase.TestData.fine,[],struct('ell',ell,'rMax',3));
verifyEqual(testCase,bf.lambda,bc.lambda,'RelTol',0.03)
% compare mode 1 at the coarse centroids (nearest fine centroid)
[~,inear] = min(sum((permute(testCase.TestData.coarse.xc,[1,3,2]) - ...
    permute(testCase.TestData.fine.xc,[3,1,2])).^2,3),[],2);
for k = 1:3
    cc = corrcoef(bc.phi(:,k),bf.phi(inear,k));
    verifyGreaterThan(testCase,abs(cc(1,2)),0.99)
end
end

function testSubsetOfPatches(testCase)
% basis on creeping patches only (e.g. lin.iC)
rcv = testCase.TestData.coarse;
iC = (1:2:rcv.N)';
basis = unicycle.inverse.klBasis(rcv,iC,struct('rMax',4));
verifyEqual(testCase,size(basis.Phi),[numel(iC),4])
verifyEqual(testCase,basis.xc,rcv.xc(iC,:))
end
