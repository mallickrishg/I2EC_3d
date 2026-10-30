function tests = testOdeRateDependent
% physics of unicycle.ode.odeRateDependent
%
% single-patch tests use a spring slider with isotropic stiffness k,
% τ̇ = -k (v - Vpl), for which the ODE has closed-form solutions
tests = functiontests(localfunctions);
end

function setupOnce(testCase)
addpath(fileparts(fileparts(mfilename('fullpath'))));
testCase.TestData.yr = 3.1536e7;
end

function testOutputIsSlipVelocity(testCase)
sys = oneSliderSystem(1,1e-9,deg2rad(30));
y = [1;2;log(3e-9/sys.v0);deg2rad(50)];
dydt = unicycle.ode.odeRateDependent(0,y,sys);
verifyEqual(testCase,dydt(1:2),3e-9*[cosd(50);sind(50)],'RelTol',1e-12)
end

function testPlateRateIsFixedPoint(testCase)
% zero to round-off relative to the natural rate scale k Vpl/((a-b)σ)
for law = {'log','asinh'}
    sys = oneSliderSystem(1,1e-9,deg2rad(30),law{1});
    dydt = unicycle.ode.odeRateDependent(0,sys.y0,sys);
    scale = 1*1e-9/sys.abSigma;
    verifyEqual(testCase,dydt(3:4),[0;0],'AbsTol',1e-14*scale)
end
end

function testProjectionOntoSlipDirection(testCase)
% ζ̇ = τ̇·v̂ / (d|τ|/dζ + η v),  θ̇ = τ̇·v̂⊥ / (|τ| + η v)
eta = 5;
sys = oneSliderSystem(0.7,1e-9,deg2rad(30),'log',0.6,eta);
v = 4e-9; th = deg2rad(80);
y = [0;0;log(v/sys.v0);th];
dydt = unicycle.ode.odeRateDependent(0,y,sys);
taudot = sys.Kcc*v*[cos(th);sin(th)] + sys.taudotLoad;
par = [cos(th),sin(th)]*taudot;
perp = [-sin(th),cos(th)]*taudot;
tau = sys.sigma*(sys.f0 + sys.ab*log(v/sys.v0));
verifyEqual(testCase,dydt(3),par/(sys.abSigma + eta*v),'RelTol',1e-12)
verifyEqual(testCase,dydt(4),perp/(tau + eta*v),'RelTol',1e-12)
end

function testAccelerationIndependentOfF0(testCase)
% load along the slip direction: no rotation, and for the log law ζ̇ does
% not depend on f0 or v0. The asinh law departs from this only through
% d|τ|/dζ = (a-b)σ/sqrt(1 + 4 exp(-2(ζ + f0/(a-b)))), which matters when
% v approaches v0 exp(-f0/(a-b))
y = @(sys) [0;0;log(5e-9/sys.v0);deg2rad(30)];
for law = {'log','asinh'}
    s1 = oneSliderSystem(1,1e-9,deg2rad(30),law{1},0.6);
    s2 = oneSliderSystem(1,1e-9,deg2rad(30),law{1},0.1);
    s3 = s2;
    s3.v0 = 1e-3;
    d1 = unicycle.ode.odeRateDependent(0,y(s1),s1);
    d2 = unicycle.ode.odeRateDependent(0,y(s2),s2);
    d3 = unicycle.ode.odeRateDependent(0,y(s3),s3);
    scale = abs(d1(3));
    verifyEqual(testCase,[d1(4),d2(4),d3(4)],[0,0,0],'AbsTol',1e-14*scale)
    verifyEqual(testCase,d2(3),d1(3),'RelTol',1e-12)
    if strcmp(law{1},'log')
        verifyEqual(testCase,d3(3),d1(3),'RelTol',1e-12)
    else
        [~,dt1] = unicycle.ode.frictionLaw(s1,[0,0,1,0]*y(s1));
        [~,dt3] = unicycle.ode.frictionLaw(s3,[0,0,1,0]*y(s3));
        verifyEqual(testCase,d3(3)/d1(3),dt1/dt3,'RelTol',1e-12)
    end
end
end

function testRotationRateScalesWithTraction(testCase)
% load normal to the slip direction: θ̇ ∝ 1/|τ| = 1/(σ (f0 + (a-b) ζ))
v = 1e-9;
y = @(sys) [0;0;log(v/sys.v0);0];
s1 = oneSliderSystem(1,1e-9,deg2rad(30),'log',0.6);
s2 = oneSliderSystem(1,1e-9,deg2rad(30),'log',0.1);
d1 = unicycle.ode.odeRateDependent(0,y(s1),s1);
d2 = unicycle.ode.odeRateDependent(0,y(s2),s2);
tau1 = unicycle.ode.frictionLaw(s1,log(v/s1.v0));
tau2 = unicycle.ode.frictionLaw(s2,log(v/s2.v0));
verifyEqual(testCase,d2(4)/d1(4),tau1/tau2,'RelTol',1e-12)
% acceleration is the same
verifyEqual(testCase,d2(3),d1(3),'RelTol',1e-12)
end

function testSpringSliderRelaxation(testCase)
% aligned load: (a-b)σ ζ̇ = k (Vpl - v), so v(t) is a logistic curve
%   v(t) = Vpl / (1 + (Vpl/vi - 1) exp(-k Vpl t / ((a-b)σ)))
% independent of f0
yr = testCase.TestData.yr;
k = 1e-3; Vpl = 1e-9; vi = 1e3*Vpl;
t = linspace(0,50,101)'*yr;
opts = odeset('RelTol',1e-10,'AbsTol',1e-12);
for law = {'log','asinh'}
    for f0 = [0.1,0.6]
        sys = oneSliderSystem(k,Vpl,deg2rad(30),law{1},f0);
        y0 = sys.y0;
        y0(3) = log(vi/sys.v0);
        [~,y] = ode45(@(t,y) unicycle.ode.odeRateDependent(t,y,sys),t,y0,opts);
        v = sys.v0*exp(y(:,3));
        vexact = Vpl./(1 + (Vpl/vi - 1)*exp(-k*Vpl*t/sys.abSigma));
        verifyEqual(testCase,v,vexact,'RelTol',1e-6)
        % slip = ∫ v dt along the fixed rake
        slip = sys.abSigma/k*log(1 + vi/Vpl*(exp(k*Vpl*t/sys.abSigma) - 1));
        verifyEqual(testCase,hypot(y(:,1),y(:,2)),slip,'RelTol',1e-6,'AbsTol',1e-9)
        verifyEqual(testCase,y(:,4),deg2rad(30)*ones(size(t)),'AbsTol',1e-12)
    end
end
end

function testRakeRelaxationTimeScale(testCase)
% start at |v| = Vpl rotated by a small angle δ from the plate rake:
%   δ(t) = δ0 exp(-k Vpl t / |τ|),  |τ| = σ (f0 + (a-b) log(Vpl/v0))
% so rotation is f0/(a-b) + log(Vpl/v0) times slower than relaxation of |v|
yr = testCase.TestData.yr;
k = 1e-3; Vpl = 1e-9; d0 = 1e-4;
opts = odeset('RelTol',1e-10,'AbsTol',1e-14);
for f0 = [0.1,0.6]
    sys = oneSliderSystem(k,Vpl,deg2rad(30),'log',f0);
    tau = unicycle.ode.frictionLaw(sys,log(Vpl/sys.v0));
    T = tau/(k*Vpl);
    t = linspace(0,3*T,31)';
    y0 = sys.y0;
    y0(4) = y0(4) + d0;
    [~,y] = ode45(@(t,y) unicycle.ode.odeRateDependent(t,y,sys),t,y0,opts);
    verifyEqual(testCase,y(:,4) - deg2rad(30),d0*exp(-t/T),'AbsTol',1e-4*d0)
    verifyLessThan(testCase,T/yr,1e6)
end
end

function testRadiationDamping(testCase)
% after a large stress step v is seismic, η|v| dominates (a-b)σ and the
% deceleration is slowed by the factor (a-b)σ/((a-b)σ + η v)
yr = testCase.TestData.yr;
eta = 30e3/(2*3e3);
k = 0.1; Vpl = 1e-9;
sys = oneSliderSystem(k,Vpl,0,'asinh',0.6,eta);
% traction 8 MPa above steady state, with friction alone carrying it
% (undamped inverse), v ~ Vpl exp(8/0.25) ~ 80 m/s
y0 = sys.y0;
y0(3) = unicycle.ode.frictionInverse(sys,unicycle.ode.frictionLaw(sys,y0(3)) + 8,false);
v = sys.v0*exp(y0(3));
verifyGreaterThan(testCase,v,1)
sysNoDamp = sys;
sysNoDamp.eta = 0;
dydt = unicycle.ode.odeRateDependent(0,y0,sys);
dydt0 = unicycle.ode.odeRateDependent(0,y0,sysNoDamp);
[~,dtau] = unicycle.ode.frictionLaw(sys,y0(3));
verifyEqual(testCase,dydt(3)/dydt0(3),dtau/(dtau + eta*v),'RelTol',1e-12)
% the solution still relaxes to the plate rate; relaxation time
% (a-b)σ/(k Vpl) is about 80 yr
opts = odeset('RelTol',1e-8,'AbsTol',1e-10);
[~,y] = ode45(@(t,y) unicycle.ode.odeRateDependent(t,y,sys),[0 2e3*yr],y0,opts);
verifyEqual(testCase,sys.v0*exp(y(end,3)),Vpl,'RelTol',1e-6)
end

function testMeshSteadyStateMatchesBEM(testCase)
% on a small mesh with one pinned patch, long integration converges to the
% linear BEM solution K_cc v_c = -τ̇_load (the friction law drops out)
d = tempname;
mkdir(d);
cleanup = onCleanup(@() rmdir(d,'s'));
rcv = unicycle.geometry.triangleReceiver(makeTestMesh(d,'plane'),unicycle.greens.nikkhoo15(30e3,0.25));
rcv = setTestProperties(rcv);
% pin a patch near the centre
[~,ic] = min(hypot(rcv.xc(:,1),rcv.xc(:,2)-15e3));
rcv.pinnedPosition = ic;
bigK = unicycle.greens.computeKernels(rcv,'plane',fullfile(d,'kernels'));
% small f0 keeps the rake relaxation time short
rcv.mu0 = 0.01;
sys = unicycle.ode.assembleODE(rcv,bigK,'asinh');
Nc = sys.Nc;
vbem = sys.Kcc \ (-sys.taudotLoad);
opts = odeset('RelTol',1e-9,'AbsTol',1e-9);
[~,y] = ode45(@(t,y) unicycle.ode.odeRateDependent(t,y,sys),[0 1e4*testCase.TestData.yr],sys.y0,opts);
v = sys.v0.*exp(y(end,2*Nc+1:3*Nc)');
th = y(end,3*Nc+1:end)';
verifyEqual(testCase,[v.*cos(th);v.*sin(th)],vbem,'AbsTol',1e-3*rcv.Vpl)
end

function testDampedTractionBalanceIsConserved(testCase)
% start from the damped coseismic state, |T0| = F(v) + η v, and integrate
% a spring slider with aligned loading. The elastic traction is then
%   |T(t)| = |T0| + k (Vpl t - slip(t))
% and the ODE must keep F(v(t)) + η v(t) = |T(t)| at all times
yr = testCase.TestData.yr;
eta = 30e3/(2*3e3);
k = 0.1; Vpl = 1e-9;
opts = odeset('RelTol',1e-10,'AbsTol',1e-12);
for law = {'log','asinh'}
    sys = oneSliderSystem(k,Vpl,0,law{1},0.6,eta);
    T0 = unicycle.ode.frictionLaw(sys,sys.y0(3)) + 8;
    y0 = sys.y0;
    y0(3) = unicycle.ode.frictionInverse(sys,T0);
    t = [0,logspace(-6,1,50)*yr]';
    [~,y] = ode45(@(t,y) unicycle.ode.odeRateDependent(t,y,sys),t,y0,opts);
    v = sys.v0*exp(y(:,3));
    T = T0 + k*(Vpl*t - y(:,1));
    verifyEqual(testCase,unicycle.ode.frictionLaw(sys,y(:,3)) + eta*v,T,'AbsTol',1e-6)
end
end
