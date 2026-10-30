function tests = testFrictionLaw
% steady-state friction laws unicycle.ode.frictionLaw and frictionInverse
tests = functiontests(localfunctions);
end

function setupOnce(testCase)
addpath(fileparts(fileparts(mfilename('fullpath'))));
end

function sys = params(law,f0)
sys.law = law;
sys.sigma = 50;
sys.ab = 0.005;
sys.abSigma = sys.ab*sys.sigma;
sys.f0 = f0;
sys.rho = f0/sys.ab;
sys.v0 = 1e-6;
end

function testLogLaw(testCase)
% |τ| = σ (f0 + (a-b) log(v/v0)) and d|τ|/dζ = (a-b)σ
sys = params('log',0.6);
zeta = linspace(-20,20,41)';
[tau,dtau] = unicycle.ode.frictionLaw(sys,zeta);
verifyEqual(testCase,tau,50*(0.6 + 0.005*zeta),'AbsTol',1e-12)
verifyEqual(testCase,dtau,sys.abSigma*ones(size(zeta)),'AbsTol',1e-15)
end

function testAsinhMatchesLogForRealisticF0(testCase)
% laws are identical when v/(2v0) exp(f0/(a-b)) >> 1
zeta = linspace(-20,20,201)';
[ta,da] = unicycle.ode.frictionLaw(params('asinh',0.6),zeta);
[tl,dl] = unicycle.ode.frictionLaw(params('log',0.6),zeta);
verifyEqual(testCase,ta,tl,'AbsTol',1e-12)
verifyEqual(testCase,da,dl,'AbsTol',1e-12)
end

function testAsinhRegularizesSmallVelocity(testCase)
% with small f0/(a-b) the log law turns negative at low v, asinh does not
zeta = [-30;-10;-2;0;2;10];
ta = unicycle.ode.frictionLaw(params('asinh',0.01),zeta);
tl = unicycle.ode.frictionLaw(params('log',0.01),zeta);
verifyGreaterThan(testCase,ta,0)
verifyLessThan(testCase,tl(1),0)
% and they still agree at high v
verifyEqual(testCase,ta(end),tl(end),'RelTol',1e-6)
end

function testAsinhDerivative(testCase)
for f0 = [0.01,0.6]
    sys = params('asinh',f0);
    zeta = linspace(-10,10,41)';
    [~,da] = unicycle.ode.frictionLaw(sys,zeta);
    h = 1e-5;
    fd = (unicycle.ode.frictionLaw(sys,zeta+h) - unicycle.ode.frictionLaw(sys,zeta-h))/(2*h);
    verifyEqual(testCase,da,fd,'AbsTol',1e-7*sys.abSigma)
end
end

function testInverseRoundTrip(testCase)
zeta = linspace(-20,20,41)';
for law = {'log','asinh'}
    for f0 = [0.01,0.6]
        sys = params(law{1},f0);
        tau = unicycle.ode.frictionLaw(sys,zeta);
        verifyEqual(testCase,unicycle.ode.frictionInverse(sys,tau),zeta,'AbsTol',1e-9)
    end
end
end

function testNoOverflowForLargeF0OverAB(testCase)
% f0/(a-b) = 1200 would overflow exp() in a naive asinh
sys = params('asinh',6);
zeta = linspace(-20,20,41)';
[tau,dtau] = unicycle.ode.frictionLaw(sys,zeta);
verifyTrue(testCase,all(isfinite([tau;dtau])))
verifyEqual(testCase,tau,50*(6 + 0.005*zeta),'RelTol',1e-12)
verifyEqual(testCase,unicycle.ode.frictionInverse(sys,tau),zeta,'AbsTol',1e-9)
end

function testDampedInverse(testCase)
% |τ| = F(v) + η v, solved implicitly for v
eta = 30e3/(2*3e3);
for law = {'log','asinh'}
    for f0 = [0.01,0.6]
        sys = params(law{1},f0);
        sys.eta = eta;
        % from creep to seismic slip: traction steps of 0 to 15 MPa above
        % steady state at 1e-9 m/s
        tau = unicycle.ode.frictionLaw(sys,log(1e-9/sys.v0)) + [0;0.1;1;3;5;8;15];
        zeta = unicycle.ode.frictionInverse(sys,tau);
        v = sys.v0*exp(zeta);
        verifyEqual(testCase,unicycle.ode.frictionLaw(sys,zeta) + eta*v,tau,'RelTol',1e-12)
        % damping lowers v below the undamped value, and below |τ|/η where
        % |τ| > 0 (the log law gives |τ| < 0 for v < v0 exp(-f0/(a-b)))
        zetaUndamped = unicycle.ode.frictionInverse(sys,tau,false);
        verifyLessThanOrEqual(testCase,zeta,zetaUndamped)
        verifyLessThanOrEqual(testCase,v(tau > 0),tau(tau > 0)/eta)
        % at creep rates the shift is ηv/(dF/dζ) to first order
        [~,dF] = unicycle.ode.frictionLaw(sys,zeta(1));
        verifyEqual(testCase,zetaUndamped(1) - zeta(1),eta*v(1)/dF,'RelTol',1e-3)
        % at seismic rates damping dominates
        verifyLessThan(testCase,v(end),1e-2*sys.v0*exp(zetaUndamped(end)))
    end
end
end

function testDampedInverseWithoutEtaIsUndamped(testCase)
for law = {'log','asinh'}
    sys = params(law{1},0.6);
    tau = [20;30;35];
    z1 = unicycle.ode.frictionInverse(sys,tau);
    sys.eta = 0;
    z2 = unicycle.ode.frictionInverse(sys,tau);
    verifyEqual(testCase,z2,z1)
    verifyEqual(testCase,unicycle.ode.frictionLaw(sys,z1),tau,'RelTol',1e-12)
end
end
