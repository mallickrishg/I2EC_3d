function rcv = setTestProperties(rcv)
% SETTESTPROPERTIES writes a standard set of user-provided properties
% (m, s, MPa) into a triangleReceiver for tests

yr = 3.1536e7;
rcv.a = 0.015;
rcv.b = 0.010;
rcv.sigma = 50;
rcv.mu0 = 0.6;
rcv.Vo = 1e-6;
rcv.Vs = 3e3;
rcv.Vpl = 0.05/yr;
rcv.Vrake = 70;
rcv.pinnedPosition = [];

end
