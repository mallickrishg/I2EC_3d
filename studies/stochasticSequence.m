% Stochastic imposed earthquake sequence on a fully velocity-strengthening
% fault (no locked zone)
%
% generic_sz15_dx20 mesh, asinh law, a-b = 0.005, sigma = 50 MPa, v0 = Vpl,
% plate rake 68 deg, Vpl = 54 mm/yr. Every 150 yr cycle has three events,
% one per segment (west, centre, east) at 0, 20 and 30 yr into the cycle.
% Each event is a square of uniform slip along the plate rake. Per event,
% the square's centre is jittered around the segment's nominal centre, its
% edge length is Gaussian and its slip is Gaussian about Vpl x 150 yr =
% 8.1 m (std 1.5 m), so on average over many cycles the events alone close
% the slip budget at the nominal segment centres. The events
% are imposed on creeping patches (unicycle.ode.sequenceEvent with no
% pinned patches) and the creep in between evolves with RSF. The whole
% catalogue is one sequence of length nCycles x 150 yr, run once
% (unicycle.ode.spinUp with repeat = false) from v = Vpl at the plate rake.
% Patches near the jittered square edges are not ruptured every cycle, so
% creep still supplies part of their budget.
%
% parameters that can be set before running (defaults in brackets):
%   seed     - random seed; 0 = deterministic reference with the mean event
%              of every segment (nominal centre, mean edge, mean slip) [1]
%   nCycles  - number of 150 yr cycles [300]
%   f0       - [0.3]
% e.g.  seed = 0; run('studies/stochasticSequence.m')
% figures: studies/stochasticSequenceCompare.m

clearvars -except seed nCycles f0
if ~exist('seed','var'), seed = 1; end
if ~exist('nCycles','var'), nCycles = 300; end
if ~exist('f0','var'), f0 = 0.3; end
yr = 3.1536e7;
projectDir = fileparts(fileparts(mfilename('fullpath')));
cd(projectDir)
addpath(fullfile(projectDir,'studies'))
outDir = fullfile('studies','output');
if ~exist(outDir,'dir'), mkdir(outDir), end

meshName = 'generic_sz15_dx20';
Tcyc = 150*yr;
dip = 15;

% sequence design (m, yr)
seg.name = {'west','centre','east'};
seg.x0 = [-150e3,0,150e3];        % nominal centre along strike
seg.w0 = 90e3*[1,1,1];            % nominal centre down-dip distance
seg.tCycle = [0,20,30];           % time of the event within the cycle
jitterX = 25e3;                   % std of the centre along strike
jitterW = 15e3;                   % std of the centre down dip
Lmean = 100e3; Lstd = 25e3;       % edge length, truncated to Llim
Llim = [40e3,160e3];
slipMean = 0.054*150;             % Gaussian slip (m): Vpl x cycle length
slipStd = 1.5;                    % truncated to positive values

eM = unicycle.greens.nikkhoo15(30e3,0.25);
rcv = unicycle.geometry.triangleReceiver(fullfile('fault',meshName),eM);
bigK = unicycle.greens.computeKernels(rcv,meshName);

% user-provided properties (m, s, MPa)
rcv.a = 0.015;
rcv.b = 0.010;
rcv.sigma = 50;
rcv.mu0 = f0;
rcv.Vs = 3e3;
rcv.Vpl = 0.054/yr;
rcv.Vo = rcv.Vpl;
rcv.Vrake = 68;
rcv.pinnedPosition = [];

sys = unicycle.ode.assembleODE(rcv,bigK,'asinh');
Nc = sys.Nc;
assert(isequal(sys.iC(:),(1:rcv.N)'),'no patch may be pinned')
rakeVec = [cosd(rcv.Vrake),sind(rcv.Vrake)];

%% catalogue
nSeg = numel(seg.x0);
nEv = nCycles*nSeg;
catalog = table(zeros(nEv,1),zeros(nEv,1),zeros(nEv,1),zeros(nEv,1),zeros(nEv,1),zeros(nEv,1),...
    zeros(nEv,1),zeros(nEv,1),'VariableNames',{'cycle','segment','t','xc','wc','L','slip','nPatches'});
if seed == 0
    tag = sprintf('_f0_%g_mean',f0);
else
    tag = sprintf('_f0_%g_seed%d',f0,seed);
    rng(seed,'twister')
end
k = 0;
for n = 1:nCycles
    for j = 1:nSeg
        k = k + 1;
        if seed == 0
            x = seg.x0(j); w = seg.w0(j); L = Lmean; s = slipMean;
        else
            x = seg.x0(j) + jitterX*randn;
            w = seg.w0(j) + jitterW*randn;
            L = Lmean + Lstd*randn;
            while L < Llim(1) || L > Llim(2), L = Lmean + Lstd*randn; end
            s = slipMean + slipStd*randn;
            while s <= 0, s = slipMean + slipStd*randn; end
        end
        catalog{k,:} = [n,j,((n-1)*Tcyc + seg.tCycle(j)*yr)/yr,x,w,L,s,0];
    end
end

events = repmat(struct('t',0,'dtau',zeros(2*Nc,1),'slipP',zeros(0,2),'slipC',zeros(Nc,2)),nEv,1);
inSquare = false(Nc,nEv);
for k = 1:nEv
    in = rectangularAsperity(rcv,catalog.xc(k) + catalog.L(k)/2*[-1,1],catalog.wc(k) + catalog.L(k)/2*[-1,1],dip);
    inSquare(in,k) = true;
    catalog.nPatches(k) = numel(in);
    slip = zeros(rcv.N,2);
    slip(in,:) = catalog.slip(k)*repmat(rakeVec,numel(in),1);
    events(k) = unicycle.ode.sequenceEvent(sys,catalog.t(k)*yr,slip);
end
fprintf('%s: %d cycles, %d events, f0 = %g; mean slip %.2f m (Vpl Tcyc = %.2f m), mean edge %.0f km, %d to %d patches\n',...
    tag(2:end),nCycles,nEv,f0,mean(catalog.slip),rcv.Vpl*Tcyc,mean(catalog.L)/1e3,min(catalog.nPatches),max(catalog.nPatches))

%% run
% one pass through the catalogue, dense output after the three events of
% the last cycle
seq = struct('T',nCycles*Tcyc,'events',events);
wall = tic;
[yEnd,out] = unicycle.ode.spinUp(sys,seq,struct('repeat',false,'dense',nEv-nSeg+1:nEv,...
    'nOut',81,'verbose',90));
fprintf('done in %.0f s (%.2f s per cycle)\n',toc(wall),toc(wall)/nCycles)

% state at the end of every cycle (just before the next west event); the
% last one is the end of the run
iEnd = find(catalog.segment == 1 & catalog.cycle > 1);
zetaEnd = [out.zetaPre(:,iEnd),yEnd(2*Nc+1:3*Nc)];
thetaEnd = [out.thetaPre(:,iEnd),yEnd(3*Nc+1:end)];
slipEnd = cat(3,out.slipPre(:,:,iEnd),out.slipSince);
slipEndPar = squeeze(slipEnd(:,1,:)*rakeVec(1) + slipEnd(:,2,:)*rakeVec(2));   % Nc x nCycles

% coseismic slip along the rake received by each patch in each cycle
seisCycle = zeros(Nc,nCycles);
for k = 1:nEv
    seisCycle(:,catalog.cycle(k)) = seisCycle(:,catalog.cycle(k)) + catalog.slip(k)*inSquare(:,k);
end

coupling = 1 - sys.v0.*exp(zetaEnd)/rcv.Vpl;
Vpl = rcv.Vpl;
xc = rcv.xc;
dense = out.dense;
save(fullfile(outDir,['stochasticSequence',tag,'.mat']),'catalog','seg','jitterX','jitterW','Lmean','Lstd',...
    'Llim','slipMean','slipStd','Tcyc','yr','f0','seed','nCycles','Nc','Vpl','xc','dip','rakeVec','inSquare',...
    'zetaEnd','thetaEnd','slipEndPar','seisCycle','coupling','dense','meshName')
fprintf('saved %s\n',fullfile(outDir,['stochasticSequence',tag,'.mat']))
