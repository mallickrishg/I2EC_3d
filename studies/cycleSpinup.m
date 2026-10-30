% Spin up periodic imposed earthquake cycles on the rectangular asperity
% until the cycle is periodic and the slip budget closes
%
% generic_sz15_dx20 mesh, 200 km x 100 km locked asperity, asinh law,
% a-b = 0.005, v0 = Vpl, plate rake 68 deg, Trec = 150 yr. The asperity
% releases exactly Vpl Trec in every event (unicycle.ode.periodicEvent).
% Spin-up starts from v = Vpl at the plate rake everywhere.
%
% parameters that can be set before running (defaults in brackets):
%   f0            - friction coefficient at plate rate [0.6]
%   ruptureExtend - along-strike extension of the rupture beyond the
%                   asperity on each side (m) [0]
%   dtauTarget    - traction change (MPa) on every creeping patch inside
%                   the extension [-3], see unicycle.ode.ruptureSlipForStressChange
%   maxCycles     - [2000]
%   startFrom     - initial pre-event state: 'steady' (BEM steady state with
%                   the asperity locked forever) or 'plate' (v = Vpl at the
%                   plate rake everywhere) ['plate']
% e.g.  f0 = 0.1; ruptureExtend = 20e3; run('studies/cycleSpinup.m')

clearvars -except f0 ruptureExtend dtauTarget maxCycles startFrom
if ~exist('f0','var'), f0 = 0.6; end
if ~exist('ruptureExtend','var'), ruptureExtend = 0; end
if ~exist('dtauTarget','var'), dtauTarget = -3; end
if ~exist('maxCycles','var'), maxCycles = 2000; end
if ~exist('startFrom','var'), startFrom = 'plate'; end
yr = 3.1536e7;
projectDir = fileparts(fileparts(mfilename('fullpath')));
cd(projectDir)
addpath(fullfile(projectDir,'studies'))
outDir = fullfile('studies','output');
if ~exist(outDir,'dir'), mkdir(outDir), end

meshName = 'generic_sz15_dx20';
Trec = 150*yr;

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
dip = 15;
aspX = [-110e3,90e3];
aspW = [40e3,140e3];
rcv.pinnedPosition = rectangularAsperity(rcv,aspX,aspW,dip);

sys = unicycle.ode.assembleODE(rcv,bigK,'asinh');
Nc = sys.Nc;
iC = sys.iC;
Vpl = rcv.Vpl;
[ySS,vSS] = unicycle.ode.steadyState(sys);

% rupture: the asperity, optionally extended along strike onto creeping
% patches with a prescribed stress drop
tag = sprintf('_f0_%g_ext%gkm',f0,ruptureExtend/1e3);
weightC = zeros(Nc,1);
ruptX = aspX + ruptureExtend*[-1,1];
if ruptureExtend > 0
    inRupture = false(rcv.N,1);
    inRupture(rectangularAsperity(rcv,ruptX,aspW,dip)) = true;
    weightC = unicycle.ode.ruptureSlipForStressChange(sys,Trec,inRupture(iC),dtauTarget,...
        unicycle.ode.stateTraction(sys,ySS));
    tag = sprintf('%s_drop%g',tag,-dtauTarget);
end
event = unicycle.ode.periodicEvent(sys,Trec,weightC);
switch startFrom
    case 'steady', y0 = ySS;
    case 'plate', y0 = sys.y0;
end
if ~strcmp(startFrom,'plate'), tag = [tag,'_start',startFrom]; end
fprintf('f0 = %g, rupture %g km x %g km, Trec = %g yr, asperity slip %.2f m, start from %s\n',...
    f0,diff(ruptX)/1e3,diff(aspW)/1e3,Trec/yr,max(hypot(event.slipP(:,1),event.slipP(:,2))),startFrom)

%% spin-up
wall = tic;
out = unicycle.ode.runCycles(sys,event,struct('maxCycles',maxCycles,'verbose',100,'y0',y0));
fprintf('%d cycles (%g kyr) in %.0f s, converged = %d\n',out.nCycles,out.nCycles*Trec/yr/1e3,toc(wall),out.converged)

% slip budget of the stored cycle
VplC = sys.Vpl(iC,:)*Trec;
budget = hypot(out.slipCycle(:,1),out.slipCycle(:,2))./hypot(VplC(:,1),VplC(:,2));
fprintf('stored cycle: slip / (Vpl Trec) on creeping patches in [%.5f, %.5f]\n',min(budget),max(budget))

% periodic pre-event state against the interseismic steady state
yP = out.yStart(:,end);
vPre = sys.v0.*exp(yP(2*Nc+1:3*Nc));
dRakePre = rad2deg(angle(exp(1i*(yP(3*Nc+1:end) - ySS(3*Nc+1:end)))));
vssMag = hypot(vSS(:,1),vSS(:,2));
fprintf('periodic pre-event state vs steady state: |v|/|v_ss| in [%.3f, %.3f], rake diff in [%.2f, %.2f] deg\n',...
    min(vPre./vssMag),max(vPre./vssMag),min(dRakePre),max(dRakePre))
H = out.history;
fprintf('over the spin-up: min |T+| %.3g MPa, max reversed patches %d\n',min(H.minTplus),max(H.nReversed))

zetaS = out.yStart(2*Nc+1:3*Nc,:);
thetaS = out.yStart(3*Nc+1:end,:);
save(fullfile(outDir,['cycleSpinup',tag,'.mat']),'H','zetaS','thetaS','vSS','ySS','yr','Trec',...
    'f0','ruptureExtend','dtauTarget','event','iC')
tCyc = out.t; yCyc = out.y; slipCycle = out.slipCycle; converged = out.converged; nCycles = out.nCycles;
save(fullfile(outDir,['cycleSpinup',tag,'.mat']),'tCyc','yCyc','slipCycle','converged','nCycles','-append')

%% figure
xc = rcv.xc;
isPinned = true(rcv.N,1);
isPinned(iC) = false;
verts = rcv.x/1e3.*[1,1,0];
drawPatches = @(val) patch('Vertices',verts,'Faces',rcv.vertices(iC,:),...
    'FaceVertexCData',val,'FaceColor','flat','EdgeColor',[1,1,1]*0.75);
drawAsperity = @() patch('Vertices',verts,'Faces',rcv.vertices(isPinned,:),...
    'FaceColor',[1,1,1]*0.6,'EdgeColor',[1,1,1]*0.4);
ruptBox = @() plot(ruptX([1,2,2,1,1])/1e3,aspW([1,1,2,2,1])*cosd(dip)/1e3,'k--','LineWidth',1.5);
dmap = interp1([0;0.5;1],[0.15,0.3,0.75;1,1,1;0.75,0.15,0.15],linspace(0,1,24)');
mapAxes = @() set(gca,'FontSize',12,'LineWidth',1.2,'DataAspectRatio',[1,1,1]);
probes = [0,25;-150,90;-120,90;0,165];
probeNames = {'up-dip centre','west flank','west edge','down-dip'};
iProbe = zeros(size(probes,1),1);
for k = 1:size(probes,1)
    [~,iProbe(k)] = min(hypot(xc(iC,1)/1e3 - probes(k,1),xc(iC,2)/1e3 - probes(k,2)));
end
pcols = lines(numel(iProbe));
markProbes = @() scatter(xc(iC(iProbe),1)/1e3,xc(iC(iProbe),2)/1e3,60,pcols,'filled','MarkerEdgeColor','k');

figure(1),clf
set(gcf,'Color','w','Position',[50,50,1800,900])

subplot(2,3,1)
semilogy(1:nCycles,H.slipErr,'LineWidth',2), hold on
semilogy(1:nCycles,H.dzeta,'LineWidth',2)
semilogy(1:nCycles,H.drakeDeg,'LineWidth',2)
yline([1e-3],'k:')
grid on, box on
xlabel('cycle'), ylabel('change per cycle')
legend({'max |\Deltas - V_{pl}T| / V_{pl}T','max |\Delta\zeta|','max |\Deltarake| (deg)'},'Location','northeast')
title(sprintf('spin-up: %d cycles (%.1f kyr)',nCycles,nCycles*Trec/yr/1e3))
set(gca,'FontSize',12,'LineWidth',1.2)

subplot(2,3,2)
drawPatches(budget), hold on
drawAsperity(), ruptBox()
clim(1 + 2e-3*[-1,1]), colormap(gca,dmap)
cb = colorbar; cb.Label.String = 'slip per cycle / (V_{pl} T_{rec})';
title('slip budget of the periodic cycle')
xlabel('East (km)'), ylabel('North (km)'), axis tight, box on, mapAxes()

subplot(2,3,3)
drawPatches(dRakePre), hold on
drawAsperity(), ruptBox(), markProbes()
clim(ceil(max(abs(dRakePre)))*[-1,1]), colormap(gca,dmap)
cb = colorbar; cb.Label.String = 'rake^- - rake_{ss} (deg)';
title('pre-event rake: periodic cycle - steady state')
xlabel('East (km)'), ylabel('North (km)'), axis tight, box on, mapAxes()

subplot(2,3,4)
vCyc = sys.v0.*exp(yCyc(:,2*Nc+1:3*Nc))';
for k = 1:numel(iProbe)
    loglog(tCyc(2:end)/yr,vCyc(iProbe(k),2:end)/Vpl,'-','Color',pcols(k,:),'LineWidth',2), hold on
end
yline(1,'k:')
grid on, box on
xlabel('time after event (yr)'), ylabel('|v| / V_{pl}')
legend(probeNames,'Location','southwest')
title('periodic cycle')
set(gca,'FontSize',12,'LineWidth',1.2)

subplot(2,3,5)
for k = 1:numel(iProbe)
    semilogx(tCyc(2:end)/yr,rad2deg(unwrap(yCyc(2:end,3*Nc+iProbe(k)))),'-','Color',pcols(k,:),'LineWidth',2), hold on
end
yline(rcv.Vrake,'k:')
grid on, box on
xlabel('time after event (yr)'), ylabel('rake (deg)')
set(gca,'FontSize',12,'LineWidth',1.2)

subplot(2,3,6)
drawPatches(1 - vPre./Vpl), hold on
drawAsperity(), ruptBox(), markProbes()
clim([0,1]), colormap(gca,flipud(gray(20)))
cb = colorbar; cb.Label.String = '1 - |v|/V_{pl} just before the event';
title('late-interseismic coupling')
xlabel('East (km)'), ylabel('North (km)'), axis tight, box on, mapAxes()

sgtitle(sprintf('periodic cycles, T_{rec} = %g yr, rupture %g km x %g km, f_0 = %g, a-b = 0.005, plate rake %g°',...
    Trec/yr,diff(ruptX)/1e3,diff(aspW)/1e3,f0,rcv.Vrake))
print(fullfile(outDir,['cycleSpinup',tag,'.png']),'-dpng','-r150')
