% One imposed earthquake on the rectangular asperity and the postseismic
% response of the surrounding creeping fault
%
% generic_sz15_dx20 mesh, 200 km x 100 km locked asperity, asinh law,
% a-b = 0.005, v0 = Vpl, plate rake 68 deg. The fault starts from the
% interseismic steady state, the asperity releases the slip deficit of one
% recurrence interval, and the creeping fault is integrated through one
% full interseismic period.
%
% parameters that can be set before running (defaults in brackets):
%   f0            - friction coefficient at plate rate [0.6]
%   ruptureExtend - along-strike extension of the rupture beyond the
%                   asperity on each side (m) [0]. The creeping patches
%                   inside the extension slip by Vpl Trec coseismically,
%                   which gives them a stress drop. 20e3 (one grid column
%                   on each side) makes the rupture 20 % larger in area.
% e.g.  f0 = 0.1; ruptureExtend = 20e3; run('studies/postseismicExample.m')

clearvars -except f0 ruptureExtend
if ~exist('f0','var'), f0 = 0.6; end
if ~exist('ruptureExtend','var'), ruptureExtend = 0; end
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

% rupture: the asperity plus creeping patches in the along-strike extension
ruptX = aspX + ruptureExtend*[-1,1];
inRupture = false(rcv.N,1);
inRupture(rectangularAsperity(rcv,ruptX,aspW,dip)) = true;
weightC = double(inRupture(iC));
tag = sprintf('_f0_%g_ext%gkm',f0,ruptureExtend/1e3);
fprintf('f0 = %g, rupture %g km x %g km (%d pinned + %d creeping patches, %.0f%% of the asperity area)\n',...
    f0,diff(ruptX)/1e3,diff(aspW)/1e3,numel(sys.iP),sum(weightC),...
    100*sum(rcv.area(inRupture))/sum(rcv.area(sys.iP)))

%% event
[ySS,vSS] = unicycle.ode.steadyState(sys);
[dtau,slipP,slipC] = unicycle.ode.eventStressChange(sys,Trec,1,weightC);
y0 = unicycle.ode.applyEvent(sys,ySS,dtau,slipC);
dtauV = [dtau(1:Nc),dtau(Nc+1:end)];
plateDir = [cosd(rcv.Vrake),sind(rcv.Vrake)];
vPlus = sys.v0.*exp(y0(2*Nc+1:3*Nc));
fprintf('coseismic slip %.2f m, Δτ along plate rake in [%.2f, %.2f] MPa, max |Δτ| %.2f MPa\n',...
    max(hypot(slipP(:,1),slipP(:,2))),min(dtauV*plateDir'),max(dtauV*plateDir'),max(hypot(dtauV(:,1),dtauV(:,2))))
Tplus = unicycle.ode.stateTraction(sys,y0);
Tminus = unicycle.ode.stateTraction(sys,ySS);
reversed = sum(Tplus.*Tminus,2) < 0;
fprintf('|T-| = %.1f MPa; |T+| in [%.3g, %.3g] MPa; traction reversed on %d patches\n',...
    median(hypot(Tminus(:,1),Tminus(:,2))),min(hypot(Tplus(:,1),Tplus(:,2))),...
    max(hypot(Tplus(:,1),Tplus(:,2))),sum(reversed))
[~,iMax] = max(vPlus);
fprintf('postseismic v+ up to %.3g m/s = %.3g Vpl (damping cap |T+|/eta there %.3g m/s)\n',...
    vPlus(iMax),vPlus(iMax)/Vpl,norm(Tplus(iMax,:))/sys.eta(iMax))

%% postseismic and interseismic period
tout = [0,logspace(-8,log10(Trec/yr),161)]*yr;
% stop if the integration takes longer than 10 min of wall-clock time
wallStart = tic;
opts = odeset('RelTol',1e-8,'AbsTol',1e-8,...
    'OutputFcn',@(t,y,flag) ~isempty(t) && toc(wallStart) > 600);
lastwarn('');
[t,y] = ode45(@(t,y) unicycle.ode.odeRateDependent(t,y,sys),tout,y0,opts);
[warnMsg,warnId] = lastwarn;
runtime = toc(wallStart);
nSteps = numel(t);
fprintf('integrated to %g yr in %.1f s (%d outputs)\n',t(end)/yr,runtime,nSteps)
if ~isempty(warnMsg), fprintf('ode45 warning [%s]: %s\n',warnId,warnMsg), end
completed = abs(t(end) - Trec) < 1e-6*Trec;

zeta = y(:,2*Nc+1:3*Nc)';
theta = y(:,3*Nc+1:end)';
vMag = sys.v0.*exp(zeta);
slipS = y(:,1:Nc)';
slipD = y(:,Nc+1:2*Nc)';
% afterslip relative to steady creep, excluding imposed coseismic slip
afterS = slipS - slipC(:,1) - vSS(:,1).*t';
afterD = slipD - slipC(:,2) - vSS(:,2).*t';

% accuracy check: the traction implied by the state must equal the elastic
% traction T+ + K_cc (slip(t) - slip+) + τ̇_load t
Tstate = (unicycle.ode.frictionLaw(sys,zeta) + sys.eta.*vMag).*cat(3,cos(theta),sin(theta));
dslip = [slipS - y0(1:Nc);slipD - y0(Nc+1:2*Nc)];
Tel = [Tplus(:,1);Tplus(:,2)] + sys.Kcc*dslip + sys.taudotLoad*t';
tracErr = max(abs([Tstate(:,:,1);Tstate(:,:,2)] - Tel),[],1);
fprintf('max traction mismatch over the run: %.3g MPa = %.3g (a-b)sigma\n',max(tracErr),max(tracErr)/0.25)
% backward slip
vPar = vMag.*cos(theta - deg2rad(rcv.Vrake));
fprintf('patches slipping backward (v . plate dir < 0) at some time: %d, min |T| over the run %.3g MPa\n',...
    sum(any(vPar < 0,2)),min(hypot(Tstate(:,:,1),Tstate(:,:,2)),[],'all'))

% preview of the cycle periodicity check (Task 6): state at the end of the
% cycle compared to the steady state the cycle started from
dRakeEnd = rad2deg(angle(exp(1i*(theta(:,end) - ySS(3*Nc+1:end)))));
fprintf('after one cycle: max |v/v_ss - 1| = %.3g, max |Δrake| = %.2f deg\n',...
    max(abs(vMag(:,end)./hypot(vSS(:,1),vSS(:,2)) - 1)),max(abs(dRakeEnd)))
fprintf('coseismic rotation of the slip-rate direction: max %.2f deg\n',...
    max(abs(rad2deg(angle(exp(1i*(theta(:,1) - ySS(3*Nc+1:end))))))))
fprintf('total afterslip / coseismic slip: max %.2f\n',max(hypot(afterS(:,end),afterD(:,end)))/max(hypot(slipP(:,1),slipP(:,2))))

save(fullfile(outDir,['postseismicExample',tag,'.mat']),'t','vMag','theta','afterS','afterD',...
    'dtau','slipP','slipC','vSS','ySS','yr','Trec','f0','ruptureExtend','tracErr','runtime',...
    'completed','warnMsg','reversed','Tplus','Tminus','iC','inRupture')

%% figure
xc = rcv.xc;
isPinned = true(rcv.N,1);
isPinned(iC) = false;
verts = rcv.x/1e3.*[1,1,0];
mapE = @(us,ud) us.*rcv.sv(iC,1) + ud.*rcv.dv(iC,1);
mapN = @(us,ud) us.*rcv.sv(iC,2) + ud.*rcv.dv(iC,2);
drawPatches = @(val) patch('Vertices',verts,'Faces',rcv.vertices(iC,:),...
    'FaceVertexCData',val,'FaceColor','flat','EdgeColor',[1,1,1]*0.75);
drawAsperity = @() patch('Vertices',verts,'Faces',rcv.vertices(isPinned,:),...
    'FaceColor',[1,1,1]*0.6,'EdgeColor',[1,1,1]*0.4);
dmap = interp1([0;0.5;1],[0.15,0.3,0.75;1,1,1;0.75,0.15,0.15],linspace(0,1,24)');
mapAxes = @() set(gca,'FontSize',12,'LineWidth',1.2,'DataAspectRatio',[1,1,1]);

% patches for time series
probes = [-130,25;0,25;110,25;0,165;-150,90;-120,90];
probeNames = {'up-dip west','up-dip centre','up-dip east','down-dip','west flank','west edge (in extended rupture)'};
iProbe = zeros(size(probes,1),1);
for k = 1:size(probes,1)
    [~,iProbe(k)] = min(hypot(xc(iC,1)/1e3 - probes(k,1),xc(iC,2)/1e3 - probes(k,2)));
end
pcols = lines(size(probes,1));
markProbes = @() scatter(xc(iC(iProbe),1)/1e3,xc(iC(iProbe),2)/1e3,60,pcols,'filled','MarkerEdgeColor','k');
ruptBox = @() plot(ruptX([1,2,2,1,1])/1e3,aspW([1,1,2,2,1])*cosd(dip)/1e3,'k--','LineWidth',1.5);

figure(1),clf
set(gcf,'Color','w','Position',[50,50,1700,900])

subplot(2,3,1)
drawPatches(dtauV*plateDir'), hold on
drawAsperity()
quiver(xc(iC,1)/1e3,xc(iC,2)/1e3,4*mapE(dtauV(:,1),dtauV(:,2)),4*mapN(dtauV(:,1),dtauV(:,2)),0,'k')
markProbes(), ruptBox()
clim([-5,5]), colormap(gca,dmap)
cb = colorbar; cb.Label.String = '\Delta\tau along plate rake (MPa)';
title(sprintf('coseismic \\Delta\\tau, slip %.1f m, rupture dashed',max(hypot(slipP(:,1),slipP(:,2)))))
xlabel('East (km)'), ylabel('North (km)'), axis tight, box on, mapAxes()

for p = 1:2
    subplot(2,3,1+p)
    tShow = [1,Trec/yr];
    [~,it] = min(abs(t/yr - tShow(p)));
    aS = afterS(:,it); aD = afterD(:,it);
    drawPatches(hypot(aS,aD)), hold on
    drawAsperity()
    quiver(xc(iC,1)/1e3,xc(iC,2)/1e3,4*mapE(aS,aD),4*mapN(aS,aD),0,'k')
    markProbes(), ruptBox()
    clim([0,max(hypot(afterS(:,end),afterD(:,end)))]), colormap(gca,parula(20))
    cb = colorbar; cb.Label.String = 'afterslip (m)';
    title(sprintf('afterslip relative to steady creep at %.3g yr',t(it)/yr))
    xlabel('East (km)'), ylabel('North (km)'), axis tight, box on, mapAxes()
end

subplot(2,3,4)
for k = 1:numel(iProbe)
    loglog(t(2:end)/yr,vMag(iProbe(k),2:end)/Vpl,'-','Color',pcols(k,:),'LineWidth',2), hold on
end
yline(1,'k--')
grid on, box on
xlabel('time after event (yr)'), ylabel('|v| / V_{pl}')
legend(probeNames,'Location','southwest')
set(gca,'FontSize',12,'LineWidth',1.2)

subplot(2,3,5)
for k = 1:numel(iProbe)
    j = iProbe(k);
    semilogx(t(2:end)/yr,rad2deg(unwrap(theta(j,2:end))),'-','Color',pcols(k,:),'LineWidth',2), hold on
    yline(rad2deg(ySS(3*Nc+j)),':','Color',pcols(k,:),'LineWidth',1.5)
end
grid on, box on
xlabel('time after event (yr)'), ylabel('rake (deg), dotted: steady state')
set(gca,'FontSize',12,'LineWidth',1.2)

subplot(2,3,6)
dRakePlus = rad2deg(angle(exp(1i*(theta(:,1) - ySS(3*Nc+1:end)))));
drawPatches(dRakePlus), hold on
drawAsperity()
markProbes(), ruptBox()
clim(ceil(max(abs(dRakePlus)))*[-1,1]), colormap(gca,dmap)
cb = colorbar; cb.Label.String = 'rake^+ - rake_{ss} (deg)';
title('coseismic rotation of the slip-rate direction')
xlabel('East (km)'), ylabel('North (km)'), axis tight, box on, mapAxes()

sgtitle(sprintf('imposed event, T_{rec} = %g yr, rupture %g km x %g km, f_0 = %g, a-b = 0.005, plate rake %g°',...
    Trec/yr,diff(ruptX)/1e3,diff(aspW)/1e3,f0,rcv.Vrake))
print(fullfile(outDir,['postseismicExample',tag,'.png']),'-dpng','-r150')
