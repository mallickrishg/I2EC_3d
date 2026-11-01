% Figures for studies/stochasticSequence.m: spin-up and variability of the
% late-interseismic coupling under a stochastic event sequence, against
% the deterministic sequence of mean events
%
% Coupling C = 1 - |v|/Vpl is taken at the end of every 150 yr cycle, just
% before the next west event (120 yr after the east event).
%
% parameters that can be set before running (defaults in brackets):
%   f0     - [0.3]
%   seeds  - stochastic runs to load [1 2 3]
%   nBurn  - cycles discarded as spin-up in the statistics [50]

clearvars -except f0 seeds nBurn
if ~exist('f0','var'), f0 = 0.3; end
if ~exist('seeds','var'), seeds = [1,2,3]; end
if ~exist('nBurn','var'), nBurn = 50; end
projectDir = fileparts(fileparts(mfilename('fullpath')));
cd(projectDir)
outDir = fullfile('studies','output');

ref = load(fullfile(outDir,sprintf('stochasticSequence_f0_%g_mean.mat',f0)));
nS = numel(seeds);
for i = 1:nS
    R(i) = load(fullfile(outDir,sprintf('stochasticSequence_f0_%g_seed%d.mat',f0,seeds(i)))); %#ok<SAGROW>
end
yr = ref.yr; Tcyc = ref.Tcyc; Vpl = ref.Vpl; dip = ref.dip; seg = ref.seg;
nCycles = R(1).nCycles;
Nc = ref.Nc;
rcv = unicycle.geometry.triangleReceiver(fullfile('fault',ref.meshName),unicycle.greens.nikkhoo15(30e3,0.25));
xc = rcv.xc;
w = xc(:,2)/cosd(dip);
win = nBurn+1:nCycles;
nWin = numel(win);

% deterministic reference: periodicity of the last cycles
dRef = max(abs(ref.coupling(:,end) - ref.coupling(:,end-1)));
fprintf('reference: max |dC| between the last two cycles %.2e\n',dRef)
Cref = ref.coupling(:,end);

% pooled statistics over the stationary window
Cpool = cell2mat(arrayfun(@(r) r.coupling(:,win),R,'UniformOutput',false));
Cmean = mean(Cpool,2);
Cstd = std(Cpool,0,2);
chi = mean(cell2mat(arrayfun(@(r) sum(r.seisCycle(:,win),2),R,'UniformOutput',false)),2)/(Vpl*Tcyc*nWin);
chiRef = ref.seisCycle(:,end)/(Vpl*Tcyc);

% nominal segment boxes and probes
inSeg = false(Nc,numel(seg.x0));
for j = 1:numel(seg.x0)
    inSeg(:,j) = abs(xc(:,1) - seg.x0(j)) < ref.Lmean/2 & abs(w - seg.w0(j)) < ref.Lmean/2;
end
probes = [seg.x0(:),seg.w0(:);-75e3,90e3;0,25e3;0,170e3];
probeNames = [strcat(seg.name,' centre'),{'W-C gap','up-dip','down-dip'}];
iProbe = zeros(size(probes,1),1);
for k = 1:size(probes,1)
    [~,iProbe(k)] = min(hypot(xc(:,1) - probes(k,1),w - probes(k,2)));
end
pcols = lines(numel(iProbe));

fprintf('stationary window: cycles %d-%d, %d seeds\n',win(1),nCycles,nS)
for j = 1:numel(seg.x0)
    Cs = mean(Cpool(inSeg(:,j),:),1);
    fprintf('%6s segment: C ref %.3f, stochastic mean %.3f, std %.3f, 5-95%% [%.3f, %.3f]; seismic fraction %.2f (ref %.2f)\n',...
        seg.name{j},mean(Cref(inSeg(:,j))),mean(Cs),std(Cs),prctile(Cs,5),prctile(Cs,95),...
        mean(chi(inSeg(:,j))),mean(chiRef(inSeg(:,j))))
end

% memory: correlation of C(n) with the coseismic slip a patch received in
% cycle n - lag (lag -1 = the next cycle, not yet happened), and with the
% slip deficit Vpl t - s at the end of cycle n
lags = -1:6;
active = chi > 0.2;
rho = nan(nnz(active),numel(lags));
rhoDef = nan(nnz(active),1);
iA = find(active);
for p = 1:numel(iA)
    for l = 1:numel(lags)
        c = []; s = [];
        for i = 1:nS
            n = win(win - lags(l) >= 1 & win - lags(l) <= nCycles);
            c = [c,R(i).coupling(iA(p),n)]; %#ok<AGROW>
            s = [s,R(i).seisCycle(iA(p),n - lags(l))]; %#ok<AGROW>
        end
        if std(s) > 0, rho(p,l) = corr(c(:),s(:)); end
    end
    c = []; d = [];
    for i = 1:nS
        c = [c,R(i).coupling(iA(p),win)]; %#ok<AGROW>
        d = [d,Vpl*Tcyc*win - R(i).slipEndPar(iA(p),win)]; %#ok<AGROW>
    end
    rhoDef(p) = corr(c(:),d(:));
end
fprintf('correlation of C with the slip deficit at the end of the cycle: median %.2f over %d patches\n',...
    median(rhoDef,'omitnan'),numel(iA))

save(fullfile(outDir,sprintf('stochasticSequenceCompare_f0_%g.mat',f0)),'Cref','Cmean','Cstd','chi','chiRef',...
    'lags','rho','rhoDef','iA','nBurn','seeds','inSeg','iProbe')

%% figure 1: maps
verts = rcv.x/1e3.*[1,1,0];
drawPatches = @(val) patch('Vertices',verts,'Faces',rcv.vertices,...
    'FaceVertexCData',val,'FaceColor','flat','EdgeColor',[1,1,1]*0.75);
box2 = @(x,wc,L,varargin) plot((x + L/2*[-1,1,1,-1,-1])/1e3,(wc + L/2*[-1,-1,1,1,-1])*cosd(dip)/1e3,varargin{:});
segBoxes = @() arrayfun(@(j) box2(seg.x0(j),seg.w0(j),ref.Lmean,'k--','LineWidth',1.2),1:numel(seg.x0));
mapAxes = @() set(gca,'FontSize',12,'LineWidth',1.2,'DataAspectRatio',[1,1,1]);
dmap = interp1([0;0.5;1],[0.15,0.3,0.75;1,1,1;0.75,0.15,0.15],linspace(0,1,24)');
cmapC = flipud(gray(20));
markProbes = @() scatter(xc(iProbe,1)/1e3,xc(iProbe,2)/1e3,60,pcols,'filled','MarkerEdgeColor','k');

figure(1),clf
set(gcf,'Color','w','Position',[50,50,1800,900])
subplot(2,3,1)
drawPatches(Cref), hold on, segBoxes(), markProbes()
clim([0,1]), colormap(gca,cmapC)
cb = colorbar; cb.Label.String = 'C = 1 - |v|/V_{pl}';
title('mean-event sequence (periodic)')
xlabel('East (km)'), ylabel('North (km)'), axis tight, box on, mapAxes()

subplot(2,3,2)
drawPatches(Cmean), hold on, segBoxes()
clim([0,1]), colormap(gca,cmapC)
cb = colorbar; cb.Label.String = 'mean C';
title(sprintf('stochastic: mean over cycles %d-%d, %d seeds',win(1),nCycles,nS))
xlabel('East (km)'), ylabel('North (km)'), axis tight, box on, mapAxes()

subplot(2,3,3)
drawPatches(Cstd), hold on, segBoxes()
colormap(gca,parula(20))
cb = colorbar; cb.Label.String = 'std of C';
title('stochastic: cycle-to-cycle std of C')
xlabel('East (km)'), ylabel('North (km)'), axis tight, box on, mapAxes()

subplot(2,3,4)
drawPatches(Cmean - Cref), hold on, segBoxes()
m = max(abs(Cmean - Cref)); clim(m*[-1,1]), colormap(gca,dmap)
cb = colorbar; cb.Label.String = 'mean C - C_{ref}';
title('stochastic mean - mean-event sequence')
xlabel('East (km)'), ylabel('North (km)'), axis tight, box on, mapAxes()

subplot(2,3,5)
drawPatches(chi), hold on, segBoxes()
clim([0,1]), colormap(gca,flipud(hot(20)))
cb = colorbar; cb.Label.String = 'seismic slip / V_{pl} t';
title('stochastic: seismic fraction of the slip budget')
xlabel('East (km)'), ylabel('North (km)'), axis tight, box on, mapAxes()

subplot(2,3,6)
r1 = R(1);
drawPatches(r1.coupling(:,end)), hold on
evLast = find(r1.catalog.cycle == nCycles);
for k = evLast'
    box2(r1.catalog.xc(k),r1.catalog.wc(k),r1.catalog.L(k),'-','Color',[0.85,0.1,0.1],'LineWidth',1.5);
    text(r1.catalog.xc(k)/1e3,(r1.catalog.wc(k) + r1.catalog.L(k)/2)*cosd(dip)/1e3 + 8,...
        sprintf('%.1f m',r1.catalog.slip(k)),'HorizontalAlignment','center','Color',[0.85,0.1,0.1],'FontSize',11)
end
clim([0,1]), colormap(gca,cmapC)
cb = colorbar; cb.Label.String = 'C';
title(sprintf('seed %d, end of cycle %d (red: its events)',seeds(1),nCycles))
xlabel('East (km)'), ylabel('North (km)'), axis tight, box on, mapAxes()
sgtitle(sprintf(['stochastic sequence, no locked zone: 3 square events per %g yr cycle at 0/20/30 yr, ',...
    'slip N(%.1f,%g) m, edge N(%g,%g) km, f_0 = %g; coupling at the end of each cycle'],...
    Tcyc/yr,ref.slipMean,ref.slipStd,ref.Lmean/1e3,ref.Lstd/1e3,f0))
print(fullfile(outDir,sprintf('stochasticSequenceMaps_f0_%g.png',f0)),'-dpng','-r150')

%% figure 2: spin-up, variability, memory
figure(2),clf
set(gcf,'Color','w','Position',[50,50,1800,900])
scols = lines(numel(seg.x0));

subplot(2,3,1)
for j = 1:numel(seg.x0)
    for i = 1:nS
        plot(1:nCycles,mean(R(i).coupling(inSeg(:,j),:),1),'-','Color',[scols(j,:),0.5],'LineWidth',1), hold on
    end
    plot(1:nCycles,mean(ref.coupling(inSeg(:,j),:),1),'--','Color',scols(j,:)*0.6,'LineWidth',2)
end
xline(nBurn,'k:')
grid on, box on
xlabel('cycle'), ylabel('C averaged over the nominal segment')
title('spin-up (thin: seeds, dashed: mean events)')
set(gca,'FontSize',12,'LineWidth',1.2)

subplot(2,3,2)
for k = 1:numel(iProbe)
    d = (Vpl*Tcyc*(1:nCycles) - R(1).slipEndPar(iProbe(k),:));
    plot((1:nCycles)*Tcyc/yr/1e3,d,'-','Color',pcols(k,:),'LineWidth',1.5), hold on
end
grid on, box on
xlabel('time (kyr)'), ylabel('slip deficit V_{pl}t - s (m)')
legend(probeNames,'Location','best')
title(sprintf('seed %d: slip deficit at the end of each cycle',seeds(1)))
set(gca,'FontSize',12,'LineWidth',1.2)

subplot(2,3,3)
edges = linspace(0,1,41);
for j = 1:numel(seg.x0)
    Cs = mean(Cpool(inSeg(:,j),:),1);
    hh(j) = histogram(Cs,edges,'Normalization','pdf','FaceColor',scols(j,:),'FaceAlpha',0.4,'EdgeColor','none'); hold on
    xline(mean(Cref(inSeg(:,j))),'--','Color',scols(j,:)*0.6,'LineWidth',2)
end
grid on, box on
xlabel('C averaged over the nominal segment'), ylabel('pdf')
legend(hh,seg.name,'Location','northwest')
title('stationary distribution (dashed: mean events)')
set(gca,'FontSize',12,'LineWidth',1.2)

subplot(2,3,4)
plot(lags,median(rho,1,'omitnan'),'ko-','LineWidth',2,'MarkerFaceColor','k'), hold on
plot(lags,prctile(rho,[25,75],1),'k:','LineWidth',1)
yline(0,'k-')
grid on, box on
xlabel('lag (cycles; -1 = next cycle)'), ylabel('corr(C_n, seismic slip in cycle n - lag)')
title(sprintf('memory: median and IQR over %d ruptured patches',numel(iA)))
set(gca,'FontSize',12,'LineWidth',1.2)

subplot(2,3,5)
d1 = R(1); dn = d1.dense;
tEv = d1.catalog.t(d1.catalog.cycle == nCycles);
for k = 1:numel(iProbe)
    tt = []; vv = [];
    for e = 1:numel(dn)
        tt = [tt;dn(e).t(2:end)/yr - tEv(1)]; %#ok<AGROW>
        vv = [vv;exp(dn(e).y(2:end,2*Nc+iProbe(k)))]; %#ok<AGROW>
    end
    semilogy(tt,vv,'-','Color',pcols(k,:),'LineWidth',1.5), hold on
end
yline(1,'k:')
xlim([0,Tcyc/yr]), grid on, box on
xlabel('time in cycle (yr)'), ylabel('|v| / V_{pl}')
title(sprintf('seed %d, cycle %d',seeds(1),nCycles))
set(gca,'FontSize',12,'LineWidth',1.2)

subplot(2,3,6)
c = []; d = [];
j = 2; p = iProbe(j);
for i = 1:nS
    c = [c,R(i).coupling(p,win)]; %#ok<AGROW>
    d = [d,Vpl*Tcyc*win - R(i).slipEndPar(p,win)]; %#ok<AGROW>
end
scatter(d,c,12,'filled','MarkerFaceAlpha',0.4), hold on
grid on, box on
xlabel('slip deficit V_{pl}t - s (m)'), ylabel('C at the end of the cycle')
title(sprintf('%s: C vs slip deficit (median corr over patches %.2f)',probeNames{j},median(rhoDef,'omitnan')))
set(gca,'FontSize',12,'LineWidth',1.2)
print(fullfile(outDir,sprintf('stochasticSequenceStats_f0_%g.png',f0)),'-dpng','-r150')
