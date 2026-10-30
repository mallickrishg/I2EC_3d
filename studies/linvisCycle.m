% Analytic periodic cycle of the linear viscous model, compared with the
% rate-and-state periodic cycle
%
% generic_sz15_dx20 mesh, 200 km x 100 km asperity, rake 68 deg,
% Trec = 150 yr. Linear viscous friction τ = τ0 + η v with
%   (1) uniform η = (a-b)σ/Vpl, the tangent viscosity of RSF at plate rate
%   (2) η increasing tenfold with depth
% against the asinh RSF cycle with a-b = 0.005, σ = 50 MPa, f0 = 0.1
% (closest to an isotropic viscosity) spun up with unicycle.ode.runCycles.

clear
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
rcv.a = 0.015;
rcv.b = 0.010;
rcv.sigma = 50;
rcv.mu0 = 0.1;
rcv.Vs = 3e3;
rcv.Vpl = 0.054/yr;
rcv.Vo = rcv.Vpl;
rcv.Vrake = 68;
dip = 15;
aspX = [-110e3,90e3];
aspW = [40e3,140e3];
rcv.pinnedPosition = rectangularAsperity(rcv,aspX,aspW,dip);
xc = rcv.xc;

%% rate-and-state reference
sys = unicycle.ode.assembleODE(rcv,bigK,'asinh');
Nc = sys.Nc;
iC = sys.iC;
tic
out = unicycle.ode.runCycles(sys,unicycle.ode.periodicEvent(sys,Trec),struct('verbose',0));
fprintf('RSF (f0 = %g): %d cycles to periodic, %.1f s\n',rcv.mu0,out.nCycles,toc)
tR = out.t;
vR = sys.v0.*exp(out.y(:,2*Nc+1:3*Nc)').*cat(3,cos(out.y(:,3*Nc+1:end)'),sin(out.y(:,3*Nc+1:end)'));
vR = [vR(:,:,1);vR(:,:,2)];

%% linear viscous, analytic
eta0 = 0.25/rcv.Vpl;
z = -xc(iC,3);
etaCases = {eta0*ones(Nc,1),eta0*10.^(z/max(z) - 0.5)};
names = {'linear, uniform \eta','linear, \eta \times 10 with depth'};
L = cell(2,1);
for k = 1:2
    tic
    lin = unicycle.ode.linvisAssemble(rcv,bigK,etaCases{k});
    tA = toc;
    tic
    event = unicycle.ode.periodicEvent(lin,Trec);
    [~,info] = unicycle.ode.linvisSpinUp(lin,event);
    [v,s] = unicycle.ode.linvisEvolve(lin,info.vPost,tR);
    tS = toc;
    VplC = lin.Vpl(iC,:)*Trec;
    budget = hypot(info.slipCycle(:,1),info.slipCycle(:,2))./hypot(VplC(:,1),VplC(:,2));
    fprintf('%s: assemble %.2f s, spin-up + %d output times %.3f s, cond(V) %.0f, spectral radius %.3f, slip/(Vpl T) in [%.6f, %.6f], min v.vplhat %.2g Vpl\n',...
        names{k},tA,numel(tR),tS,lin.condV,info.spectralRadius,min(budget),max(budget),...
        min(sum(reshape(v(:,end),[],2).*lin.Vpl(iC,:),2)./sum(lin.Vpl(iC,:).^2,2)))
    L{k} = struct('lin',lin,'info',info,'v',v,'s',s,'event',event);
end

%% surface stations
north = [25;60;150;250]*1e3;
xObs = [-10e3*ones(size(north)),north,zeros(size(north))];
G = unicycle.greens.computeSurfaceKernels(rcv,xObs);
o = struct('removeLongTerm',true,'pinnedSlip',out.event.slipP);
[~,uR] = unicycle.ode.surfaceDeformation(G,sys,tR,out.y,o);
uL = cell(2,1);
for k = 1:2
    [~,uL{k}] = unicycle.ode.surfaceDeformation(G,L{k}.lin,tR,struct('v',L{k}.v,'s',L{k}.s),o);
end

%% figure
probes = [0,25;-150,90;0,165];
probeNames = {'up-dip centre','west flank','down-dip'};
iProbe = zeros(size(probes,1),1);
for k = 1:size(probes,1)
    [~,iProbe(k)] = min(hypot(xc(iC,1)/1e3 - probes(k,1),xc(iC,2)/1e3 - probes(k,2)));
end
pcols = lines(4);
styles = {'-','--',':'};
isPinned = true(rcv.N,1);
isPinned(iC) = false;
verts = rcv.x/1e3.*[1,1,0];
spd = @(v,i) hypot(v(i,:),v(Nc+i,:));

figure(1),clf
set(gcf,'Color','w','Position',[50,50,1800,900])

subplot(2,3,1)
for k = 1:numel(iProbe)
    loglog(tR(2:end)/yr,spd(vR(:,2:end),iProbe(k))/rcv.Vpl,styles{1},'Color',pcols(k,:),'LineWidth',2), hold on
    for c = 1:2
        loglog(tR(2:end)/yr,spd(L{c}.v(:,2:end),iProbe(k))/rcv.Vpl,styles{c+1},'Color',pcols(k,:),'LineWidth',2)
    end
end
yline(1,'k:')
grid on, box on
xlabel('time after event (yr)'), ylabel('|v| / V_{pl}')
title('slip rate: blue up-dip centre, red west flank, yellow down-dip')
h = [plot(nan,nan,'k-','LineWidth',2);plot(nan,nan,'k--','LineWidth',2);plot(nan,nan,'k:','LineWidth',2)];
legend(h,{sprintf('RSF, f_0 = %g',rcv.mu0),names{:}},'Location','southwest')
set(gca,'FontSize',12,'LineWidth',1.2)

subplot(2,3,2)
patch('Vertices',verts,'Faces',rcv.vertices(iC,:),'FaceVertexCData',log10(etaCases{2}),...
    'FaceColor','flat','EdgeColor',[1,1,1]*0.75), hold on
patch('Vertices',verts,'Faces',rcv.vertices(isPinned,:),'FaceColor',[1,1,1]*0.6,'EdgeColor',[1,1,1]*0.4)
scatter(xc(iC(iProbe),1)/1e3,xc(iC(iProbe),2)/1e3,60,pcols(1:numel(iProbe),:),'filled','MarkerEdgeColor','k')
cb = colorbar; cb.Label.String = 'log_{10} \eta (MPa s/m)';
axis equal tight, box on
xlabel('East (km)'), ylabel('North (km)')
title('depth-dependent \eta')
set(gca,'FontSize',12,'LineWidth',1.2)

subplot(2,3,3)
for k = 1:2
    tau = -1./real(L{k}.lin.lambda)/yr;
    histogram(log10(tau),40,'DisplayStyle','stairs','LineWidth',2), hold on
end
xline(log10(Trec/yr),'k--')
grid on, box on
xlabel('log_{10} relaxation time -1/Re(\lambda) (yr)'), ylabel('number of modes')
legend([names,{'T_{rec}'}],'Location','northwest')
title('spectrum of A = D^{-1}K_{cc}')
set(gca,'FontSize',12,'LineWidth',1.2)

comps = {'east','north','up'};
for c = 1:3
    subplot(2,3,3+c)
    for k = 1:numel(north)
        plot(tR/yr,uR(:,k,c),styles{1},'Color',pcols(k,:),'LineWidth',2), hold on
        for q = 1:2
            plot(tR/yr,uL{q}(:,k,c),styles{q+1},'Color',pcols(k,:),'LineWidth',2)
        end
    end
    grid on, box on
    xlabel('time after event (yr)'), ylabel(sprintf('%s displacement (m)',comps{c}))
    title(sprintf('%s, relative to long-term motion',comps{c}))
    if c == 1
        hs = gobjects(numel(north),1);
        for k = 1:numel(north), hs(k) = plot(nan,nan,'-','Color',pcols(k,:),'LineWidth',2); end
        legend(hs,arrayfun(@(n) sprintf('north = %g km',n),north/1e3,'UniformOutput',false),'Location','best')
    end
    set(gca,'FontSize',12,'LineWidth',1.2)
end
sgtitle(sprintf('analytic linear viscous cycle vs RSF, T_{rec} = %g yr, rake %g°, \\eta_0 = (a-b)\\sigma/V_{pl} = %.2g MPa s/m',...
    Trec/yr,rcv.Vrake,eta0))
print(fullfile(outDir,'linvisCycle.png'),'-dpng','-r150')
