% Surface velocity and displacement time series over the periodic
% earthquake cycle, at a trench-normal profile of stations
%
% generic_sz15_dx20 mesh, 200 km x 100 km asperity, rake 68 deg,
% Trec = 150 yr, asinh law, a-b = 0.005, v0 = Vpl, f0 = 0.1 and 0.6.
% The cycle is spun up with unicycle.ode.runCycles; surface deformation is
% computed relative to the long-term motion (removeLongTerm) and includes
% the coseismic step at t = 0.

clear
yr = 3.1536e7;
projectDir = fileparts(fileparts(mfilename('fullpath')));
cd(projectDir)
addpath(fullfile(projectDir,'studies'))
outDir = fullfile('studies','output');
if ~exist(outDir,'dir'), mkdir(outDir), end

meshName = 'generic_sz15_dx20';
Trec = 150*yr;
f0list = [0.6,0.1];

eM = unicycle.greens.nikkhoo15(30e3,0.25);
rcv = unicycle.geometry.triangleReceiver(fullfile('fault',meshName),eM);
bigK = unicycle.greens.computeKernels(rcv,meshName);
rcv.a = 0.015;
rcv.b = 0.010;
rcv.sigma = 50;
rcv.Vs = 3e3;
rcv.Vpl = 0.054/yr;
rcv.Vo = rcv.Vpl;
rcv.Vrake = 68;
dip = 15;
aspX = [-110e3,90e3];
aspW = [40e3,140e3];
rcv.pinnedPosition = rectangularAsperity(rcv,aspX,aspW,dip);

% stations: trench-normal profile through the middle of the asperity
% (no station on the fault trace at north = 0, where the kernels are singular)
north = (-25:10:285)'*1e3;
xObs = [-10e3*ones(size(north)),north,zeros(size(north))];
G = unicycle.greens.computeSurfaceKernels(rcv,xObs);
probeN = [40,90,160,250]*1e3;
[~,iSta] = min(abs(north - probeN),[],1);

R = cell(numel(f0list),1);
for j = 1:numel(f0list)
    rcv.mu0 = f0list(j);
    sys = unicycle.ode.assembleODE(rcv,bigK,'asinh');
    event = unicycle.ode.periodicEvent(sys,Trec);
    out = unicycle.ode.runCycles(sys,event,struct('verbose',0));
    o = struct('removeLongTerm',true,'pinnedSlip',out.event.slipP);
    [vRel,uRel] = unicycle.ode.surfaceDeformation(G,sys,out.t,out.y,o);
    [vAbs,~] = unicycle.ode.surfaceDeformation(G,sys,out.t,out.y);
    % over a periodic cycle the slip budget closes, so the displacement
    % relative to the long-term motion returns to zero (coseismic +
    % postseismic + interseismic = G Vpl Trec)
    fprintf('f0 = %.1f: %d cycles; |u_rel(Trec)| max %.2e m vs coseismic step max %.2f m\n',...
        f0list(j),out.nCycles,max(abs(uRel(end,:,:)),[],'all'),...
        max(abs(uRel(1,:,:)),[],'all'))
    R{j} = struct('t',out.t,'vRel',vRel,'uRel',uRel,'vAbs',vAbs,'f0',f0list(j));
end
save(fullfile(outDir,'surfaceTimeseries.mat'),'R','xObs','yr','Trec')

%% figure
comps = {'east','north','up'};
styles = {'-','--'};
pcols = lines(numel(iSta));
tShow = [1,10,50,149]*yr;
tcols = parula(numel(tShow)+1);

figure(1),clf
set(gcf,'Color','w','Position',[50,50,1800,950])
for c = 1:3
    % displacement time series at probe stations
    subplot(3,3,c)
    for k = 1:numel(iSta)
        for j = 1:numel(R)
            r = R{j};
            plot(r.t/yr,r.uRel(:,iSta(k),c),styles{j},'Color',pcols(k,:),'LineWidth',2), hold on
        end
    end
    grid on, box on
    xlabel('time after event (yr)'), ylabel(sprintf('%s displacement (m)',comps{c}))
    title(sprintf('%s, relative to long-term motion',comps{c}))
    set(gca,'FontSize',12,'LineWidth',1.2)
    if c == 1
        h = gobjects(numel(iSta)+2,1);
        for k = 1:numel(iSta), h(k) = plot(nan,nan,'-','Color',pcols(k,:),'LineWidth',2); end
        h(end-1) = plot(nan,nan,'k-','LineWidth',2);
        h(end) = plot(nan,nan,'k--','LineWidth',2);
        legend(h,[arrayfun(@(n) sprintf('north = %g km',n),north(iSta)/1e3,'UniformOutput',false)',...
            {'f_0 = 0.6','f_0 = 0.1'}],'Location','best','FontSize',9)
    end

    % velocity time series (log time) at probe stations
    subplot(3,3,3+c)
    for k = 1:numel(iSta)
        for j = 1:numel(R)
            r = R{j};
            semilogx(r.t(2:end)/yr,r.vRel(2:end,iSta(k),c)*yr*1e3,styles{j},'Color',pcols(k,:),'LineWidth',2), hold on
        end
    end
    grid on, box on
    xlabel('time after event (yr)'), ylabel(sprintf('%s velocity (mm/yr)',comps{c}))
    set(gca,'FontSize',12,'LineWidth',1.2,'YScale','linear')
    late = R{1}.t > 0.1*yr;
    xlim([1e-3,Trec/yr])
    ylim([min(R{1}.vRel(late,iSta,c),[],'all'),max(R{1}.vRel(late,iSta,c),[],'all')]*yr*1e3 + [-5,5])

    % velocity profiles at selected times
    subplot(3,3,6+c)
    for q = 1:numel(tShow)
        for j = 1:numel(R)
            r = R{j};
            [~,it] = min(abs(r.t - tShow(q)));
            plot(north/1e3,r.vRel(it,:,c)*yr*1e3,styles{j},'Color',tcols(q,:),'LineWidth',2), hold on
        end
    end
    xline(aspW*cosd(dip)/1e3,'k:')
    grid on, box on
    xlabel('north (km), trench at 0'), ylabel(sprintf('%s velocity (mm/yr)',comps{c}))
    if c == 1
        legend(arrayfun(@(s) sprintf('%g yr',s),reshape(repmat(tShow/yr,numel(R),1),1,[]),'UniformOutput',false),...
            'Location','best','FontSize',9)
    end
    set(gca,'FontSize',12,'LineWidth',1.2)
end
sgtitle(sprintf('surface deformation over the periodic cycle (relative to long-term motion), profile at east = -10 km, T_{rec} = %g yr, rake %g°',...
    Trec/yr,rcv.Vrake))
print(fullfile(outDir,'surfaceTimeseries.png'),'-dpng','-r150')
