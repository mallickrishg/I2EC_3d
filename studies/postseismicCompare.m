% Compare the postseismicExample runs for f0 = 0.1 and 0.6, with the rupture
% equal to the asperity, extended 20 km along strike on each side with
% uniform slip, and extended with slip solved for a 3 MPa stress drop on
% every creeping patch inside the extension
%
% run studies/postseismicExample.m for the six cases first

clear
projectDir = fileparts(fileparts(mfilename('fullpath')));
cd(projectDir)
outDir = fullfile('studies','output');

f0s = [0.6,0.1];
% rows: {extension (m), file suffix, label}
cases = {0,'','rupture = asperity';
         20e3,'','+20 km each side, uniform slip';
         20e3,'_drop3','+20 km each side, 3 MPa stress drop on extension'};
nr = size(cases,1);
meshName = 'generic_sz15_dx20';
rcv = unicycle.geometry.triangleReceiver(fullfile('fault',meshName),unicycle.greens.nikkhoo15(30e3,0.25));
xc = rcv.xc;
Vpl = 0.054/3.1536e7;
plateRake = 68;

R = cell(nr,2);
for i = 1:nr
    for j = 1:2
        R{i,j} = load(fullfile(outDir,sprintf('postseismicExample_f0_%g_ext%gkm%s.mat',f0s(j),cases{i,1}/1e3,cases{i,2})));
    end
end
iC = R{1,1}.iC;
xcC = xc(iC,:)/1e3;
yr = R{1,1}.yr;
Nc = numel(iC);

% summary table
fprintf('%10s %6s %10s %10s %10s %9s %10s %10s %12s\n','case','f0','|T-| MPa','min|T+|','reversed',...
    'backward','max dθ+','runtime s','max T err')
for i = 1:nr
    for j = 1:2
        r = R{i,j};
        dthPlus = rad2deg(angle(exp(1i*(r.theta(:,1) - r.ySS(3*Nc+1:end)))));
        vPar = r.vMag.*cos(r.theta - deg2rad(plateRake));
        fprintf('%10d %6.1f %10.1f %10.3g %10d %9d %9.1f° %10.1f %12.3g\n',i,f0s(j),...
            median(hypot(r.Tminus(:,1),r.Tminus(:,2))),min(hypot(r.Tplus(:,1),r.Tplus(:,2))),...
            sum(r.reversed),sum(any(vPar < 0,2)),max(abs(dthPlus)),r.runtime,max(r.tracErr))
    end
end

%% figure: time series at probes and profiles
probes = [0,25;-150,90;-120,90];
probeNames = {'up-dip centre','west flank','west edge'};
iProbe = zeros(size(probes,1),1);
for k = 1:size(probes,1)
    [~,iProbe(k)] = min(hypot(xcC(:,1) - probes(k,1),xcC(:,2) - probes(k,2)));
end
pcols = lines(numel(iProbe));
styles = {'-','--'};
% along-strike profiles through the middle of the asperity and up-dip of it
profN = [90,25];
profName = {'north ≈ 90 km (mid-asperity)','north ≈ 25 km (up-dip)'};

figure(1),clf
set(gcf,'Color','w','Position',[50,50,2200,450*nr])
for i = 1:nr
    % velocity
    subplot(nr,5,5*(i-1)+1)
    for k = 1:numel(iProbe)
        for j = 1:2
            r = R{i,j};
            loglog(r.t(2:end)/yr,r.vMag(iProbe(k),2:end)/Vpl,styles{j},'Color',pcols(k,:),'LineWidth',2), hold on
        end
    end
    yline(1,'k:')
    grid on, box on, ylim([1e-6,1e9])
    xlabel('time after event (yr)'), ylabel('|v| / V_{pl}')
    title(cases{i,3})
    set(gca,'FontSize',12,'LineWidth',1.2)
    if i == 1
        h = gobjects(numel(iProbe)+2,1);
        for k = 1:numel(iProbe), h(k) = plot(nan,nan,'-','Color',pcols(k,:),'LineWidth',2); end
        h(end-1) = plot(nan,nan,'k-','LineWidth',2);
        h(end) = plot(nan,nan,'k--','LineWidth',2);
        legend(h,[probeNames,{'f_0 = 0.6','f_0 = 0.1'}],'Location','southwest','FontSize',10)
    end

    % rake relative to the plate rake
    subplot(nr,5,5*(i-1)+2)
    for k = 1:numel(iProbe)
        for j = 1:2
            r = R{i,j};
            semilogx(r.t(2:end)/yr,rad2deg(unwrap(r.theta(iProbe(k),2:end))),styles{j},'Color',pcols(k,:),'LineWidth',2), hold on
        end
    end
    yline(plateRake,'k:')
    grid on, box on
    xlabel('time after event (yr)'), ylabel('rake (deg)')
    set(gca,'FontSize',12,'LineWidth',1.2)

    % coseismic rake rotation along the mid-asperity profile and afterslip
    % along the up-dip profile
    for q = 1:2
        subplot(nr,5,5*(i-1)+2+q)
        sel = abs(xcC(:,2) - profN(q)) < 10;
        [xs,ord] = sort(xcC(sel,1));
        idx = find(sel);
        idx = idx(ord);
        for j = 1:2
            r = R{i,j};
            if q == 1
                val = rad2deg(angle(exp(1i*(r.theta(idx,1) - r.ySS(3*Nc+idx)))));
            else
                val = hypot(r.afterS(idx,end),r.afterD(idx,end));
            end
            plot(xs,val,[styles{j},'o'],'Color',[0.2,0.2,0.2]*j,'LineWidth',1.5,'MarkerSize',4), hold on
        end
        xline([-110,90],'k:'), xline([-110,90] + cases{i,1}/1e3*[-1,1],'r--')
        grid on, box on
        xlabel('East (km)')
        if q == 1
            ylabel('coseismic rake rotation (deg)')
        else
            ylabel('afterslip at 150 yr (m)')
        end
        title(profName{q})
        legend({'f_0 = 0.6','f_0 = 0.1'},'Location','best','FontSize',10)
        set(gca,'FontSize',12,'LineWidth',1.2)
    end
end
%% traction recovery at the probes
for i = 1:nr
    subplot(nr,5,5*i)
    for k = 1:numel(iProbe)
        for j = 1:2
            r = R{i,j};
            semilogx(r.t(2:end)/yr,r.Tmag(iProbe(k),2:end),styles{j},'Color',pcols(k,:),'LineWidth',2), hold on
        end
    end
    grid on, box on
    xlabel('time after event (yr)'), ylabel('|\tau| (MPa)')
    title('traction magnitude')
    set(gca,'FontSize',12,'LineWidth',1.2,'YScale','log')
end
sgtitle('f_0 = 0.6 (solid) vs 0.1 (dashed); dotted: asperity edges, red dashed: rupture edges')
print(fullfile(outDir,'postseismicCompare'),'-dpng','-r150')
