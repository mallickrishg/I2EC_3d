% Convergence of the interseismic ODE to the BEM steady state for
% different f0 at fixed a-b
%
% hypothesis: slip-rate magnitude relaxes on a time scale (dF/dζ)/(k Vpl),
% rake on |τ|/(k Vpl). With the log law, |τ|/(dF/dζ) = f0/(a-b) + log(v/v0),
% so rake convergence slows down as f0 grows while magnitude convergence
% does not depend on f0
%
% generic_sz15_dx20 mesh (20 km triangles), rectangular asperity pinned,
% oblique plate motion, asinh (or log) friction law, a-b = 0.005, v0 = Vpl so that
% f0 is the friction coefficient at plate rate, 1e4 yr from v = Vpl
% everywhere

% plate rake (deg) and friction law ('asinh' or 'log') can be set before
% running, e.g. plateRake = 90; law = 'log'; run(...)
clearvars -except plateRake law
if ~exist('plateRake','var'), plateRake = 68; end
if ~exist('law','var'), law = 'asinh'; end
yr = 3.1536e7;
projectDir = fileparts(fileparts(mfilename('fullpath')));
cd(projectDir)
outDir = fullfile('studies','output');
if ~exist(outDir,'dir'), mkdir(outDir), end

meshName = 'generic_sz15_dx20';
f0list = [0.1,0.2,0.4,0.6];
tmax = 1e4*yr;
tout = [0,logspace(-1,4,101)]*yr;

eM = unicycle.greens.nikkhoo15(30e3,0.25);
rcv = unicycle.geometry.triangleReceiver(fullfile('fault',meshName),eM);
bigK = unicycle.greens.computeKernels(rcv,meshName);

% user-provided properties (m, s, MPa)
rcv.a = 0.015;
rcv.b = 0.010;
rcv.sigma = 50;
rcv.Vs = 3e3;
rcv.Vpl = 0.054/yr;
rcv.Vo = rcv.Vpl;
rcv.Vrake = plateRake;
% locked asperity: 200 km along strike x 100 km down dip, edges on the 20 km
% grid lines (the grid has nodes at x = ±10, ±30, ... km, so the rectangle
% is offset 10 km west of x = 0)
addpath(fullfile(projectDir,'studies'))
dip = 15;
aspX = [-110e3,90e3];
aspW = [40e3,140e3];
rcv.pinnedPosition = rectangularAsperity(rcv,aspX,aspW,dip);

%% BEM steady state (independent of friction)
rcv.mu0 = f0list(1);
sys = unicycle.ode.assembleODE(rcv,bigK,law);
Nc = sys.Nc;
Vpl = rcv.Vpl;
vbem = sys.Kcc \ (-sys.taudotLoad);
vbMag = hypot(vbem(1:Nc),vbem(Nc+1:end));
vbRake = atan2(vbem(Nc+1:end),vbem(1:Nc));
% rake is only meaningful where the patch is actually sliding
wellDefined = vbMag > 0.2*Vpl;
fprintf('%s: creeping %d, rake evaluated on %d patches with |v_bem| > 0.2 Vpl\n',meshName,Nc,sum(wellDefined))
fprintf('BEM: |v|/Vpl in [%.3f, %.3f], rake in [%.1f, %.1f] deg\n',...
    min(vbMag)/Vpl,max(vbMag)/Vpl,min(rad2deg(vbRake)),max(rad2deg(vbRake)))

%% integrate for each f0
opts = odeset('RelTol',1e-8,'AbsTol',1e-8);
nf = numel(f0list);
magErr = zeros(numel(tout),nf);     % max | |v| - |v_bem| | / Vpl
magErrRms = zeros(numel(tout),nf);
rakeErr = zeros(numel(tout),nf);    % max |Δrake| (deg) where well defined
rakeErrRms = zeros(numel(tout),nf);
vecErr = zeros(numel(tout),nf);     % max |v - v_bem| / Vpl
finalMag = zeros(Nc,nf);
finalRake = zeros(Nc,nf);
runtime = zeros(1,nf);
ratio = zeros(1,nf);                % |τ|/(dF/dζ) at v = Vpl

for j = 1:nf
    rcv.mu0 = f0list(j);
    sys = unicycle.ode.assembleODE(rcv,bigK,law);
    [tauPl,dFPl] = unicycle.ode.frictionLaw(sys,log(Vpl./sys.v0));
    ratio(j) = tauPl(1)/dFPl(1);

    fprintf('f0 = %.2f: integrating ...\n',f0list(j))
    tic
    [t,y] = ode45(@(t,y) unicycle.ode.odeRateDependent(t,y,sys),tout,sys.y0,opts);
    runtime(j) = toc;

    zeta = y(:,2*Nc+1:3*Nc)';
    theta = y(:,3*Nc+1:end)';
    vMag = sys.v0.*exp(zeta);
    dRake = rad2deg(angle(exp(1i*(theta - vbRake))));
    vs = vMag.*cos(theta);
    vd = vMag.*sin(theta);

    magErr(:,j) = max(abs(vMag - vbMag),[],1)'/Vpl;
    magErrRms(:,j) = sqrt(mean((vMag - vbMag).^2,1))'/Vpl;
    rakeErr(:,j) = max(abs(dRake(wellDefined,:)),[],1)';
    rakeErrRms(:,j) = sqrt(mean(dRake(wellDefined,:).^2,1))';
    vecErr(:,j) = max(hypot(vs - vbem(1:Nc),vd - vbem(Nc+1:end)),[],1)'/Vpl;
    finalMag(:,j) = vMag(:,end);
    finalRake(:,j) = theta(:,end);

    fprintf('f0 = %.2f: %5.1f s, %d steps\n',f0list(j),runtime(j),numel(t))
end

%% time to converge
tol = [1e-1,1e-2];
tc = @(err,tl) arrayfun(@(j) firstBelow(t/yr,err(:,j),tl),1:nf);

fprintf('\n|τ|/(dF/dζ) at Vpl is the predicted ratio of rake to magnitude time scales\n')
fprintf('%6s %10s %10s %10s %10s %10s %10s %12s %12s\n','f0','|τ|/dF','dF/(a-b)σ',...
    't|v|<0.1','t|v|<0.01','tθ<1°','tθ<0.1°','|v| err 1e4','θ err 1e4')
for j = 1:nf
    rcv.mu0 = f0list(j);
    sys = unicycle.ode.assembleODE(rcv,bigK,law);
    [~,dFPl] = unicycle.ode.frictionLaw(sys,log(Vpl./sys.v0));
    fprintf('%6.2f %10.1f %10.4f %10.1f %10.1f %10.1f %10.1f %12.2e %11.3f°\n',f0list(j),ratio(j),...
        dFPl(1)/sys.abSigma(1),...
        firstBelow(t/yr,magErr(:,j),0.1),firstBelow(t/yr,magErr(:,j),0.01),...
        firstBelow(t/yr,rakeErr(:,j),1),firstBelow(t/yr,rakeErr(:,j),0.1),...
        magErr(end,j),rakeErr(end,j))
end

%% compare final states between f0 (all should equal the BEM solution)
fprintf('\nat 1e4 yr, difference from the f0 = %.2f run:\n',f0list(1))
for j = 2:nf
    fprintf('  f0 = %.2f: max ||v| diff| %.2e Vpl, max rake diff %.3f deg (well-defined patches)\n',...
        f0list(j),max(abs(finalMag(:,j)-finalMag(:,1)))/Vpl,...
        max(abs(rad2deg(angle(exp(1i*(finalRake(wellDefined,j)-finalRake(wellDefined,1))))))))
end

tag = sprintf('_rake%g',plateRake);
if ~strcmp(law,'asinh'), tag = [tag,'_',law]; end
save(fullfile(outDir,['f0ConvergenceStudy',tag,'.mat']),'f0list','t','magErr','magErrRms',...
    'rakeErr','rakeErrRms','vecErr','finalMag','finalRake','vbMag','vbRake',...
    'wellDefined','ratio','runtime','yr','plateRake','law')

%% figure
cols = parula(nf+1);
figure(1),clf
set(gcf,'Color','w','Position',[100,100,1300,800])
panels = {magErr,'max ||v| - |v_{BEM}|| / V_{pl}';...
          rakeErr,'max |rake - rake_{BEM}| (deg)';...
          magErrRms,'rms ||v| - |v_{BEM}|| / V_{pl}';...
          rakeErrRms,'rms |rake - rake_{BEM}| (deg)'};
for p = 1:4
    subplot(2,2,p)
    for j = 1:nf
        loglog(t(2:end)/yr,panels{p,1}(2:end,j),'-','Color',cols(j,:),'LineWidth',2), hold on
    end
    grid on, box on
    xlabel('time (yr)'), ylabel(panels{p,2})
    set(gca,'FontSize',13,'LineWidth',1.2)
    if p == 1
        legend(arrayfun(@(f) sprintf('f_0 = %.2f',f),f0list,'UniformOutput',false),'Location','southwest')
    end
end
sgtitle(['Convergence to the BEM steady state, ',strrep(meshName,'_','\_'),', a-b = 0.005, v_0 = V_{pl}, ',law,' law, plate rake ',num2str(plateRake),'°'])
print(fullfile(outDir,['f0ConvergenceStudy',tag]),'-dpng','-r150')

%% final slip-rate vectors on the mesh
% map-view projection of v = vs sv + vd dv, plotted at patch centroids.
% First panel: BEM steady state coloured by rake. Other panels: final state
% of each run coloured by rake - rake_BEM (blank where |v_BEM| < 0.2 Vpl and
% the rake is ill-defined), black arrows = run, grey arrows = BEM.
% Last panel: zoom on the asperity, slip-rate directions of BEM and all runs.
iC = sys.iC;
xc = rcv.xc;
isPinned = true(rcv.N,1);
isPinned(iC) = false;
mapE = @(vs,vd) vs.*rcv.sv(iC,1) + vd.*rcv.dv(iC,1);
mapN = @(vs,vd) vs.*rcv.sv(iC,2) + vd.*rcv.dv(iC,2);
arrowScale = 15e3/Vpl;   % 15 km of arrow per Vpl
arrows = @(vs,vd,col,lw) quiver(xc(iC,1)/1e3,xc(iC,2)/1e3,...
    arrowScale*mapE(vs,vd)/1e3,arrowScale*mapN(vs,vd)/1e3,0,...
    'Color',col,'LineWidth',lw,'MaxHeadSize',0.4);
verts = rcv.x/1e3.*[1,1,0];
drawPatches = @(val) patch('Vertices',verts,'Faces',rcv.vertices(iC,:),...
    'FaceVertexCData',val,'FaceColor','flat','EdgeColor',[1,1,1]*0.75);
drawAsperity = @(col) patch('Vertices',verts,'Faces',rcv.vertices(isPinned,:),...
    'FaceColor',col,'EdgeColor',[1,1,1]*0.4);
aspBox = @() plot(aspX([1,2,2,1,1])/1e3,aspW([1,1,2,2,1])*cosd(dip)/1e3,'k-','LineWidth',1.5);
% blue-white-red
dmap = interp1([0;0.5;1],[0.15,0.3,0.75;1,1,1;0.75,0.15,0.15],linspace(0,1,24)');
dRakeLim = 30;
vsF = finalMag.*cos(finalRake);
vdF = finalMag.*sin(finalRake);

figure(2),clf
set(gcf,'Color','w','Position',[100,100,1600,850])
subplot(2,3,1)
drawPatches(rad2deg(vbRake)), hold on
drawAsperity([1,1,1]*0.6)
arrows(vbem(1:Nc),vbem(Nc+1:end),'k',0.8)
aspBox()
axis equal tight, box on
clim(plateRake + [-50,60]), colormap(gca,turbo(22))
cb = colorbar; cb.Label.String = 'rake (deg)';
title('BEM steady state')
xlabel('East (km)'), ylabel('North (km)')
set(gca,'FontSize',12,'LineWidth',1.2)
for j = 1:nf
    subplot(2,3,j+1)
    dr = rad2deg(angle(exp(1i*(finalRake(:,j) - vbRake))));
    dr(~wellDefined) = NaN;
    drawPatches(dr), hold on
    drawAsperity([1,1,1]*0.6)
    arrows(vbem(1:Nc),vbem(Nc+1:end),[1,1,1]*0.55,0.8)
    arrows(vsF(:,j),vdF(:,j),'k',0.8)
    aspBox()
    axis equal tight, box on
    clim(dRakeLim*[-1,1]), colormap(gca,dmap)
    cb = colorbar; cb.Label.String = 'rake - rake_{BEM} (deg)';
    title(sprintf('f_0 = %.1f at 10^4 yr: max |\\Deltarake| = %.1f°',f0list(j),max(abs(dr),[],'omitnan')))
    xlabel('East (km)'), ylabel('North (km)')
    set(gca,'FontSize',12,'LineWidth',1.2)
end
subplot(2,3,6)
% directions only (unit vectors, no heads), so rake differences show as fans
dirLen = 14;   % km
unitArrows = @(rake,col,lw) quiver(xc(iC,1)/1e3,xc(iC,2)/1e3,...
    dirLen*mapE(cos(rake),sin(rake)),dirLen*mapN(cos(rake),sin(rake)),0,...
    'Color',col,'LineWidth',lw,'ShowArrowHead','off');
drawAsperity([1,1,1]*0.85), hold on
hq = unitArrows(vbRake,'k',2.5);
for j = 1:nf
    hq(j+1) = unitArrows(finalRake(:,j),cols(j,:),1.2);
end
aspBox()
axis equal, box on
xlim([-200,170]), ylim([0,170])
cb = colorbar; cb.Visible = 'off';   % keep the panel the same width
legend(hq,[{'BEM'},arrayfun(@(f) sprintf('f_0 = %.1f',f),f0list,'UniformOutput',false)],...
    'Location','northoutside','NumColumns',3)
title('slip-rate direction near the asperity')
xlabel('East (km)'), ylabel('North (km)')
set(gca,'FontSize',12,'LineWidth',1.2)
sgtitle(sprintf('slip-rate vectors at 10^4 yr (%g km per V_{pl}), plate rake %g°, %s law, locked asperity grey',arrowScale*Vpl/1e3,plateRake,law))
print(fullfile(outDir,['f0ConvergenceStudyVectors',tag]),'-dpng','-r150')

function tc = firstBelow(t,err,tol)
% first time after which err stays below tol (Inf if never)
above = find(err >= tol,1,'last');
if isempty(above)
    tc = t(1);
elseif above == numel(err)
    tc = Inf;
else
    tc = t(above+1);
end
end
