% Leading Karhunen-Loève modes on the creeping patches of the 20 km mesh,
% for correlation lengths L/4, L/2 (default) and L
%
% exponential (Matérn nu = 1/2) covariance, area weighted, centred on the
% mean (the mean level of log eta is a separate parameter)

clear
projectDir = fileparts(fileparts(mfilename('fullpath')));
cd(projectDir)
addpath(fullfile(projectDir,'studies'))
outDir = fullfile('studies','output');
if ~exist(outDir,'dir'), mkdir(outDir), end

meshName = 'generic_sz15_dx20';
rcv = unicycle.geometry.triangleReceiver(fullfile('fault',meshName),unicycle.greens.nikkhoo15(30e3,0.25));
dip = 15;
aspX = [-110e3,90e3];
aspW = [40e3,140e3];
pinned = rectangularAsperity(rcv,aspX,aspW,dip);
isPinned = false(rcv.N,1);
isPinned(pinned) = true;
iC = find(~isPinned);

b = unicycle.inverse.klBasis(rcv,iC);
L = b.L;
fracs = [0.25,0.5,1];
nModes = 4;
B = cell(numel(fracs),1);
for k = 1:numel(fracs)
    B{k} = unicycle.inverse.klBasis(rcv,iC,struct('ell',fracs(k)*L,'rMax',30));
    fprintf('ell = %4.0f km (%.2f L): variance captured by 1, 2, 3, 5, 10 modes: %s\n',...
        fracs(k)*L/1e3,fracs(k),mat2str(B{k}.varFraction([1,2,3,5,10])',3))
end

verts = rcv.x/1e3.*[1,1,0];
dmap = interp1([0;0.5;1],[0.15,0.3,0.75;1,1,1;0.75,0.15,0.15],linspace(0,1,24)');
figure(1),clf
set(gcf,'Color','w','Position',[50,50,1800,800])
for k = 1:numel(fracs)
    for j = 1:nModes
        subplot(numel(fracs),nModes + 1,(k-1)*(nModes + 1) + j)
        val = B{k}.phi(:,j);
        patch('Vertices',verts,'Faces',rcv.vertices(iC,:),'FaceVertexCData',val,...
            'FaceColor','flat','EdgeColor','none'), hold on
        patch('Vertices',verts,'Faces',rcv.vertices(isPinned,:),'FaceColor',[1,1,1]*0.6,'EdgeColor','none')
        clim(max(abs(val))*[-1,1]), colormap(gca,dmap)
        axis equal tight, box on
        set(gca,'FontSize',10,'XTick',[],'YTick',[])
        title(sprintf('mode %d (%.0f%%)',j,100*B{k}.lambda(j)/sum(B{k}.lambda)*B{k}.varFraction(end)))
        if j == 1
            ylabel(['corr. length ',sprintf('%.0f km = %.2g L',fracs(k)*L/1e3,fracs(k))],'FontWeight','bold')
        end
    end
    subplot(numel(fracs),nModes + 1,k*(nModes + 1))
    plot(1:numel(B{k}.varFraction),B{k}.varFraction,'o-','LineWidth',1.5)
    grid on, box on, ylim([0,1])
    xlabel('number of modes'), ylabel('variance captured')
    set(gca,'FontSize',10)
end
sgtitle(sprintf('KL modes of an exponential covariance on the creeping patches (%s, L = %.0f km), mean removed',...
    strrep(meshName,'_','\_'),L/1e3))
print(fullfile(outDir,'klModes.png'),'-dpng','-r150')
