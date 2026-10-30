% plot the generic_sz15_dx20 mesh with the pinned (locked) asperity used in
% f0ConvergenceStudy highlighted

clear
projectDir = fileparts(fileparts(mfilename('fullpath')));
cd(projectDir)
meshName = 'generic_sz15_dx20';

eM = unicycle.greens.nikkhoo15(30e3,0.25);
rcv = unicycle.geometry.triangleReceiver(fullfile('fault',meshName),eM);
% locked asperity: 200 km along strike x 100 km down dip, edges on the 20 km
% grid lines (the grid has nodes at x = ±10, ±30, ... km, so the rectangle
% is offset 10 km west of x = 0)
addpath(fullfile(projectDir,'studies'))
dip = 15;
aspX = [-110e3,90e3];
aspW = [40e3,140e3];
rcv.pinnedPosition = rectangularAsperity(rcv,aspX,aspW,dip);
isPinned = false(rcv.N,1);
isPinned(rcv.pinnedPosition) = true;
fprintf('%s: %d patches, %d pinned (%.1f%% of area)\n',meshName,rcv.N,...
    sum(isPinned),100*sum(rcv.area(isPinned))/sum(rcv.area))

figure(1),clf
set(gcf,'Color','w','Position',[100,100,1200,800])
for p = 1:2
    subplot(2,1,p)
    rcv.plotPatch(double(isPinned)), hold on
    rcv.plotPatch()
    % asperity outline
    ox = aspX([1,2,2,1,1]); ow = aspW([1,1,2,2,1]);
    plot3(ox,ow*cosd(dip),-ow*sind(dip)+50,'r-','LineWidth',2)
    colormap([0.85,0.9,1;0.8,0.15,0.15]), clim([0,1])
    axis tight equal, grid on, box on
    xlabel('East (m)'), ylabel('North (m)'), zlabel('Up (m)')
    set(gca,'FontSize',13,'LineWidth',1.2)
    if p == 1
        view(-30,30)
        title(sprintf('%s: %d triangles, locked asperity (red) %d triangles',...
            strrep(meshName,'_','\_'),rcv.N,sum(isPinned)))
    else
        view(2)
        title(sprintf('map view, red line: x in [%g, %g] km, down-dip distance in [%g, %g] km',aspX/1e3,aspW/1e3))
    end
end
print(fullfile('studies','output','lockedAsperity'),'-dpng','-r150')
