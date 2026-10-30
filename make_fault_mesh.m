% Generate a simple planar subduction megathrust as a triangular mesh
% written in unicycle .ned/.tri format
%
% .ned : vertex_index  north(y)  east(x)  depth(positive down)   [m]
% .tri : triangle_index  i1  i2  i3  rake
%
% geometry: trench runs along x (east), centred at x = 0, at y = 0 and
% the fault dips toward +y (north) at a uniform dip angle
%
% AUTHORS:
% Rishav Mallick, 2026

clear

% fault dimensions [m] and dip [degrees]
Lstrike = 500e3;
Wdip = 300e3;
dip = 15;
% depth of the up-dip edge [m]
ztop = 0;
% approximate triangle edge length [m] and output file prefix
meshes = {10e3,'fault/generic_sz15';...
          20e3,'fault/generic_sz15_dx20'};

for m = 1:size(meshes,1)
    writeMesh(meshes{m,2},meshes{m,1},Lstrike,Wdip,dip,ztop)
end
fname = meshes{end,2};

%% check that the mesh loads in unicycle and plot it
eM = unicycle.greens.nikkhoo15(30e3,0.25);
rcv = unicycle.geometry.triangleReceiver(fname,eM);

figure(1),clf
rcv.plotPatch(rcv.xc(:,3)./1e3), hold on
rcv.plotPatch()
axis tight equal, grid on, box on
cb = colorbar; cb.Label.String = 'depth (km)';
colormap(turbo(16))
xlabel('East (m)'), ylabel('North (m)'), zlabel('Up (m)')
view(-30,30)
set(gca,'FontSize',15,'LineWidth',1.5)

function writeMesh(fname,dx,Lstrike,Wdip,dip,ztop)
%% build structured grid of vertices in fault plane coordinates
nx = round(Lstrike/dx) + 1;
nw = round(Wdip/dx) + 1;
s = linspace(-Lstrike/2,Lstrike/2,nx);% along-strike
w = linspace(0,Wdip,nw);% along-dip (down-dip distance)
[S,W] = meshgrid(s,w);

% map to east, north, depth
east = S(:);
north = W(:).*cosd(dip);
depth = ztop + W(:).*sind(dip);

%% split each grid cell into 2 triangles, alternating the diagonal
% vertex index of grid node (i along-dip, j along-strike)
id = @(i,j) (j-1).*nw + i;
tri = zeros(2*(nx-1)*(nw-1),3);
k = 0;
for j = 1:nx-1
    for i = 1:nw-1
        a = id(i,j); b = id(i,j+1); c = id(i+1,j+1); d = id(i+1,j);
        if mod(i+j,2) == 0
            tri(k+1,:) = [a,b,c];
            tri(k+2,:) = [a,c,d];
        else
            tri(k+1,:) = [a,b,d];
            tri(k+2,:) = [b,c,d];
        end
        k = k + 2;
    end
end

%% write .ned and .tri files
nvert = numel(east);
ntri = size(tri,1);

fid = fopen([fname '.ned'],'w');
fprintf(fid,'%d\t%.3f\t%.3f\t%.3f\n',[(1:nvert);north';east';depth']);
fclose(fid);

fid = fopen([fname '.tri'],'w');
fprintf(fid,'%d\t%d\t%d\t%d\t%d\n',[(1:ntri);tri';zeros(1,ntri)]);
fclose(fid);

fprintf('wrote %s.ned (%d vertices) and %s.tri (%d triangles)\n',...
    fname,nvert,fname,ntri);
end
