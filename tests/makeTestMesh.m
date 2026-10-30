function basename = makeTestMesh(folder,name,Lstrike,Wdip,dip,dx)
% MAKETESTMESH writes a small planar fault in unicycle .ned/.tri format
% for tests. Same geometry convention as make_fault_mesh.m: trench along
% east at north = 0, dipping toward north, top at the surface.
%
%   basename = makeTestMesh(folder,name,Lstrike,Wdip,dip,dx)
%
% defaults: 60 km x 30 km, dip 20 degrees, dx = 10 km (36 triangles)

if nargin < 3, Lstrike = 60e3; end
if nargin < 4, Wdip = 30e3; end
if nargin < 5, dip = 20; end
if nargin < 6, dx = 10e3; end

nx = round(Lstrike/dx) + 1;
nw = round(Wdip/dx) + 1;
[S,W] = meshgrid(linspace(-Lstrike/2,Lstrike/2,nx),linspace(0,Wdip,nw));
east = S(:);
north = W(:).*cosd(dip);
depth = W(:).*sind(dip);

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

basename = fullfile(folder,name);
fid = fopen([basename '.ned'],'w');
fprintf(fid,'%d\t%.3f\t%.3f\t%.3f\n',[(1:numel(east));north';east';depth']);
fclose(fid);
fid = fopen([basename '.tri'],'w');
fprintf(fid,'%d\t%d\t%d\t%d\t%d\n',[(1:size(tri,1));tri';zeros(1,size(tri,1))]);
fclose(fid);

end
