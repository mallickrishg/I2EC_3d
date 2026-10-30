function [L11,L12,L13,L22,L23,L33]=computeDisplacementKernelsTetrahedronShearZone(shz,nu,x,vecsize)
% COMPUTEDISPLACEMENTKERNELSTETRAHEDRONSHEARZONE computes inversion kernels
% (Green's functions) to connect distributed inelastic strain to surface
% displacements (Barbot, 2018).
%
% INTERFACE:
%
%   [L11,L12,L13,L22,L23,L33] = ...
%            unicycle.greens.computeDisplacementKernelsTetrahedronShearZone(shz,nu,x,vecsize)
%
% INPUT:
%
% shz     geometry.tetrahedron object
% nu      Poisson's ratio
% x       coordinates of observation point (east, north, up)
% vecsize number of component for displacement vector
%
%
% DATA LAYOUT:
%
%            strain eij
%       /e1               \
%       |n1               |
%       |u1               |
%       |e2               |
%       |n2               |
%       |u2               |
% Lij = | .               |
%       | .               |
%       | .               |
%       | .               |
%       |en               |
%       |nn               |
%       \un               /
%
% Author: Sylvain Barbot, 06/26/2018, Los Angeles.
%
% SEE ALSO: unicycle

import unicycle.utils.*

% number of GPS stations
D=size(x,1);

% Green's function matrix
L=cell(6,1);

% individual strain components
e=eye(6);

textprogressbar('# tetrahedral shear zone displacement kernels: ');

for j=1:6

    G=zeros(vecsize*D,shz.N);
    
    for i=1:shz.N
        
        if 0==mod(i-1,2)
            textprogressbar((i/shz.N+(j-1))*100/6);
        end
        [u1,u2,u3]=unicycle.greens.computeDisplacementTetrahedronShearZoneMixedQuad( ...
            x(:,2),x(:,1),-x(:,3), ...
            shz.x(shz.vertices(k,1),:),shz.x(shz.vertices(k,2),:),shz.x(shz.vertices(k,3),:),shz.x(shz.vertices(k,4),:), ...
            e(j,1),e(j,2),e(j,3),e(j,4),e(j,5),e(j,6),nu);
            
        switch(vecsize)
            case 2
                u=[u2,u1]';
            case 3
                u=[u2,u1,-u3]';
            otherwise
                error('greens:computeDisplacementKernelsTetrahedronShearZone','invalid degrees of freedom (%d)\n',...
                    unicycle.greens.model.dgf);
        end
       
        G(:,i)=u(:);
        
    end
    
    L{j,1}=G;
    
    
end

% distribute kernels into two separate variables.
[L11,L12,L13,L22,L23,L33]=deal(L{:});

textprogressbar(100);
textprogressbar('');

end



