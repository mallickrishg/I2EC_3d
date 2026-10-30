function [L11,L12,L13,L22,L23,L33]=computeDisplacementKernelsPolyhedronShearZone(shz,nu,x,vecsize)
% COMPUTEDISPLACEMENTKERNELSPOLYHEDRONSHEARZONE computes inversion kernels
% (Green's functions) to connect distributed inelastic strain to surface
% displacements (Barbot, 2018).
%
% INTERFACE:
%
%   [L11,L12,L13,L22,L23,L33] = ...
%            unicycle.greens.computeDisplacementKernelsPolyhedronShearZone(shz,nu,x,vecsize)
%
% INPUT:
%
% shz     geometry.polyhedron object
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
% Author: Sylvain Barbot, 07/03/2018, Singapore.
%
% SEE ALSO: unicycle

import unicycle.utils.*

% number of GPS stations
D=size(x,1);

% Green's function matrix
L=cell(6,1);

% individual strain components
e=eye(6);

textprogressbar('# polyhedral shear zone displacement kernels: ');

for j=1:6

    G=zeros(vecsize*D,shz.N);
    
    for i=1:shz.N
        
        if 0==mod(i-1,2)
            textprogressbar((i/shz.N+(j-1))*100/6);
        end

        [u1,u2,u3]=unicycle.greens.computeDisplacementPolyhedronShearZoneMixedQuad( ...
            x(:,2),x(:,1),-x(:,3), ...
            [shz.x(:,2),shz.x(:,1),-shz.x(:,3)],shz.vertices(shz.faces{i},:),...
            e(j,1),e(j,2),e(j,3),e(j,4),e(j,5),e(j,6),nu);
            
        switch(vecsize)
            case 2
                u=[u2,u1]';
            case 3
                u=[u2,u1,-u3]';
            otherwise
                error('greens:computeDisplacementKernelsPolyhedronShearZone','invalid degrees of freedom (%d)\n',...
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



