function [u1,u2,u3]=computeDisplacementTetrahedronShearZoneMixedQuad( ...
    x1,x2,x3,A,B,C,D,e11,e12,e13,e22,e23,e33,nu,varargin)
% function COMPUTEDISPLACEMENTTETRAHEDRONSHEARZONEMIXEDQUAD computes the
% displacement field associated with deforming tetrahedral strain volume
% considering using either the double-exponential or the Gauss-Legendre
% quadrature using the following geometry
%
%
%                      / North (x1)
%                     /
%        surface     /
%      -------------+-------------- East (x2)
%                  /|
%                 / |     + A
%                /  |    /  .
%                   |   /     .
%                   |  /        .
%                   | /           .
%                   |/              + B
%                   /            .  |
%                  /|          /    |
%                 / :       .       |
%                /  |    /          |
%               /   : .             |
%              /   /|               |
%             / .   :               |
%            +------|---------------+
%          C        :                 D
%                   |
%                   Down (x3)
%
%
% Input:
% x1, x2, x3         north, east coordinates and depth of the observation point,
% A, B, C, D         north, east, and depth coordinates of the vertices,
% eij                source strain component 11, 12, 13, 22, 23 and 33
%                    in the strain volume,
% nu                 Poisson's ratio in the half space.
%
% Output:
% u1                 displacement component in the north direction,
% u2                 displacement component in the east direction,
% u3                 displacement component in the down direction.
%
% Author: Sylvain Barbot (sbarbot@ntu.edu.sg) - June 2, 2018, Los Angeles.

assert(min(x3(:))>=0,'depth must be positive.');

% process optional input
p = inputParser;
p.addParameter('precision',0.01,@validatePrecision);
p.addParameter('bound',3.0,@validateBound);
p.addParameter('n',7,@validateN);
p.addParameter('radius',1.1,@validateRadius);
p.parse(varargin{:});
optionStruct = p.Results;

% circumsphere of tetrahedron
a=det([[A',1];[B',1];[C',1];[D',1]]);
O(1)=+det([[sum(A.^2),A([2,3])',1];[sum(B.^2),B([2,3])',1];[sum(C.^2),C([2,3])',1];[sum(D.^2),D([2,3])',1]])/2/a;
O(2)=-det([[sum(A.^2),A([1,3])',1];[sum(B.^2),B([1,3])',1];[sum(C.^2),C([1,3])',1];[sum(D.^2),D([1,3])',1]])/2/a;
O(3)=+det([[sum(A.^2),A([1,2])',1];[sum(B.^2),B([1,2])',1];[sum(C.^2),C([1,2])',1];[sum(D.^2),D([1,2])',1]])/2/a;

% radius of circumsphere
r=norm(O(:)-A(:));

inside=sqrt((x1-O(1)).^2+(x2-O(2)).^2+(x3-O(3)).^2)<optionStruct.radius*r;
outside=~inside;

% initiate empty array
u1=zeros(size(x1));
u2=zeros(size(x1));
u3=zeros(size(x1));

% numerical solution with Gauss-Legendre quadrature for points
% outside the circumsphere
if numel(x1(outside))>0
    [u1(outside),u2(outside),u3(outside)]= ...
        computeDisplacementTetrahedronShearZoneGauss( ...
            x1(outside),x2(outside),x3(outside), ...
            A,B,C,D,e11,e12,e13,e22,e23,e33,nu, ...
            'n',optionStruct.n);
end

% numerical solution with double-exponential quadrature for points
% inside the circumsphere
if numel(x1(inside))>0
    [u1(inside),u2(inside),u3(inside)]= ...
        computeDisplacementTetrahedronShearZoneTanhSinh( ...
            x1(inside),x2(inside),x3(inside), ...
            A,B,C,D,e11,e12,e13,e22,e23,e33,nu, ...
            'precision',optionStruct.precision, ...
            'bound',optionStruct.bound);
end

end

function p = validatePrecision(x)
if ~(isreal(x) && isscalar(x) && x > 0)
    error('MATLAB:invalidPrecision','invalid precision');
end
p = true;
end

function p = validateBound(x)
if ~(isreal(x) && isscalar(x) && x > 0)
    error('MATLAB:invalidBound','invalid truncation');
end
p = true;
end

function p = validateN(x)
if ~(x > 0)
    error('MATLAB:invalid','invalid number of integration points');
end
p = true;
end

function p = validateRadius(x)
if ~(isreal(x) && isscalar(x) && x > 0)
    error('MATLAB:invalidRadius','invalid radius');
end
p = true;
end

