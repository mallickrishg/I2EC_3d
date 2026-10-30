function [s11,s12,s13,s22,s23,s33]=computeStressTetrahedronShearZoneMixedQuad( ...
    x1,x2,x3,A,B,C,D,e11,e12,e13,e22,e23,e33,G,nu,varargin)
% function COMPUTESTRESSTETRAHEDRONSHEARZONEMIXEDQUAD computes
% the stress field associated with deforming tetrahedron strain volume
% considering the following geometry using the double-exponential or
% Gauss-Legendre numerical quadrature based on the distance from the 
% circumsphere center.
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
% Input:
% x1, x2, x3         north, east coordinates and depth of the observation point,
% A, B, C, D         north, east and depth coordinates of the vertices,
% eij                source strain component in the tetrahedron strain volume,
% G, nu              shear modulus and Poisson's ratio in the half-space.
%
% Options:
%
% 'precision' [0.01] for double-exponential quadrature
% 'bound' [3.0]      for double-exponential quadrature
% 'n' [7]            for Gauss-Legendre quadrature
% 'radius' [1.1]     scaling factor for circumsphere radius
%
% Output:
% sij                stress components at coordinates x1, x2, x3.
%
% Author: Sylvain Barbot (sbarbot@ntu.edu.sg) - June 1, 2018, Los Angeles.

import unicycle.greens.*

assert(min(x3(:))>=0,'depth must be positive.');
assert(A(3)>=0,'depth must be positive.');
assert(B(3)>=0,'depth must be positive.');
assert(C(3)>=0,'depth must be positive.');
assert(D(3)>=0,'depth must be positive.');

% process optional input
p = inputParser;
p.addParameter('precision',0.01,@validatePrecision);
p.addParameter('bound',3.0,@validateBound);
p.addParameter('n',7,@validateN);
p.addParameter('radius',1.1,@validateRadius);
p.parse(varargin{:});
optionStruct = p.Results;

% vector coordinates
A=A(:);
B=B(:);
C=C(:);
D=D(:);

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
s11=zeros(size(x1));
s12=zeros(size(x1));
s13=zeros(size(x1));
s22=zeros(size(x1));
s23=zeros(size(x1));
s33=zeros(size(x1));

% numerical solution with Gauss-Legendre quadrature for points outside the circumsphere
if numel(x1(outside))>0
    [s11(outside),s12(outside),s13(outside),s22(outside),s23(outside),s33(outside)]= ...
        unicycle.greens.computeStressTetrahedronShearZoneGauss( ...
        x1(outside),x2(outside),x3(outside),A,B,C,D,e11,e12,e13,e22,e23,e33,G,nu, ...
        'n',optionStruct.n);
end

% numerical solution with double-exponential quadrature for points inside the circumsphere
if numel(x1(inside))>0
    [s11(inside),s12(inside),s13(inside),s22(inside),s23(inside),s33(inside)]= ...
        unicycle.greens.computeStressTetrahedronShearZoneTanhSinh( ...
        x1(inside),x2(inside),x3(inside),A,B,C,D,e11,e12,e13,e22,e23,e33,G,nu, ...
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
if ~(isreal(x) && isscalar(x))
    error('MATLAB:invalidRadius','invalid radius');
end
p = true;
end
