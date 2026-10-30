function [u1,u2,u3]=computeDisplacementTetrahedronShearZoneTanhSinh( ...
    x1,x2,x3,A,B,C,D,e11,e12,e13,e22,e23,e33,nu,varargin)
% function COMPUTEDISPLACEMENTTETRAHEDRONSHEARZONETANHSINH computes the
% displacement field associated with deforming tetrahedral strain volume
% considering the following geometry using the tanh-sinh numerical
% quadrature.
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
% Options:
% 'precision',real   sampling [0.001]
% 'bound',real       truncation [3.5]
%
% Output:
% u1                 displacement component in the north direction,
% u2                 displacement component in the east direction,
% u3                 displacement component in the down direction.
%
% Author: Sylvain Barbot (sbarbot@ntu.edu.sg) - Feb 9, 2018, Los Angeles.

assert(min(x3(:))>=0,'depth must be positive.');

% process optional input
p = inputParser;
p.addParameter('precision',0.001,@validatePrecision);
p.addParameter('bound',3.5,@validateBound);
p.parse(varargin{:});
optionStruct = p.Results;

% Lame parameter
lambda=2*nu/(1-2*nu);

% isotropic strain
ekk=e11+e22+e33;

% unit normal vectors
nA=cross(C-B,D-B); % BCD
nB=cross(D-C,A-C); % CDA
nC=cross(A-D,B-D); % DAB
nD=cross(B-A,C-A); % ABC

nA=nA/norm(nA);
nB=nB/norm(nB);
nC=nC/norm(nC);
nD=nD/norm(nD);

% check that unit vectors are pointing outward
if (nA'*(A(:)-(B(:)+C(:)+D(:))/3))>0
    nA=-nA;
end
if (nB'*(B(:)-(C(:)+D(:)+A(:))/3))>0
    nB=-nB;
end
if (nC'*(C(:)-(D(:)+A(:)+B(:))/3))>0
    nC=-nC;
end
if (nD'*(D(:)-(A(:)+B(:)+C(:))/3))>0
    nD=-nD;
end

% area of triangle ABC
ABC=norm(cross(C-A,B-A))/2;
% area of triangle BCD
BCD=norm(cross(D-B,C-B))/2;
% area of triangle CDA
CDA=norm(cross(A-C,D-C))/2;
% area of triangle DAB
DAB=norm(cross(B-D,A-D))/2;

% Radii
r1=@(y1,y2,y3) sqrt((x1-y1).^2+(x2-y2).^2+(x3-y3).^2);
r2=@(y1,y2,y3) sqrt((x1-y1).^2+(x2-y2).^2+(x3+y3).^2);

% Green's functions
G11=@(y1,y2,y3) 1/(16*pi*(1-nu))*( ...
    (3-4*nu)./r1(y1,y2,y3)+1./r2(y1,y2,y3)+(x1-y1).^2./r1(y1,y2,y3).^3 ...
    +(3-4*nu)*(x1-y1).^2./r2(y1,y2,y3).^3+2*x3.*y3.*(r2(y1,y2,y3).^2-3*(x1-y1).^2)./r2(y1,y2,y3).^5 ...
    +4*(1-2*nu)*(1-nu)*(r2(y1,y2,y3).^2-(x1-y1).^2+r2(y1,y2,y3).*(x3+y3))./(r2(y1,y2,y3).*(r2(y1,y2,y3)+x3+y3).^2)...
    );
G12=@(y1,y2,y3) (x1-y1).*(x2-y2)/(16*pi*(1-nu)).*( ...
    1./r1(y1,y2,y3).^3+(3-4*nu)./r2(y1,y2,y3).^3-6*x3.*y3./r2(y1,y2,y3).^5 ...
    -4*(1-2*nu)*(1-nu)./(r2(y1,y2,y3).*(r2(y1,y2,y3)+x3+y3).^2) ...
    );
G13=@(y1,y2,y3) (x1-y1)/(16*pi*(1-nu)).*( ...
    (x3-y3)./r1(y1,y2,y3).^3+(3-4*nu)*(x3-y3)./r2(y1,y2,y3).^3 ...
    -6*x3.*y3.*(x3+y3)./r2(y1,y2,y3).^5+4*(1-2*nu)*(1-nu)./(r2(y1,y2,y3).*(r2(y1,y2,y3)+x3+y3)) ...
    );
G21=@(y1,y2,y3) (x1-y1).*(x2-y2)/(16*pi*(1-nu)).*( ...
    1./r1(y1,y2,y3).^3+(3-4*nu)./r2(y1,y2,y3).^3-6*x3.*y3./r2(y1,y2,y3).^5 ...
    -4*(1-2*nu)*(1-nu)./(r2(y1,y2,y3).*(r2(y1,y2,y3)+x3+y3).^2) ...
    );
G22=@(y1,y2,y3) 1/(16*pi*(1-nu))*( ...
    (3-4*nu)./r1(y1,y2,y3)+1./r2(y1,y2,y3)+(x2-y2).^2./r1(y1,y2,y3).^3 ...
    +(3-4*nu)*(x2-y2).^2./r2(y1,y2,y3).^3+2*x3.*y3.*(r2(y1,y2,y3).^2-3*(x2-y2).^2)./r2(y1,y2,y3).^5 ...
    +4*(1-2*nu)*(1-nu)*(r2(y1,y2,y3).^2-(x2-y2).^2+r2(y1,y2,y3).*(x3+y3))./(r2(y1,y2,y3).*(r2(y1,y2,y3)+x3+y3).^2)...
    );
G23=@(y1,y2,y3) (x2-y2)/(16*pi*(1-nu)).*( ...
    (x3-y3)./r1(y1,y2,y3).^3+(3-4*nu)*(x3-y3)./r2(y1,y2,y3).^3 ...
    -6*x3.*y3.*(x3+y3)./r2(y1,y2,y3).^5+4*(1-2*nu)*(1-nu)./(r2(y1,y2,y3).*(r2(y1,y2,y3)+x3+y3)) ...
    );
G31=@(y1,y2,y3) (x1-y1)/(16*pi*(1-nu)).*( ...
    (x3-y3)./r1(y1,y2,y3).^3+(3-4*nu)*(x3-y3)./r2(y1,y2,y3).^3 ...
    +6*x3.*y3.*(x3+y3)./r2(y1,y2,y3).^5-4*(1-2*nu)*(1-nu)./(r2(y1,y2,y3).*(r2(y1,y2,y3)+x3+y3)) ...
    );
G32=@(y1,y2,y3) (x2-y2)/(16*pi*(1-nu)).*( ...
    (x3-y3)./r1(y1,y2,y3).^3+(3-4*nu)*(x3-y3)./r2(y1,y2,y3).^3 ...
    +6*x3.*y3.*(x3+y3)./r2(y1,y2,y3).^5-4*(1-2*nu)*(1-nu)./(r2(y1,y2,y3).*(r2(y1,y2,y3)+x3+y3)) ...
    );
G33=@(y1,y2,y3) 1/(16*pi*(1-nu))*( ...
    (3-4*nu)./r1(y1,y2,y3)+(5-12*nu+8*nu^2)./r2(y1,y2,y3)+(x3-y3).^2./r1(y1,y2,y3).^3 ...
    +6*x3.*y3.*(x3+y3).^2./r2(y1,y2,y3).^5+((3-4*nu)*(x3+y3).^2-2*x3.*y3)./r2(y1,y2,y3).^3 ...
    );

% parameterized surface integral
y=@(u,v,A,B,C) A*(1-u)*(1-v)/4+B*(1+u)*(1-v)/4+C*(1+v)/2;

% moment density
m11=lambda*ekk+2*e11;
m12=2*e12;
m13=2*e13;
m22=lambda*ekk+2*e22;
m23=2*e23;
m33=lambda*ekk+2*e33;

% numerical solution with tanh/sinh quadrature
h=optionStruct.precision;
n=fix(1/h*optionStruct.bound);

u1=zeros(size(x1));
u2=zeros(size(x2));
u3=zeros(size(x3));

for k=-n:n
    fprintf('k=%d/%d\n',k,n);
    wk=(0.5*h*pi*cosh(k*h))./(cosh(0.5*pi*sinh(k*h))).^2;
    xk=tanh(0.5*pi*sinh(k*h));
    for j=-n:n
        wj=(0.5*h*pi*cosh(j*h))./(cosh(0.5*pi*sinh(j*h))).^2;
        xj=tanh(0.5*pi*sinh(j*h));
        u1=u1+wk*wj*(1-xk)*IU1(xj,xk);
        u2=u2+wk*wj*(1-xk)*IU2(xj,xk);
        u3=u3+wk*wj*(1-xk)*IU3(xj,xk);
    end
end

    function d = IU1(u,v)
        % function IU1 is the integrand for displacement component u1
        d=zeros(size(x1));
        if (m11*nD(1)+m12*nD(2)+m13*nD(3)) ~= 0
            d=d+ABC/4*(m11*nD(1)+m12*nD(2)+m13*nD(3))*G11(y(u,v,A(1),B(1),C(1)),y(u,v,A(2),B(2),C(2)),y(u,v,A(3),B(3),C(3)));
        end
        if (m12*nD(1)+m22*nD(2)+m23*nD(3)) ~= 0
            d=d+ABC/4*(m12*nD(1)+m22*nD(2)+m23*nD(3))*G21(y(u,v,A(1),B(1),C(1)),y(u,v,A(2),B(2),C(2)),y(u,v,A(3),B(3),C(3)));
        end
        if (m13*nD(1)+m23*nD(2)+m33*nD(3)) ~= 0
            d=d+ABC/4*(m13*nD(1)+m23*nD(2)+m33*nD(3))*G31(y(u,v,A(1),B(1),C(1)),y(u,v,A(2),B(2),C(2)),y(u,v,A(3),B(3),C(3)));
        end
        if (m11*nA(1)+m12*nA(2)+m13*nA(3)) ~= 0
            d=d+BCD/4*(m11*nA(1)+m12*nA(2)+m13*nA(3))*G11(y(u,v,B(1),C(1),D(1)),y(u,v,B(2),C(2),D(2)),y(u,v,B(3),C(3),D(3)));
        end
        if (m12*nA(1)+m22*nA(2)+m23*nA(3)) ~= 0
            d=d+BCD/4*(m12*nA(1)+m22*nA(2)+m23*nA(3))*G21(y(u,v,B(1),C(1),D(1)),y(u,v,B(2),C(2),D(2)),y(u,v,B(3),C(3),D(3)));
        end
        if (m13*nA(1)+m23*nA(2)+m33*nA(3)) ~= 0
            d=d+BCD/4*(m13*nA(1)+m23*nA(2)+m33*nA(3))*G31(y(u,v,B(1),C(1),D(1)),y(u,v,B(2),C(2),D(2)),y(u,v,B(3),C(3),D(3)));
        end
        if (m11*nB(1)+m12*nB(2)+m13*nB(3)) ~= 0
            d=d+CDA/4*(m11*nB(1)+m12*nB(2)+m13*nB(3))*G11(y(u,v,C(1),D(1),A(1)),y(u,v,C(2),D(2),A(2)),y(u,v,C(3),D(3),A(3)));
        end
        if (m12*nB(1)+m22*nB(2)+m23*nB(3)) ~= 0
            d=d+CDA/4*(m12*nB(1)+m22*nB(2)+m23*nB(3))*G21(y(u,v,C(1),D(1),A(1)),y(u,v,C(2),D(2),A(2)),y(u,v,C(3),D(3),A(3)));
        end
        if (m13*nB(1)+m23*nB(2)+m33*nB(3)) ~= 0
            d=d+CDA/4*(m13*nB(1)+m23*nB(2)+m33*nB(3))*G31(y(u,v,C(1),D(1),A(1)),y(u,v,C(2),D(2),A(2)),y(u,v,C(3),D(3),A(3)));
        end
        if (m11*nC(1)+m12*nC(2)+m13*nC(3)) ~= 0
            d=d+DAB/4*(m11*nC(1)+m12*nC(2)+m13*nC(3))*G11(y(u,v,D(1),A(1),B(1)),y(u,v,D(2),A(2),B(2)),y(u,v,D(3),A(3),B(3)));
        end
        if (m12*nC(1)+m22*nC(2)+m23*nC(3)) ~= 0
            d=d+DAB/4*(m12*nC(1)+m22*nC(2)+m23*nC(3))*G21(y(u,v,D(1),A(1),B(1)),y(u,v,D(2),A(2),B(2)),y(u,v,D(3),A(3),B(3)));
        end
        if (m13*nC(1)+m23*nC(2)+m33*nC(3)) ~= 0
            d=d+DAB/4*(m13*nC(1)+m23*nC(2)+m33*nC(3))*G31(y(u,v,D(1),A(1),B(1)),y(u,v,D(2),A(2),B(2)),y(u,v,D(3),A(3),B(3)));
        end
    end

    function d = IU2(u,v)
        % function IU2 is the integrand for displacement component u2
        d=zeros(size(x2));
        if (m11*nD(1)+m12*nD(2)+m13*nD(3)) ~= 0
            d=d+ABC/4*(m11*nD(1)+m12*nD(2)+m13*nD(3))*G12(y(u,v,A(1),B(1),C(1)),y(u,v,A(2),B(2),C(2)),y(u,v,A(3),B(3),C(3)));
        end
        if (m12*nD(1)+m22*nD(2)+m23*nD(3)) ~= 0
            d=d+ABC/4*(m12*nD(1)+m22*nD(2)+m23*nD(3))*G22(y(u,v,A(1),B(1),C(1)),y(u,v,A(2),B(2),C(2)),y(u,v,A(3),B(3),C(3)));
        end
        if (m13*nD(1)+m23*nD(2)+m33*nD(3)) ~= 0
            d=d+ABC/4*(m13*nD(1)+m23*nD(2)+m33*nD(3))*G32(y(u,v,A(1),B(1),C(1)),y(u,v,A(2),B(2),C(2)),y(u,v,A(3),B(3),C(3)));
        end
        if (m11*nA(1)+m12*nA(2)+m13*nA(3)) ~= 0
            d=d+BCD/4*(m11*nA(1)+m12*nA(2)+m13*nA(3))*G12(y(u,v,B(1),C(1),D(1)),y(u,v,B(2),C(2),D(2)),y(u,v,B(3),C(3),D(3)));
        end
        if (m12*nA(1)+m22*nA(2)+m23*nA(3)) ~= 0
            d=d+BCD/4*(m12*nA(1)+m22*nA(2)+m23*nA(3))*G22(y(u,v,B(1),C(1),D(1)),y(u,v,B(2),C(2),D(2)),y(u,v,B(3),C(3),D(3)));
        end
        if (m13*nA(1)+m23*nA(2)+m33*nA(3)) ~= 0
            d=d+BCD/4*(m13*nA(1)+m23*nA(2)+m33*nA(3))*G32(y(u,v,B(1),C(1),D(1)),y(u,v,B(2),C(2),D(2)),y(u,v,B(3),C(3),D(3)));
        end
        if (m11*nB(1)+m12*nB(2)+m13*nB(3)) ~= 0
            d=d+CDA/4*(m11*nB(1)+m12*nB(2)+m13*nB(3))*G12(y(u,v,C(1),D(1),A(1)),y(u,v,C(2),D(2),A(2)),y(u,v,C(3),D(3),A(3)));
        end
        if (m12*nB(1)+m22*nB(2)+m23*nB(3)) ~= 0
            d=d+CDA/4*(m12*nB(1)+m22*nB(2)+m23*nB(3))*G22(y(u,v,C(1),D(1),A(1)),y(u,v,C(2),D(2),A(2)),y(u,v,C(3),D(3),A(3)));
        end
        if (m13*nB(1)+m23*nB(2)+m33*nB(3)) ~= 0
            d=d+CDA/4*(m13*nB(1)+m23*nB(2)+m33*nB(3))*G32(y(u,v,C(1),D(1),A(1)),y(u,v,C(2),D(2),A(2)),y(u,v,C(3),D(3),A(3)));
        end
        if (m11*nC(1)+m12*nC(2)+m13*nC(3)) ~= 0
            d=d+DAB/4*(m11*nC(1)+m12*nC(2)+m13*nC(3))*G12(y(u,v,D(1),A(1),B(1)),y(u,v,D(2),A(2),B(2)),y(u,v,D(3),A(3),B(3)));
        end
        if (m12*nC(1)+m22*nC(2)+m23*nC(3)) ~= 0
            d=d+DAB/4*(m12*nC(1)+m22*nC(2)+m23*nC(3))*G22(y(u,v,D(1),A(1),B(1)),y(u,v,D(2),A(2),B(2)),y(u,v,D(3),A(3),B(3)));
        end
        if (m13*nC(1)+m23*nC(2)+m33*nC(3)) ~= 0
            d=d+DAB/4*(m13*nC(1)+m23*nC(2)+m33*nC(3))*G32(y(u,v,D(1),A(1),B(1)),y(u,v,D(2),A(2),B(2)),y(u,v,D(3),A(3),B(3)));
        end
    end

    function d = IU3(u,v)
        % function IU3 is the integrand for displacement component u3
        d=zeros(size(x3));
        if (m11*nD(1)+m12*nD(2)+m13*nD(3)) ~= 0
            d=d+ABC/4*(m11*nD(1)+m12*nD(2)+m13*nD(3))*G13(y(u,v,A(1),B(1),C(1)),y(u,v,A(2),B(2),C(2)),y(u,v,A(3),B(3),C(3)));
        end
        if (m12*nD(1)+m22*nD(2)+m23*nD(3)) ~= 0
            d=d+ABC/4*(m12*nD(1)+m22*nD(2)+m23*nD(3))*G23(y(u,v,A(1),B(1),C(1)),y(u,v,A(2),B(2),C(2)),y(u,v,A(3),B(3),C(3)));
        end
        if (m13*nD(1)+m23*nD(2)+m33*nD(3)) ~= 0
            d=d+ABC/4*(m13*nD(1)+m23*nD(2)+m33*nD(3))*G33(y(u,v,A(1),B(1),C(1)),y(u,v,A(2),B(2),C(2)),y(u,v,A(3),B(3),C(3)));
        end
        if (m11*nA(1)+m12*nA(2)+m13*nA(3)) ~= 0
            d=d+BCD/4*(m11*nA(1)+m12*nA(2)+m13*nA(3))*G13(y(u,v,B(1),C(1),D(1)),y(u,v,B(2),C(2),D(2)),y(u,v,B(3),C(3),D(3)));
        end
        if (m12*nA(1)+m22*nA(2)+m23*nA(3)) ~= 0
            d=d+BCD/4*(m12*nA(1)+m22*nA(2)+m23*nA(3))*G23(y(u,v,B(1),C(1),D(1)),y(u,v,B(2),C(2),D(2)),y(u,v,B(3),C(3),D(3)));
        end
        if (m13*nA(1)+m23*nA(2)+m33*nA(3)) ~= 0
            d=d+BCD/4*(m13*nA(1)+m23*nA(2)+m33*nA(3))*G33(y(u,v,B(1),C(1),D(1)),y(u,v,B(2),C(2),D(2)),y(u,v,B(3),C(3),D(3)));
        end
        if (m11*nB(1)+m12*nB(2)+m13*nB(3)) ~= 0
            d=d+CDA/4*(m11*nB(1)+m12*nB(2)+m13*nB(3))*G13(y(u,v,C(1),D(1),A(1)),y(u,v,C(2),D(2),A(2)),y(u,v,C(3),D(3),A(3)));
        end
        if (m12*nB(1)+m22*nB(2)+m23*nB(3)) ~= 0
            d=d+CDA/4*(m12*nB(1)+m22*nB(2)+m23*nB(3))*G23(y(u,v,C(1),D(1),A(1)),y(u,v,C(2),D(2),A(2)),y(u,v,C(3),D(3),A(3)));
        end
        if (m13*nB(1)+m23*nB(2)+m33*nB(3)) ~= 0
            d=d+CDA/4*(m13*nB(1)+m23*nB(2)+m33*nB(3))*G33(y(u,v,C(1),D(1),A(1)),y(u,v,C(2),D(2),A(2)),y(u,v,C(3),D(3),A(3)));
        end
        if (m11*nC(1)+m12*nC(2)+m13*nC(3)) ~= 0
            d=d+DAB/4*(m11*nC(1)+m12*nC(2)+m13*nC(3))*G13(y(u,v,D(1),A(1),B(1)),y(u,v,D(2),A(2),B(2)),y(u,v,D(3),A(3),B(3)));
        end
        if (m12*nC(1)+m22*nC(2)+m23*nC(3)) ~= 0
            d=d+DAB/4*(m12*nC(1)+m22*nC(2)+m23*nC(3))*G23(y(u,v,D(1),A(1),B(1)),y(u,v,D(2),A(2),B(2)),y(u,v,D(3),A(3),B(3)));
        end
        if (m13*nC(1)+m23*nC(2)+m33*nC(3)) ~= 0
            d=d+DAB/4*(m13*nC(1)+m23*nC(2)+m33*nC(3))*G33(y(u,v,D(1),A(1),B(1)),y(u,v,D(2),A(2),B(2)),y(u,v,D(3),A(3),B(3)));
        end
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
