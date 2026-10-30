function [s11,s12,s13,s22,s23,s33]=computeStressTetrahedronShearZoneTanhSinh( ...
    x1,x2,x3,A,B,C,D,e11,e12,e13,e22,e23,e33,G,nu,varargin)
% function COMPUTESTRESSTETRAHEDRONSHEARZONETANHSINH computes the
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
% G, nu              rigidity and Poisson's ratio in the half space.
%
% Output:
% sij                stress components ij.
%
% Author: Sylvain Barbot (sbarbot@ntu.edu.sg) - May 19, 2018, Los Angeles.

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

% parameterized surface integral
y=@(u,v,A,B,C) A*(1-u)*(1-v)/4+B*(1+u)*(1-v)/4+C*(1+v)/2;

% moment density
m11=2*e11+lambda*ekk;
m12=2*e12;
m13=2*e13;
m22=2*e22+lambda*ekk;
m23=2*e23;
m33=2*e33+lambda*ekk;

%% numerical solution with tanh/sinh quadrature
h=optionStruct.precision;
n=fix(1/h*optionStruct.bound);

u11=zeros(size(x1));
u12=zeros(size(x1));
u13=zeros(size(x1));

u21=zeros(size(x1));
u22=zeros(size(x1));
u23=zeros(size(x1));

u31=zeros(size(x1));
u32=zeros(size(x1));
u33=zeros(size(x1));

for k=-n:n
    fprintf('k=%d/%d\n',k,n);
    wk=(0.5*h*pi*cosh(k*h))./(cosh(0.5*pi*sinh(k*h))).^2;
    xk=tanh(0.5*pi*sinh(k*h));
    for j=-n:n
        wj=(0.5*h*pi*cosh(j*h))./(cosh(0.5*pi*sinh(j*h))).^2;
        xj=tanh(0.5*pi*sinh(j*h));
        
        u11=u11+wk*wj*(1-xk)*IU11(xj,xk);
        u12=u12+wk*wj*(1-xk)*IU12(xj,xk);
        u13=u13+wk*wj*(1-xk)*IU13(xj,xk);
        
        u21=u21+wk*wj*(1-xk)*IU21(xj,xk);
        u22=u22+wk*wj*(1-xk)*IU22(xj,xk);
        u23=u23+wk*wj*(1-xk)*IU23(xj,xk);
        
        u31=u31+wk*wj*(1-xk)*IU31(xj,xk);
        u32=u32+wk*wj*(1-xk)*IU32(xj,xk);
        u33=u33+wk*wj*(1-xk)*IU33(xj,xk);
    end
end

% remove anelastic strain
omega=  heaviside(((A(1)+B(1)+C(1))/3-x1)*nD(1)+((A(2)+B(2)+C(2))/3-x2)*nD(2)+((A(3)+B(3)+C(3))/3-x3)*nD(3)) ...
      .*heaviside(((B(1)+C(1)+D(1))/3-x1)*nA(1)+((B(2)+C(2)+D(2))/3-x2)*nA(2)+((B(3)+C(3)+D(3))/3-x3)*nA(3)) ...
      .*heaviside(((C(1)+D(1)+A(1))/3-x1)*nB(1)+((C(2)+D(2)+A(2))/3-x2)*nB(2)+((C(3)+D(3)+A(3))/3-x3)*nB(3)) ...
      .*heaviside(((D(1)+A(1)+B(1))/3-x1)*nC(1)+((D(2)+A(2)+B(2))/3-x2)*nC(2)+((D(3)+A(3)+B(3))/3-x3)*nC(3));

% strain components
e11=u11        -omega*e11;
e12=(u12+u21)/2-omega*e12;
e13=(u13+u31)/2-omega*e13;
e22=u22        -omega*e22;
e23=(u23+u32)/2-omega*e23;
e33=u33        -omega*e33;

% divergence
div=e11+e22+e33;

s11=2*G*e11+G*lambda*div;
s12=2*G*e12;
s13=2*G*e13;
s22=2*G*e22+G*lambda*div;
s23=2*G*e23;
s33=2*G*e33+G*lambda*div;

    %% Green's functions
    function u=G11d1(y1,y2,y3)
        lr1=r1(y1,y2,y3);
        lr2=r2(y1,y2,y3);
        
        u=1/(16*pi*(1-nu))*(x1-y1).*( ...
            -(3-4*nu)./lr1.^3 ...
            -1./lr2.^3 ...
            +(2*lr1.^2-3*(x1-y1).^2)./lr1.^5 ...
            +(3-4*nu)*(2*lr2.^2-3*(x1-y1).^2)./lr2.^5 ...
            -6*y3*x3.*(3*lr2.^2-5*(x1-y1).^2)./lr2.^7 ...
            -4*(1-2*nu)*(1-nu)./(lr2.*(lr2+x3+y3).^2) ...
            -8*(1-2*nu)*(1-nu)./(lr2.*(lr2+x3+y3).^2) ...
            +4*(1-2*nu)*(1-nu)*(x1-y1).^2./(lr2.^3.*(lr2+x3+y3).^2) ...
            +8*(1-2*nu)*(1-nu)*(x1-y1).^2./(lr2.^2.*(lr2+x3+y3).^3) ...
            );
    end

    function u=G11d2(y1,y2,y3)
        lr1=r1(y1,y2,y3);
        lr2=r2(y1,y2,y3);
        
        u=1/(16*pi*(1-nu))*(x2-y2).*( ...
            -(3-4*nu)./lr1.^3 ...
            -1./lr2.^3 ...
            -3*(x1-y1).^2./lr1.^5 ...
            -3*(3-4*nu)*(x1-y1).^2./lr2.^5 ...
            -6*y3*x3.*(lr2.^2-5*(x1-y1).^2)./lr2.^7 ...
            -4*(1-2*nu)*(1-nu)./(lr2.*(lr2+x3+y3).^2) ...
            +4*(1-2*nu)*(1-nu)*(x1-y1).^2.*(3*lr2+x3+y3)./(lr2.^3.*(lr2+x3+y3).^3) ...
            );
    end

    function u=G11d3(y1,y2,y3)
        lr1=r1(y1,y2,y3);
        lr2=r2(y1,y2,y3);
        
        u=1/(16*pi*(1-nu))*( ...
            -(3-4*nu)*(x3-y3)./lr1.^3 ...
            -(x3+y3)./lr2.^3 ...
            -3*(x1-y1).^2.*(x3-y3)./lr1.^5 ...
            -3*(3-4*nu)*(x1-y1).^2.*(x3+y3)./lr2.^5 ...
            +2*y3*(lr2.^2-3*x3.*(x3+y3))./lr2.^5 ...
            -6*y3*(x1-y1).^2.*(lr2.^2-5*x3.*(x3+y3))./lr2.^7 ...
            -4*(1-2*nu)*(1-nu)./(lr2.*(lr2+x3+y3)) ...
            +4*(1-2*nu)*(1-nu)*(x1-y1).^2.*(2*lr2+x3+y3)./(lr2.^3.*(lr2+x3+y3).^2) ...
            );
    end

    function u=G21d1(y1,y2,y3)
        lr1=r1(y1,y2,y3);
        lr2=r2(y1,y2,y3);
        
        u=1/(16*pi*(1-nu))*(x2-y2).*( ...
            +(lr1.^2-3*(x1-y1).^2)./lr1.^5 ...
            +(3-4*nu)*(lr2.^2-3*(x1-y1).^2)./lr2.^5 ...
            -6*y3*x3.*(lr2.^2-5*(x1-y1).^2)./lr2.^7 ...
            -4*(1-2*nu)*(1-nu)./(lr2.*(lr2+x3+y3).^2) ...
            +4*(1-2*nu)*(1-nu)*(x1-y1).^2.*(3*lr2+x3+y3)./(lr2.^3.*(lr2+x3+y3).^3) ...
            );
    end

    function u=G21d2(y1,y2,y3)
        lr1=r1(y1,y2,y3);
        lr2=r2(y1,y2,y3);
        
        u=1/(16*pi*(1-nu))*(x1-y1).*( ...
            +(lr1.^2-3*(x2-y2).^2)./lr1.^5 ...
            +(3-4*nu)*(lr2.^2-3*(x2-y2).^2)./lr2.^5 ...
            -6*y3*x3.*(lr2.^2-5*(x2-y2).^2)./lr2.^7 ...
            -4*(1-2*nu)*(1-nu)./(lr2.*(lr2+x3+y3).^2) ...
            +4*(1-2*nu)*(1-nu)*(x2-y2).^2.*(3*lr2+x3+y3)./(lr2.^3.*(lr2+x3+y3).^3) ...
            );
    end

    function u=G21d3(y1,y2,y3)
        lr1=r1(y1,y2,y3);
        lr2=r2(y1,y2,y3);
        
        u=1/(16*pi*(1-nu))*(x1-y1).*(x2-y2).*( ...
            -3*(x3-y3)./lr1.^5 ...
            -3*(3-4*nu)*(x3+y3)./lr2.^5 ...
            -6*y3*(lr2.^2-5*x3.*(x3+y3))./lr2.^7 ...
            +4*(1-2*nu)*(1-nu)*(2*lr2+x3+y3)./(lr2.^3.*(lr2+x3+y3).^2) ...
            );
    end

    function u=G31d1(y1,y2,y3)
        lr1=r1(y1,y2,y3);
        lr2=r2(y1,y2,y3);
        
        u=1/(16*pi*(1-nu))*( ...
            +(x3-y3).*(lr1.^2-3*(x1-y1).^2)./lr1.^5 ...
            +(3-4*nu)*(x3-y3).*(lr2.^2-3*(x1-y1).^2)./lr2.^5 ...
            +6*y3*x3.*(x3+y3).*(lr2.^2-5*(x1-y1).^2)./lr2.^7 ...
            -4*(1-2*nu)*(1-nu)./(lr2.*(lr2+x3+y3)) ...
            +4*(1-2*nu)*(1-nu)*(x1-y1).^2.*(2*lr2+x3+y3)./(lr2.^3.*(lr2+x3+y3).^2) ...
            );
    end

    function u=G31d2(y1,y2,y3)
        lr1=r1(y1,y2,y3);
        lr2=r2(y1,y2,y3);
        
        u=1/(16*pi*(1-nu))*(x1-y1).*(x2-y2).*( ...
            -3*(x3-y3)./lr1.^5 ...
            -3*(3-4*nu)*(x3-y3)./lr2.^5 ...
            -30*y3*x3.*(x3+y3)./lr2.^7 ...
            +4*(1-2*nu)*(1-nu)*(2*lr2+x3+y3)./(lr2.^3.*(lr2+x3+y3).^2) ...
            );
    end

    function u=G31d3(y1,y2,y3)
        lr1=r1(y1,y2,y3);
        lr2=r2(y1,y2,y3);
        
        u=1/(16*pi*(1-nu))*(x1-y1).*( ...
            +(lr1.^2-3*(x3-y3).^2)./lr1.^5 ...
            +(3-4*nu)*(lr2.^2-3*(x3.^2-y3.^2))./lr2.^5 ...
            +6*y3*(2*x3+y3)./lr2.^5 ...
            -30*y3*x3.*(x3+y3).^2./lr2.^7 ...
            +4*(1-2*nu)*(1-nu)./lr2.^3 ...
            );
    end

    function u=G12d1(y1,y2,y3)
        lr1=r1(y1,y2,y3);
        lr2=r2(y1,y2,y3);
        
        u=1/(16*pi*(1-nu))*(x2-y2).*( ...
            +(lr1.^2-3*(x1-y1).^2)./lr1.^5 ...
            +(3-4*nu)*(lr2.^2-3*(x1-y1).^2)./lr2.^5 ...
            -6*y3*x3.*(lr2.^2-5*(x1-y1).^2)./lr2.^7 ...
            -4*(1-nu)*(1-2*nu)./(lr2.*(lr2+x3+y3).^2) ...
            +4*(1-nu)*(1-2*nu)*(x1-y1).^2.*(3*lr2+x3+y3)./(lr2.^3.*(lr2+x3+y3).^3) ...
            );
    end

    function u=G12d2(y1,y2,y3)
        lr1=r1(y1,y2,y3);
        lr2=r2(y1,y2,y3);
        
        u=1/(16*pi*(1-nu))*(x1-y1).*( ...
            +(lr1.^2-3*(x2-y2).^2)./lr1.^5 ...
            +(3-4*nu)*(lr2.^2-3*(x2-y2).^2)./lr2.^5 ...
            -6*y3*x3.*(lr2.^2-5*(x2-y2).^2)./lr2.^7 ...
            -4*(1-nu)*(1-2*nu)./(lr2.*(lr2+x3+y3).^2) ...
            +4*(1-nu)*(1-2*nu)*(x2-y2).^2.*(3*lr2+x3+y3)./(lr2.^3.*(lr2+x3+y3).^3) ...
            );
    end

    function u=G12d3(y1,y2,y3)
        lr1=r1(y1,y2,y3);
        lr2=r2(y1,y2,y3);
        
        u=1/(16*pi*(1-nu))*(x1-y1).*(x2-y2).*( ...
            -3*(x3-y3)./lr1.^5 ...
            -3*(3-4*nu)*(x3+y3)./lr2.^5 ...
            -6*y3*(lr2.^2-5*x3.*(x3+y3))./lr2.^7 ...
            +4*(1-2*nu)*(1-nu)*(2*lr2+x3+y3)./(lr2.^3.*(lr2+x3+y3).^2) ...
            );
    end

    function u=G22d1(y1,y2,y3)
        lr1=r1(y1,y2,y3);
        lr2=r2(y1,y2,y3);
        
        u=1/(16*pi*(1-nu))*(x1-y1).*( ...
            -(3-4*nu)./lr1.^3 ...
            -1./lr2.^3 ...
            -3*(x2-y2).^2./lr1.^5 ...
            -3*(3-4*nu)*(x2-y2).^2./lr2.^5 ...
            -6*y3*x3.*(lr2.^2-5*(x2-y2).^2)./lr2.^7 ...
            -4*(1-2*nu)*(1-nu)./(lr2.*(lr2+x3+y3).^2) ...
            +4*(1-2*nu)*(1-nu)*(x2-y2).^2.*(3*lr2+x3+y3)./(lr2.^3.*(lr2+x3+y3).^3) ...
            );
    end

    function u=G22d2(y1,y2,y3)
        lr1=r1(y1,y2,y3);
        lr2=r2(y1,y2,y3);
        
        u=1/(16*pi*(1-nu))*(x2-y2).*( ...
            -(3-4*nu)./lr1.^3 ...
            -1./lr2.^3 ...
            +(2*lr1.^2-3*(x2-y2).^2)./lr1.^5 ...
            +(3-4*nu)*(2*lr2.^2-3*(x2-y2).^2)./lr2.^5 ...
            -6*y3*x3.*(3*lr2.^2-5*(x2-y2).^2)./lr2.^7 ...
            -12*(1-2*nu)*(1-nu)./(lr2.*(lr2+x3+y3).^2) ...
            +4*(1-2*nu)*(1-nu)*(x2-y2).^2.*(3*lr2+x3+y3)./(lr2.^3.*(lr2+x3+y3).^3) ...
            );
    end

    function u=G22d3(y1,y2,y3)
        lr1=r1(y1,y2,y3);
        lr2=r2(y1,y2,y3);
        
        u=1/(16*pi*(1-nu))*( ...
            -(3-4*nu)*(x3-y3)./lr1.^3 ...
            -(x3+y3)./lr2.^3 ...
            -3*(x2-y2).^2.*(x3-y3)./lr1.^5 ...
            -3*(3-4*nu)*(x2-y2).^2.*(x3+y3)./lr2.^5 ...
            +2*y3*(lr2.^2-3*x3.*(x3+y3))./lr2.^5 ...
            -6*y3*(x2-y2).^2.*(lr2.^2-5*x3.*(x3+y3))./lr2.^7 ...
            -4*(1-2*nu)*(1-nu)./(lr2.*(lr2+x3+y3)) ...
            +4*(1-2*nu)*(1-nu)*(x2-y2).^2.*(2*lr2+x3+y3)./(lr2.^3.*(lr2+x3+y3).^2) ...
            );
    end

    function u=G32d1(y1,y2,y3)
        lr1=r1(y1,y2,y3);
        lr2=r2(y1,y2,y3);
        
        u=1/(16*pi*(1-nu))*(x1-y1).*(x2-y2).*( ...
            -3*(x3-y3)./lr1.^5 ...
            -3*(3-4*nu)*(x3-y3)./lr2.^5 ...
            -30*y3*x3.*(x3+y3)./lr2.^7 ...
            +4*(1-2*nu)*(1-nu)*(2*lr2+x3+y3)./(lr2.^3.*(lr2+x3+y3).^2) ...
            );
    end

    function u=G32d2(y1,y2,y3)
        lr1=r1(y1,y2,y3);
        lr2=r2(y1,y2,y3);
        
        u=1/(16*pi*(1-nu))*( ...
            +(x3-y3).*(lr1.^2-3*(x2-y2).^2)./lr1.^5 ...
            +(3-4*nu)*(x3-y3).*(lr2.^2-3*(x2-y2).^2)./lr2.^5 ...
            +6*y3*x3.*(x3+y3).*(lr2.^2-5*(x2-y2).^2)./lr2.^7 ...
            -4*(1-2*nu)*(1-nu)./(lr2.*(lr2+x3+y3)) ...
            +4*(1-2*nu)*(1-nu)*(x2-y2).^2.*(2*lr2+x3+y3)./(lr2.^3.*(lr2+x3+y3).^2) ...
            );
    end

    function u=G32d3(y1,y2,y3)
        lr1=r1(y1,y2,y3);
        lr2=r2(y1,y2,y3);
        
        u=1/(16*pi*(1-nu))*(x2-y2).*( ...
            +(lr1.^2-3*(x3-y3).^2)./lr1.^5 ...
            +(3-4*nu)*(lr2.^2-3*(x3.^2-y3^2))./lr2.^5 ...
            +6*y3*(2*x3+y3)./lr2.^5 ...
            -30*y3*x3.*(x3+y3).^2./lr2.^7 ...
            +4*(1-2*nu)*(1-nu)./lr2.^3 ...
            );
    end

    function u=G13d1(y1,y2,y3)
        lr1=r1(y1,y2,y3);
        lr2=r2(y1,y2,y3);
        
        u=1/(16*pi*(1-nu))*( ...
            +(x3-y3).*(lr1.^2-3*(x1-y1).^2)./lr1.^5 ...
            +(3-4*nu)*(x3-y3).*(lr2.^2-3*(x1-y1).^2)./lr2.^5 ...
            -6*y3*x3.*(x3+y3).*(lr2.^2-5*(x1-y1).^2)./lr2.^7 ...
            +4*(1-2*nu)*(1-nu)./(lr2.*(lr2+x3+y3)) ...
            -4*(1-2*nu)*(1-nu)*(x1-y1).^2.*(2*lr2+x3+y3)./(lr2.^3.*(lr2+x3+y3).^2) ...
            );
    end

    function u=G13d2(y1,y2,y3)
        lr1=r1(y1,y2,y3);
        lr2=r2(y1,y2,y3);
        
        u=1/(16*pi*(1-nu))*(x1-y1).*(x2-y2).*( ...
            -3*(x3-y3)./lr1.^5 ...
            -3*(3-4*nu)*(x3-y3)./lr2.^5 ...
            +30*y3*x3.*(x3+y3)./lr2.^7 ...
            -4*(1-2*nu)*(1-nu)*(2*lr2+x3+y3)./(lr2.^3.*(lr2+x3+y3).^2) ...
            );
    end

    function u=G13d3(y1,y2,y3)
        lr1=r1(y1,y2,y3);
        lr2=r2(y1,y2,y3);
        
        u=1/(16*pi*(1-nu))*(x1-y1).*( ...
            +(lr1.^2-3*(x3-y3).^2)./lr1.^5 ...
            +(3-4*nu)*(lr2.^2-3*(x3.^2-y3^2))./lr2.^5 ...
            -6*y3*(2*x3+y3)./lr2.^5 ...
            +30*y3*x3.*(x3+y3).^2./lr2.^7 ...
            -4*(1-2*nu)*(1-nu)./lr2.^3 ...
            );
    end

    function u=G23d1(y1,y2,y3)
        lr1=r1(y1,y2,y3);
        lr2=r2(y1,y2,y3);
        
        u=1/(16*pi*(1-nu))*(x1-y1).*(x2-y2).*( ...
            -3*(x3-y3)./lr1.^5 ...
            -3*(3-4*nu)*(x3-y3)./lr2.^5 ...
            +30*y3*x3.*(x3+y3)./lr2.^7 ...
            -4*(1-2*nu)*(1-nu)*(2*lr2+x3+y3)./(lr2.^3.*(lr2+x3+y3).^2) ...
            );
    end

    function u=G23d2(y1,y2,y3)
        lr1=r1(y1,y2,y3);
        lr2=r2(y1,y2,y3);
        
        u=1/(16*pi*(1-nu))*( ...
            +(x3-y3).*(lr1.^2-3*(x2-y2).^2)./lr1.^5 ...
            +(3-4*nu)*(x3-y3).*(lr2.^2-3*(x2-y2).^2)./lr2.^5 ...
            -6*y3*x3.*(x3+y3).*(lr2.^2-5*(x2-y2).^2)./lr2.^7 ...
            +4*(1-2*nu)*(1-nu)./(lr2.*(lr2+x3+y3)) ...
            -4*(1-2*nu)*(1-nu)*(x2-y2).^2.*(2*lr2+x3+y3)./(lr2.^3.*(lr2+x3+y3).^2) ...
            );
    end

    function u=G23d3(y1,y2,y3)
        lr1=r1(y1,y2,y3);
        lr2=r2(y1,y2,y3);
        
        u=1/(16*pi*(1-nu))*(x2-y2).*( ...
            +(lr1.^2-3*(x3-y3).^2)./lr1.^5 ...
            +(3-4*nu)*(lr2.^2-3*(x3.^2-y3^2))./lr2.^5 ...
            -6*y3*(2*x3+y3)./lr2.^5 ...
            +30*y3*x3.*(x3+y3).^2./lr2.^7 ...
            -4*(1-2*nu)*(1-nu)./lr2.^3 ...
            );
    end

    function u=G33d1(y1,y2,y3)
        lr1=r1(y1,y2,y3);
        lr2=r2(y1,y2,y3);
        
        u=1/(16*pi*(1-nu))*(x1-y1).*( ...
            -(3-4*nu)./lr1.^3 ...
            -(5-12*nu+8*nu^2)./lr2.^3 ...
            -3*(x3-y3).^2./lr1.^5 ...
            -30*y3*x3.*(x3+y3).^2./lr2.^7 ...
            -3*(3-4*nu)*(x3+y3).^2./lr2.^5 ...
            +6*y3*x3./lr2.^5 ...
            );
    end

    function u=G33d2(y1,y2,y3)
        lr1=r1(y1,y2,y3);
        lr2=r2(y1,y2,y3);
        
        u=1/(16*pi*(1-nu))*(x2-y2).*( ...
            -(3-4*nu)./lr1.^3 ...
            -(5-12*nu+8*nu^2)./lr2.^3 ...
            -3*(x3-y3).^2./lr1.^5 ...
            -30*y3*x3.*(x3+y3).^2./lr2.^7 ...
            -3*(3-4*nu)*(x3+y3).^2./lr2.^5 ...
            +6*y3*x3./lr2.^5 ...
            );
    end

    function u=G33d3(y1,y2,y3)
        lr1=r1(y1,y2,y3);
        lr2=r2(y1,y2,y3);
        
        u=1/(16*pi*(1-nu))*( ...
            -(3-4*nu)*(x3-y3)./lr1.^3 ...
            -(5-12*nu+8*nu^2)*(x3+y3)./lr2.^3 ...
            +(x3-y3).*(2*lr1.^2-3*(x3-y3).^2)./lr1.^5 ...
            +6*y3*(x3+y3).^2./lr2.^5 ...
            +6*y3*x3.*(x3+y3).*(2*lr2.^2-5*(x3+y3).^2)./lr2.^7 ...
            +(3-4*nu)*(x3+y3).*(2*lr2.^2-3*(x3+y3).^2)./lr2.^5 ...
            -2*y3*(lr2.^2-3*x3.*(x3+y3))./lr2.^5 ...
            );
    end

    function d = IU11(u,v)
        % function IU11 is the integrand for displacement gradient u1,1
        d=zeros(size(x1));
        if (m11*nD(1)+m12*nD(2)+m13*nD(3)) ~= 0
            d=d+ABC/4*(m11*nD(1)+m12*nD(2)+m13*nD(3))*G11d1(y(u,v,A(1),B(1),C(1)),y(u,v,A(2),B(2),C(2)),y(u,v,A(3),B(3),C(3)));
        end
        if (m12*nD(1)+m22*nD(2)+m23*nD(3)) ~= 0
            d=d+ABC/4*(m12*nD(1)+m22*nD(2)+m23*nD(3))*G21d1(y(u,v,A(1),B(1),C(1)),y(u,v,A(2),B(2),C(2)),y(u,v,A(3),B(3),C(3)));
        end
        if (m13*nD(1)+m23*nD(2)+m33*nD(3)) ~= 0
            d=d+ABC/4*(m13*nD(1)+m23*nD(2)+m33*nD(3))*G31d1(y(u,v,A(1),B(1),C(1)),y(u,v,A(2),B(2),C(2)),y(u,v,A(3),B(3),C(3)));
        end
        if (m11*nA(1)+m12*nA(2)+m13*nA(3)) ~= 0
            d=d+BCD/4*(m11*nA(1)+m12*nA(2)+m13*nA(3))*G11d1(y(u,v,B(1),C(1),D(1)),y(u,v,B(2),C(2),D(2)),y(u,v,B(3),C(3),D(3)));
        end
        if (m12*nA(1)+m22*nA(2)+m23*nA(3)) ~= 0
            d=d+BCD/4*(m12*nA(1)+m22*nA(2)+m23*nA(3))*G21d1(y(u,v,B(1),C(1),D(1)),y(u,v,B(2),C(2),D(2)),y(u,v,B(3),C(3),D(3)));
        end
        if (m13*nA(1)+m23*nA(2)+m33*nA(3)) ~= 0
            d=d+BCD/4*(m13*nA(1)+m23*nA(2)+m33*nA(3))*G31d1(y(u,v,B(1),C(1),D(1)),y(u,v,B(2),C(2),D(2)),y(u,v,B(3),C(3),D(3)));
        end
        if (m11*nB(1)+m12*nB(2)+m13*nB(3)) ~= 0
            d=d+CDA/4*(m11*nB(1)+m12*nB(2)+m13*nB(3))*G11d1(y(u,v,C(1),D(1),A(1)),y(u,v,C(2),D(2),A(2)),y(u,v,C(3),D(3),A(3)));
        end
        if (m12*nB(1)+m22*nB(2)+m23*nB(3)) ~= 0
            d=d+CDA/4*(m12*nB(1)+m22*nB(2)+m23*nB(3))*G21d1(y(u,v,C(1),D(1),A(1)),y(u,v,C(2),D(2),A(2)),y(u,v,C(3),D(3),A(3)));
        end
        if (m13*nB(1)+m23*nB(2)+m33*nB(3)) ~= 0
            d=d+CDA/4*(m13*nB(1)+m23*nB(2)+m33*nB(3))*G31d1(y(u,v,C(1),D(1),A(1)),y(u,v,C(2),D(2),A(2)),y(u,v,C(3),D(3),A(3)));
        end
        if (m11*nC(1)+m12*nC(2)+m13*nC(3)) ~= 0
            d=d+DAB/4*(m11*nC(1)+m12*nC(2)+m13*nC(3))*G11d1(y(u,v,D(1),A(1),B(1)),y(u,v,D(2),A(2),B(2)),y(u,v,D(3),A(3),B(3)));
        end
        if (m12*nC(1)+m22*nC(2)+m23*nC(3)) ~= 0
            d=d+DAB/4*(m12*nC(1)+m22*nC(2)+m23*nC(3))*G21d1(y(u,v,D(1),A(1),B(1)),y(u,v,D(2),A(2),B(2)),y(u,v,D(3),A(3),B(3)));
        end
        if (m13*nC(1)+m23*nC(2)+m33*nC(3)) ~= 0
            d=d+DAB/4*(m13*nC(1)+m23*nC(2)+m33*nC(3))*G31d1(y(u,v,D(1),A(1),B(1)),y(u,v,D(2),A(2),B(2)),y(u,v,D(3),A(3),B(3)));
        end
    end

    function d = IU12(u,v)
        % function IU12 is the integrand for displacement gradient u1,2
        d=zeros(size(x1));
        if (m11*nD(1)+m12*nD(2)+m13*nD(3)) ~= 0
            d=d+ABC/4*(m11*nD(1)+m12*nD(2)+m13*nD(3))*G11d2(y(u,v,A(1),B(1),C(1)),y(u,v,A(2),B(2),C(2)),y(u,v,A(3),B(3),C(3)));
        end
        if (m12*nD(1)+m22*nD(2)+m23*nD(3)) ~= 0
            d=d+ABC/4*(m12*nD(1)+m22*nD(2)+m23*nD(3))*G21d2(y(u,v,A(1),B(1),C(1)),y(u,v,A(2),B(2),C(2)),y(u,v,A(3),B(3),C(3)));
        end
        if (m13*nD(1)+m23*nD(2)+m33*nD(3)) ~= 0
            d=d+ABC/4*(m13*nD(1)+m23*nD(2)+m33*nD(3))*G31d2(y(u,v,A(1),B(1),C(1)),y(u,v,A(2),B(2),C(2)),y(u,v,A(3),B(3),C(3)));
        end
        if (m11*nA(1)+m12*nA(2)+m13*nA(3)) ~= 0
            d=d+BCD/4*(m11*nA(1)+m12*nA(2)+m13*nA(3))*G11d2(y(u,v,B(1),C(1),D(1)),y(u,v,B(2),C(2),D(2)),y(u,v,B(3),C(3),D(3)));
        end
        if (m12*nA(1)+m22*nA(2)+m23*nA(3)) ~= 0
            d=d+BCD/4*(m12*nA(1)+m22*nA(2)+m23*nA(3))*G21d2(y(u,v,B(1),C(1),D(1)),y(u,v,B(2),C(2),D(2)),y(u,v,B(3),C(3),D(3)));
        end
        if (m13*nA(1)+m23*nA(2)+m33*nA(3)) ~= 0
            d=d+BCD/4*(m13*nA(1)+m23*nA(2)+m33*nA(3))*G31d2(y(u,v,B(1),C(1),D(1)),y(u,v,B(2),C(2),D(2)),y(u,v,B(3),C(3),D(3)));
        end
        if (m11*nB(1)+m12*nB(2)+m13*nB(3)) ~= 0
            d=d+CDA/4*(m11*nB(1)+m12*nB(2)+m13*nB(3))*G11d2(y(u,v,C(1),D(1),A(1)),y(u,v,C(2),D(2),A(2)),y(u,v,C(3),D(3),A(3)));
        end
        if (m12*nB(1)+m22*nB(2)+m23*nB(3)) ~= 0
            d=d+CDA/4*(m12*nB(1)+m22*nB(2)+m23*nB(3))*G21d2(y(u,v,C(1),D(1),A(1)),y(u,v,C(2),D(2),A(2)),y(u,v,C(3),D(3),A(3)));
        end
        if (m13*nB(1)+m23*nB(2)+m33*nB(3)) ~= 0
            d=d+CDA/4*(m13*nB(1)+m23*nB(2)+m33*nB(3))*G31d2(y(u,v,C(1),D(1),A(1)),y(u,v,C(2),D(2),A(2)),y(u,v,C(3),D(3),A(3)));
        end
        if (m11*nC(1)+m12*nC(2)+m13*nC(3)) ~= 0
            d=d+DAB/4*(m11*nC(1)+m12*nC(2)+m13*nC(3))*G11d2(y(u,v,D(1),A(1),B(1)),y(u,v,D(2),A(2),B(2)),y(u,v,D(3),A(3),B(3)));
        end
        if (m12*nC(1)+m22*nC(2)+m23*nC(3)) ~= 0
            d=d+DAB/4*(m12*nC(1)+m22*nC(2)+m23*nC(3))*G21d2(y(u,v,D(1),A(1),B(1)),y(u,v,D(2),A(2),B(2)),y(u,v,D(3),A(3),B(3)));
        end
        if (m13*nC(1)+m23*nC(2)+m33*nC(3)) ~= 0
            d=d+DAB/4*(m13*nC(1)+m23*nC(2)+m33*nC(3))*G31d2(y(u,v,D(1),A(1),B(1)),y(u,v,D(2),A(2),B(2)),y(u,v,D(3),A(3),B(3)));
        end
    end

    function d = IU13(u,v)
        % function IU13 is the integrand for displacement gradient u1,3
        d=zeros(size(x1));
        if (m11*nD(1)+m12*nD(2)+m13*nD(3)) ~= 0
            d=d+ABC/4*(m11*nD(1)+m12*nD(2)+m13*nD(3))*G11d3(y(u,v,A(1),B(1),C(1)),y(u,v,A(2),B(2),C(2)),y(u,v,A(3),B(3),C(3)));
        end
        if (m12*nD(1)+m22*nD(2)+m23*nD(3)) ~= 0
            d=d+ABC/4*(m12*nD(1)+m22*nD(2)+m23*nD(3))*G21d3(y(u,v,A(1),B(1),C(1)),y(u,v,A(2),B(2),C(2)),y(u,v,A(3),B(3),C(3)));
        end
        if (m13*nD(1)+m23*nD(2)+m33*nD(3)) ~= 0
            d=d+ABC/4*(m13*nD(1)+m23*nD(2)+m33*nD(3))*G31d3(y(u,v,A(1),B(1),C(1)),y(u,v,A(2),B(2),C(2)),y(u,v,A(3),B(3),C(3)));
        end
        if (m11*nA(1)+m12*nA(2)+m13*nA(3)) ~= 0
            d=d+BCD/4*(m11*nA(1)+m12*nA(2)+m13*nA(3))*G11d3(y(u,v,B(1),C(1),D(1)),y(u,v,B(2),C(2),D(2)),y(u,v,B(3),C(3),D(3)));
        end
        if (m12*nA(1)+m22*nA(2)+m23*nA(3)) ~= 0
            d=d+BCD/4*(m12*nA(1)+m22*nA(2)+m23*nA(3))*G21d3(y(u,v,B(1),C(1),D(1)),y(u,v,B(2),C(2),D(2)),y(u,v,B(3),C(3),D(3)));
        end
        if (m13*nA(1)+m23*nA(2)+m33*nA(3)) ~= 0
            d=d+BCD/4*(m13*nA(1)+m23*nA(2)+m33*nA(3))*G31d3(y(u,v,B(1),C(1),D(1)),y(u,v,B(2),C(2),D(2)),y(u,v,B(3),C(3),D(3)));
        end
        if (m11*nB(1)+m12*nB(2)+m13*nB(3)) ~= 0
            d=d+CDA/4*(m11*nB(1)+m12*nB(2)+m13*nB(3))*G11d3(y(u,v,C(1),D(1),A(1)),y(u,v,C(2),D(2),A(2)),y(u,v,C(3),D(3),A(3)));
        end
        if (m12*nB(1)+m22*nB(2)+m23*nB(3)) ~= 0
            d=d+CDA/4*(m12*nB(1)+m22*nB(2)+m23*nB(3))*G21d3(y(u,v,C(1),D(1),A(1)),y(u,v,C(2),D(2),A(2)),y(u,v,C(3),D(3),A(3)));
        end
        if (m13*nB(1)+m23*nB(2)+m33*nB(3)) ~= 0
            d=d+CDA/4*(m13*nB(1)+m23*nB(2)+m33*nB(3))*G31d3(y(u,v,C(1),D(1),A(1)),y(u,v,C(2),D(2),A(2)),y(u,v,C(3),D(3),A(3)));
        end
        if (m11*nC(1)+m12*nC(2)+m13*nC(3)) ~= 0
            d=d+DAB/4*(m11*nC(1)+m12*nC(2)+m13*nC(3))*G11d3(y(u,v,D(1),A(1),B(1)),y(u,v,D(2),A(2),B(2)),y(u,v,D(3),A(3),B(3)));
        end
        if (m12*nC(1)+m22*nC(2)+m23*nC(3)) ~= 0
            d=d+DAB/4*(m12*nC(1)+m22*nC(2)+m23*nC(3))*G21d3(y(u,v,D(1),A(1),B(1)),y(u,v,D(2),A(2),B(2)),y(u,v,D(3),A(3),B(3)));
        end
        if (m13*nC(1)+m23*nC(2)+m33*nC(3)) ~= 0
            d=d+DAB/4*(m13*nC(1)+m23*nC(2)+m33*nC(3))*G31d3(y(u,v,D(1),A(1),B(1)),y(u,v,D(2),A(2),B(2)),y(u,v,D(3),A(3),B(3)));
        end
    end

    function d = IU21(u,v)
        % function IU21 is the integrand for displacement gradient u2,1
        d=zeros(size(x2));
        if (m11*nD(1)+m12*nD(2)+m13*nD(3)) ~= 0
            d=d+ABC/4*(m11*nD(1)+m12*nD(2)+m13*nD(3))*G12d1(y(u,v,A(1),B(1),C(1)),y(u,v,A(2),B(2),C(2)),y(u,v,A(3),B(3),C(3)));
        end
        if (m12*nD(1)+m22*nD(2)+m23*nD(3)) ~= 0
            d=d+ABC/4*(m12*nD(1)+m22*nD(2)+m23*nD(3))*G22d1(y(u,v,A(1),B(1),C(1)),y(u,v,A(2),B(2),C(2)),y(u,v,A(3),B(3),C(3)));
        end
        if (m13*nD(1)+m23*nD(2)+m33*nD(3)) ~= 0
            d=d+ABC/4*(m13*nD(1)+m23*nD(2)+m33*nD(3))*G32d1(y(u,v,A(1),B(1),C(1)),y(u,v,A(2),B(2),C(2)),y(u,v,A(3),B(3),C(3)));
        end
        if (m11*nA(1)+m12*nA(2)+m13*nA(3)) ~= 0
            d=d+BCD/4*(m11*nA(1)+m12*nA(2)+m13*nA(3))*G12d1(y(u,v,B(1),C(1),D(1)),y(u,v,B(2),C(2),D(2)),y(u,v,B(3),C(3),D(3)));
        end
        if (m12*nA(1)+m22*nA(2)+m23*nA(3)) ~= 0
            d=d+BCD/4*(m12*nA(1)+m22*nA(2)+m23*nA(3))*G22d1(y(u,v,B(1),C(1),D(1)),y(u,v,B(2),C(2),D(2)),y(u,v,B(3),C(3),D(3)));
        end
        if (m13*nA(1)+m23*nA(2)+m33*nA(3)) ~= 0
            d=d+BCD/4*(m13*nA(1)+m23*nA(2)+m33*nA(3))*G32d1(y(u,v,B(1),C(1),D(1)),y(u,v,B(2),C(2),D(2)),y(u,v,B(3),C(3),D(3)));
        end
        if (m11*nB(1)+m12*nB(2)+m13*nB(3)) ~= 0
            d=d+CDA/4*(m11*nB(1)+m12*nB(2)+m13*nB(3))*G12d1(y(u,v,C(1),D(1),A(1)),y(u,v,C(2),D(2),A(2)),y(u,v,C(3),D(3),A(3)));
        end
        if (m12*nB(1)+m22*nB(2)+m23*nB(3)) ~= 0
            d=d+CDA/4*(m12*nB(1)+m22*nB(2)+m23*nB(3))*G22d1(y(u,v,C(1),D(1),A(1)),y(u,v,C(2),D(2),A(2)),y(u,v,C(3),D(3),A(3)));
        end
        if (m13*nB(1)+m23*nB(2)+m33*nB(3)) ~= 0
            d=d+CDA/4*(m13*nB(1)+m23*nB(2)+m33*nB(3))*G32d1(y(u,v,C(1),D(1),A(1)),y(u,v,C(2),D(2),A(2)),y(u,v,C(3),D(3),A(3)));
        end
        if (m11*nC(1)+m12*nC(2)+m13*nC(3)) ~= 0
            d=d+DAB/4*(m11*nC(1)+m12*nC(2)+m13*nC(3))*G12d1(y(u,v,D(1),A(1),B(1)),y(u,v,D(2),A(2),B(2)),y(u,v,D(3),A(3),B(3)));
        end
        if (m12*nC(1)+m22*nC(2)+m23*nC(3)) ~= 0
            d=d+DAB/4*(m12*nC(1)+m22*nC(2)+m23*nC(3))*G22d1(y(u,v,D(1),A(1),B(1)),y(u,v,D(2),A(2),B(2)),y(u,v,D(3),A(3),B(3)));
        end
        if (m13*nC(1)+m23*nC(2)+m33*nC(3)) ~= 0
            d=d+DAB/4*(m13*nC(1)+m23*nC(2)+m33*nC(3))*G32d1(y(u,v,D(1),A(1),B(1)),y(u,v,D(2),A(2),B(2)),y(u,v,D(3),A(3),B(3)));
        end
    end

    function d = IU22(u,v)
        % function IU22 is the integrand for displacement gradient u2,2
        d=zeros(size(x2));
        if (m11*nD(1)+m12*nD(2)+m13*nD(3)) ~= 0
            d=d+ABC/4*(m11*nD(1)+m12*nD(2)+m13*nD(3))*G12d2(y(u,v,A(1),B(1),C(1)),y(u,v,A(2),B(2),C(2)),y(u,v,A(3),B(3),C(3)));
        end
        if (m12*nD(1)+m22*nD(2)+m23*nD(3)) ~= 0
            d=d+ABC/4*(m12*nD(1)+m22*nD(2)+m23*nD(3))*G22d2(y(u,v,A(1),B(1),C(1)),y(u,v,A(2),B(2),C(2)),y(u,v,A(3),B(3),C(3)));
        end
        if (m13*nD(1)+m23*nD(2)+m33*nD(3)) ~= 0
            d=d+ABC/4*(m13*nD(1)+m23*nD(2)+m33*nD(3))*G32d2(y(u,v,A(1),B(1),C(1)),y(u,v,A(2),B(2),C(2)),y(u,v,A(3),B(3),C(3)));
        end
        if (m11*nA(1)+m12*nA(2)+m13*nA(3)) ~= 0
            d=d+BCD/4*(m11*nA(1)+m12*nA(2)+m13*nA(3))*G12d2(y(u,v,B(1),C(1),D(1)),y(u,v,B(2),C(2),D(2)),y(u,v,B(3),C(3),D(3)));
        end
        if (m12*nA(1)+m22*nA(2)+m23*nA(3)) ~= 0
            d=d+BCD/4*(m12*nA(1)+m22*nA(2)+m23*nA(3))*G22d2(y(u,v,B(1),C(1),D(1)),y(u,v,B(2),C(2),D(2)),y(u,v,B(3),C(3),D(3)));
        end
        if (m13*nA(1)+m23*nA(2)+m33*nA(3)) ~= 0
            d=d+BCD/4*(m13*nA(1)+m23*nA(2)+m33*nA(3))*G32d2(y(u,v,B(1),C(1),D(1)),y(u,v,B(2),C(2),D(2)),y(u,v,B(3),C(3),D(3)));
        end
        if (m11*nB(1)+m12*nB(2)+m13*nB(3)) ~= 0
            d=d+CDA/4*(m11*nB(1)+m12*nB(2)+m13*nB(3))*G12d2(y(u,v,C(1),D(1),A(1)),y(u,v,C(2),D(2),A(2)),y(u,v,C(3),D(3),A(3)));
        end
        if (m12*nB(1)+m22*nB(2)+m23*nB(3)) ~= 0
            d=d+CDA/4*(m12*nB(1)+m22*nB(2)+m23*nB(3))*G22d2(y(u,v,C(1),D(1),A(1)),y(u,v,C(2),D(2),A(2)),y(u,v,C(3),D(3),A(3)));
        end
        if (m13*nB(1)+m23*nB(2)+m33*nB(3)) ~= 0
            d=d+CDA/4*(m13*nB(1)+m23*nB(2)+m33*nB(3))*G32d2(y(u,v,C(1),D(1),A(1)),y(u,v,C(2),D(2),A(2)),y(u,v,C(3),D(3),A(3)));
        end
        if (m11*nC(1)+m12*nC(2)+m13*nC(3)) ~= 0
            d=d+DAB/4*(m11*nC(1)+m12*nC(2)+m13*nC(3))*G12d2(y(u,v,D(1),A(1),B(1)),y(u,v,D(2),A(2),B(2)),y(u,v,D(3),A(3),B(3)));
        end
        if (m12*nC(1)+m22*nC(2)+m23*nC(3)) ~= 0
            d=d+DAB/4*(m12*nC(1)+m22*nC(2)+m23*nC(3))*G22d2(y(u,v,D(1),A(1),B(1)),y(u,v,D(2),A(2),B(2)),y(u,v,D(3),A(3),B(3)));
        end
        if (m13*nC(1)+m23*nC(2)+m33*nC(3)) ~= 0
            d=d+DAB/4*(m13*nC(1)+m23*nC(2)+m33*nC(3))*G32d2(y(u,v,D(1),A(1),B(1)),y(u,v,D(2),A(2),B(2)),y(u,v,D(3),A(3),B(3)));
        end
    end

    function d = IU23(u,v)
        % function IU23 is the integrand for displacement gradient u2,3
        d=zeros(size(x2));
        if (m11*nD(1)+m12*nD(2)+m13*nD(3)) ~= 0
            d=d+ABC/4*(m11*nD(1)+m12*nD(2)+m13*nD(3))*G12d3(y(u,v,A(1),B(1),C(1)),y(u,v,A(2),B(2),C(2)),y(u,v,A(3),B(3),C(3)));
        end
        if (m12*nD(1)+m22*nD(2)+m23*nD(3)) ~= 0
            d=d+ABC/4*(m12*nD(1)+m22*nD(2)+m23*nD(3))*G22d3(y(u,v,A(1),B(1),C(1)),y(u,v,A(2),B(2),C(2)),y(u,v,A(3),B(3),C(3)));
        end
        if (m13*nD(1)+m23*nD(2)+m33*nD(3)) ~= 0
            d=d+ABC/4*(m13*nD(1)+m23*nD(2)+m33*nD(3))*G32d3(y(u,v,A(1),B(1),C(1)),y(u,v,A(2),B(2),C(2)),y(u,v,A(3),B(3),C(3)));
        end
        if (m11*nA(1)+m12*nA(2)+m13*nA(3)) ~= 0
            d=d+BCD/4*(m11*nA(1)+m12*nA(2)+m13*nA(3))*G12d3(y(u,v,B(1),C(1),D(1)),y(u,v,B(2),C(2),D(2)),y(u,v,B(3),C(3),D(3)));
        end
        if (m12*nA(1)+m22*nA(2)+m23*nA(3)) ~= 0
            d=d+BCD/4*(m12*nA(1)+m22*nA(2)+m23*nA(3))*G22d3(y(u,v,B(1),C(1),D(1)),y(u,v,B(2),C(2),D(2)),y(u,v,B(3),C(3),D(3)));
        end
        if (m13*nA(1)+m23*nA(2)+m33*nA(3)) ~= 0
            d=d+BCD/4*(m13*nA(1)+m23*nA(2)+m33*nA(3))*G32d3(y(u,v,B(1),C(1),D(1)),y(u,v,B(2),C(2),D(2)),y(u,v,B(3),C(3),D(3)));
        end
        if (m11*nB(1)+m12*nB(2)+m13*nB(3)) ~= 0
            d=d+CDA/4*(m11*nB(1)+m12*nB(2)+m13*nB(3))*G12d3(y(u,v,C(1),D(1),A(1)),y(u,v,C(2),D(2),A(2)),y(u,v,C(3),D(3),A(3)));
        end
        if (m12*nB(1)+m22*nB(2)+m23*nB(3)) ~= 0
            d=d+CDA/4*(m12*nB(1)+m22*nB(2)+m23*nB(3))*G22d3(y(u,v,C(1),D(1),A(1)),y(u,v,C(2),D(2),A(2)),y(u,v,C(3),D(3),A(3)));
        end
        if (m13*nB(1)+m23*nB(2)+m33*nB(3)) ~= 0
            d=d+CDA/4*(m13*nB(1)+m23*nB(2)+m33*nB(3))*G32d3(y(u,v,C(1),D(1),A(1)),y(u,v,C(2),D(2),A(2)),y(u,v,C(3),D(3),A(3)));
        end
        if (m11*nC(1)+m12*nC(2)+m13*nC(3)) ~= 0
            d=d+DAB/4*(m11*nC(1)+m12*nC(2)+m13*nC(3))*G12d3(y(u,v,D(1),A(1),B(1)),y(u,v,D(2),A(2),B(2)),y(u,v,D(3),A(3),B(3)));
        end
        if (m12*nC(1)+m22*nC(2)+m23*nC(3)) ~= 0
            d=d+DAB/4*(m12*nC(1)+m22*nC(2)+m23*nC(3))*G22d3(y(u,v,D(1),A(1),B(1)),y(u,v,D(2),A(2),B(2)),y(u,v,D(3),A(3),B(3)));
        end
        if (m13*nC(1)+m23*nC(2)+m33*nC(3)) ~= 0
            d=d+DAB/4*(m13*nC(1)+m23*nC(2)+m33*nC(3))*G32d3(y(u,v,D(1),A(1),B(1)),y(u,v,D(2),A(2),B(2)),y(u,v,D(3),A(3),B(3)));
        end
    end

    function d = IU31(u,v)
        % function IU31 is the integrand for displacement gradient u3,1
        d=zeros(size(x3));
        if (m11*nD(1)+m12*nD(2)+m13*nD(3)) ~= 0
            d=d+ABC/4*(m11*nD(1)+m12*nD(2)+m13*nD(3))*G13d1(y(u,v,A(1),B(1),C(1)),y(u,v,A(2),B(2),C(2)),y(u,v,A(3),B(3),C(3)));
        end
        if (m12*nD(1)+m22*nD(2)+m23*nD(3)) ~= 0
            d=d+ABC/4*(m12*nD(1)+m22*nD(2)+m23*nD(3))*G23d1(y(u,v,A(1),B(1),C(1)),y(u,v,A(2),B(2),C(2)),y(u,v,A(3),B(3),C(3)));
        end
        if (m13*nD(1)+m23*nD(2)+m33*nD(3)) ~= 0
            d=d+ABC/4*(m13*nD(1)+m23*nD(2)+m33*nD(3))*G33d1(y(u,v,A(1),B(1),C(1)),y(u,v,A(2),B(2),C(2)),y(u,v,A(3),B(3),C(3)));
        end
        if (m11*nA(1)+m12*nA(2)+m13*nA(3)) ~= 0
            d=d+BCD/4*(m11*nA(1)+m12*nA(2)+m13*nA(3))*G13d1(y(u,v,B(1),C(1),D(1)),y(u,v,B(2),C(2),D(2)),y(u,v,B(3),C(3),D(3)));
        end
        if (m12*nA(1)+m22*nA(2)+m23*nA(3)) ~= 0
            d=d+BCD/4*(m12*nA(1)+m22*nA(2)+m23*nA(3))*G23d1(y(u,v,B(1),C(1),D(1)),y(u,v,B(2),C(2),D(2)),y(u,v,B(3),C(3),D(3)));
        end
        if (m13*nA(1)+m23*nA(2)+m33*nA(3)) ~= 0
            d=d+BCD/4*(m13*nA(1)+m23*nA(2)+m33*nA(3))*G33d1(y(u,v,B(1),C(1),D(1)),y(u,v,B(2),C(2),D(2)),y(u,v,B(3),C(3),D(3)));
        end
        if (m11*nB(1)+m12*nB(2)+m13*nB(3)) ~= 0
            d=d+CDA/4*(m11*nB(1)+m12*nB(2)+m13*nB(3))*G13d1(y(u,v,C(1),D(1),A(1)),y(u,v,C(2),D(2),A(2)),y(u,v,C(3),D(3),A(3)));
        end
        if (m12*nB(1)+m22*nB(2)+m23*nB(3)) ~= 0
            d=d+CDA/4*(m12*nB(1)+m22*nB(2)+m23*nB(3))*G23d1(y(u,v,C(1),D(1),A(1)),y(u,v,C(2),D(2),A(2)),y(u,v,C(3),D(3),A(3)));
        end
        if (m13*nB(1)+m23*nB(2)+m33*nB(3)) ~= 0
            d=d+CDA/4*(m13*nB(1)+m23*nB(2)+m33*nB(3))*G33d1(y(u,v,C(1),D(1),A(1)),y(u,v,C(2),D(2),A(2)),y(u,v,C(3),D(3),A(3)));
        end
        if (m11*nC(1)+m12*nC(2)+m13*nC(3)) ~= 0
            d=d+DAB/4*(m11*nC(1)+m12*nC(2)+m13*nC(3))*G13d1(y(u,v,D(1),A(1),B(1)),y(u,v,D(2),A(2),B(2)),y(u,v,D(3),A(3),B(3)));
        end
        if (m12*nC(1)+m22*nC(2)+m23*nC(3)) ~= 0
            d=d+DAB/4*(m12*nC(1)+m22*nC(2)+m23*nC(3))*G23d1(y(u,v,D(1),A(1),B(1)),y(u,v,D(2),A(2),B(2)),y(u,v,D(3),A(3),B(3)));
        end
        if (m13*nC(1)+m23*nC(2)+m33*nC(3)) ~= 0
            d=d+DAB/4*(m13*nC(1)+m23*nC(2)+m33*nC(3))*G33d1(y(u,v,D(1),A(1),B(1)),y(u,v,D(2),A(2),B(2)),y(u,v,D(3),A(3),B(3)));
        end
    end

    function d = IU32(u,v)
        % function IU32 is the integrand for displacement gradient u3,2
        d=zeros(size(x3));
        if (m11*nD(1)+m12*nD(2)+m13*nD(3)) ~= 0
            d=d+ABC/4*(m11*nD(1)+m12*nD(2)+m13*nD(3))*G13d2(y(u,v,A(1),B(1),C(1)),y(u,v,A(2),B(2),C(2)),y(u,v,A(3),B(3),C(3)));
        end
        if (m12*nD(1)+m22*nD(2)+m23*nD(3)) ~= 0
            d=d+ABC/4*(m12*nD(1)+m22*nD(2)+m23*nD(3))*G23d2(y(u,v,A(1),B(1),C(1)),y(u,v,A(2),B(2),C(2)),y(u,v,A(3),B(3),C(3)));
        end
        if (m13*nD(1)+m23*nD(2)+m33*nD(3)) ~= 0
            d=d+ABC/4*(m13*nD(1)+m23*nD(2)+m33*nD(3))*G33d2(y(u,v,A(1),B(1),C(1)),y(u,v,A(2),B(2),C(2)),y(u,v,A(3),B(3),C(3)));
        end
        if (m11*nA(1)+m12*nA(2)+m13*nA(3)) ~= 0
            d=d+BCD/4*(m11*nA(1)+m12*nA(2)+m13*nA(3))*G13d2(y(u,v,B(1),C(1),D(1)),y(u,v,B(2),C(2),D(2)),y(u,v,B(3),C(3),D(3)));
        end
        if (m12*nA(1)+m22*nA(2)+m23*nA(3)) ~= 0
            d=d+BCD/4*(m12*nA(1)+m22*nA(2)+m23*nA(3))*G23d2(y(u,v,B(1),C(1),D(1)),y(u,v,B(2),C(2),D(2)),y(u,v,B(3),C(3),D(3)));
        end
        if (m13*nA(1)+m23*nA(2)+m33*nA(3)) ~= 0
            d=d+BCD/4*(m13*nA(1)+m23*nA(2)+m33*nA(3))*G33d2(y(u,v,B(1),C(1),D(1)),y(u,v,B(2),C(2),D(2)),y(u,v,B(3),C(3),D(3)));
        end
        if (m11*nB(1)+m12*nB(2)+m13*nB(3)) ~= 0
            d=d+CDA/4*(m11*nB(1)+m12*nB(2)+m13*nB(3))*G13d2(y(u,v,C(1),D(1),A(1)),y(u,v,C(2),D(2),A(2)),y(u,v,C(3),D(3),A(3)));
        end
        if (m12*nB(1)+m22*nB(2)+m23*nB(3)) ~= 0
            d=d+CDA/4*(m12*nB(1)+m22*nB(2)+m23*nB(3))*G23d2(y(u,v,C(1),D(1),A(1)),y(u,v,C(2),D(2),A(2)),y(u,v,C(3),D(3),A(3)));
        end
        if (m13*nB(1)+m23*nB(2)+m33*nB(3)) ~= 0
            d=d+CDA/4*(m13*nB(1)+m23*nB(2)+m33*nB(3))*G33d2(y(u,v,C(1),D(1),A(1)),y(u,v,C(2),D(2),A(2)),y(u,v,C(3),D(3),A(3)));
        end
        if (m11*nC(1)+m12*nC(2)+m13*nC(3)) ~= 0
            d=d+DAB/4*(m11*nC(1)+m12*nC(2)+m13*nC(3))*G13d2(y(u,v,D(1),A(1),B(1)),y(u,v,D(2),A(2),B(2)),y(u,v,D(3),A(3),B(3)));
        end
        if (m12*nC(1)+m22*nC(2)+m23*nC(3)) ~= 0
            d=d+DAB/4*(m12*nC(1)+m22*nC(2)+m23*nC(3))*G23d2(y(u,v,D(1),A(1),B(1)),y(u,v,D(2),A(2),B(2)),y(u,v,D(3),A(3),B(3)));
        end
        if (m13*nC(1)+m23*nC(2)+m33*nC(3)) ~= 0
            d=d+DAB/4*(m13*nC(1)+m23*nC(2)+m33*nC(3))*G33d2(y(u,v,D(1),A(1),B(1)),y(u,v,D(2),A(2),B(2)),y(u,v,D(3),A(3),B(3)));
        end
    end

    function d = IU33(u,v)
        % function IU33 is the integrand for displacement gradient u3,3
        d=zeros(size(x3));
        if (m11*nD(1)+m12*nD(2)+m13*nD(3)) ~= 0
            d=d+ABC/4*(m11*nD(1)+m12*nD(2)+m13*nD(3))*G13d3(y(u,v,A(1),B(1),C(1)),y(u,v,A(2),B(2),C(2)),y(u,v,A(3),B(3),C(3)));
        end
        if (m12*nD(1)+m22*nD(2)+m23*nD(3)) ~= 0
            d=d+ABC/4*(m12*nD(1)+m22*nD(2)+m23*nD(3))*G23d3(y(u,v,A(1),B(1),C(1)),y(u,v,A(2),B(2),C(2)),y(u,v,A(3),B(3),C(3)));
        end
        if (m13*nD(1)+m23*nD(2)+m33*nD(3)) ~= 0
            d=d+ABC/4*(m13*nD(1)+m23*nD(2)+m33*nD(3))*G33d3(y(u,v,A(1),B(1),C(1)),y(u,v,A(2),B(2),C(2)),y(u,v,A(3),B(3),C(3)));
        end
        if (m11*nA(1)+m12*nA(2)+m13*nA(3)) ~= 0
            d=d+BCD/4*(m11*nA(1)+m12*nA(2)+m13*nA(3))*G13d3(y(u,v,B(1),C(1),D(1)),y(u,v,B(2),C(2),D(2)),y(u,v,B(3),C(3),D(3)));
        end
        if (m12*nA(1)+m22*nA(2)+m23*nA(3)) ~= 0
            d=d+BCD/4*(m12*nA(1)+m22*nA(2)+m23*nA(3))*G23d3(y(u,v,B(1),C(1),D(1)),y(u,v,B(2),C(2),D(2)),y(u,v,B(3),C(3),D(3)));
        end
        if (m13*nA(1)+m23*nA(2)+m33*nA(3)) ~= 0
            d=d+BCD/4*(m13*nA(1)+m23*nA(2)+m33*nA(3))*G33d3(y(u,v,B(1),C(1),D(1)),y(u,v,B(2),C(2),D(2)),y(u,v,B(3),C(3),D(3)));
        end
        if (m11*nB(1)+m12*nB(2)+m13*nB(3)) ~= 0
            d=d+CDA/4*(m11*nB(1)+m12*nB(2)+m13*nB(3))*G13d3(y(u,v,C(1),D(1),A(1)),y(u,v,C(2),D(2),A(2)),y(u,v,C(3),D(3),A(3)));
        end
        if (m12*nB(1)+m22*nB(2)+m23*nB(3)) ~= 0
            d=d+CDA/4*(m12*nB(1)+m22*nB(2)+m23*nB(3))*G23d3(y(u,v,C(1),D(1),A(1)),y(u,v,C(2),D(2),A(2)),y(u,v,C(3),D(3),A(3)));
        end
        if (m13*nB(1)+m23*nB(2)+m33*nB(3)) ~= 0
            d=d+CDA/4*(m13*nB(1)+m23*nB(2)+m33*nB(3))*G33d3(y(u,v,C(1),D(1),A(1)),y(u,v,C(2),D(2),A(2)),y(u,v,C(3),D(3),A(3)));
        end
        if (m11*nC(1)+m12*nC(2)+m13*nC(3)) ~= 0
            d=d+DAB/4*(m11*nC(1)+m12*nC(2)+m13*nC(3))*G13d3(y(u,v,D(1),A(1),B(1)),y(u,v,D(2),A(2),B(2)),y(u,v,D(3),A(3),B(3)));
        end
        if (m12*nC(1)+m22*nC(2)+m23*nC(3)) ~= 0
            d=d+DAB/4*(m12*nC(1)+m22*nC(2)+m23*nC(3))*G23d3(y(u,v,D(1),A(1),B(1)),y(u,v,D(2),A(2),B(2)),y(u,v,D(3),A(3),B(3)));
        end
        if (m13*nC(1)+m23*nC(2)+m33*nC(3)) ~= 0
            d=d+DAB/4*(m13*nC(1)+m23*nC(2)+m33*nC(3))*G33d3(y(u,v,D(1),A(1),B(1)),y(u,v,D(2),A(2),B(2)),y(u,v,D(3),A(3),B(3)));
        end
    end

end

function y=heaviside(x)
y=x>0;
end

function p = validatePrecision(x)
if ~(isreal(x) && isscalar(x) && x > 0)
    error(message('MATLAB:invalidPrecision','invalid precision'));
end
p = true;
end

function p = validateBound(x)
if ~(isreal(x) && isscalar(x) && x > 0)
    error(message('MATLAB:invalidBound','invalid trunction'));
end
p = true;
end

