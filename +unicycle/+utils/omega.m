function y=omega(x,beta)

y=zeros(size(x));
pos=abs(x) < 1/(1-beta)/2;
y(pos)=cos(pi*((1-beta)*abs(x(pos))-0.5+beta)/2/beta).^2;
y(abs(x) < (1-2*beta)/(1-beta)/2)=1;
    
    
