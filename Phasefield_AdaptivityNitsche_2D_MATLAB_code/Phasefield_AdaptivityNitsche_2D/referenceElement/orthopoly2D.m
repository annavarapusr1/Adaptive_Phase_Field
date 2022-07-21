function p = orthopoly2D(x,n)
% p = orthopoly2D(x,n)
% Computes the ortogonal base of 2D polynomials of degree less 
% or equal to n at the point x=(xi,eta) in the reference triangle
%

xi = x(1); eta = x(2); 

if eta==1 
    r = -1; s=1;
else
    r = 2*(1+xi)/(1-eta)-1;
    s = eta;
end

p = orthopoly2D_rst([r,s],n);