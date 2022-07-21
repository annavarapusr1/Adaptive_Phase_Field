function [p,dp_dxi,dp_deta] = orthopoly2D_deriv_rst(x,n)
%
% [p,dp_dxi,dp_deta] = orthopoly2D_deriv_rst(x,n)
% Computes the ortogonal base of 2D polynomials of degree less 
% or equal to n at the point x=(r,s) in [-1,1]^2
%

N = (n+1)*(n+2)/2 ;%number of nodes/polynomials
p = zeros(N,1);
dp_dxi  = zeros(N,1);
dp_deta = zeros(N,1);

r = x(1); s = x(2); 

xi = (1+r)*(1-s)/2-1;
eta = s;

dr_dxi  = 2/(1-eta);
dr_deta = 2*(1+xi)/(1-eta)^2;

%Ordering: 1st incresing the degree and 2nd lexicogafic order
ncount = 0; %counter for the polynomials order
%Loop on degree
for nDeg = 0:n
  %Loop increasing i
  for i = 0:nDeg
     if i==0
         p_i = 1;  q_i = 1;  dp_i = 0; dq_i = 0;
     else
         p_i = jacobiP(i,0,0,r); dp_i = jacobiP(i-1,1,1,r)*(i+1)/2;    
         q_i = q_i*(1-s)/2; dq_i = q_i*(-i)/(1-s);
     end
     %Value for j
     j = nDeg-i;
     if j==0
        p_j = 1;  dp_j = 0; 
     else
        p_j = jacobiP(j,2*i+1,0,s); dp_j = jacobiP(j-1,2*i+2,1,s)*(j+2*i+2)/2;  
     end
     ncount= ncount+1;
     factor = sqrt( (2*i+1)*(i+j+1)/2 );
     %Normalized polinomial
     p(ncount)    = ( p_i*q_i*p_j )*factor;
     %Derivatives with respect to (r,s)
     dp_dr = ( (dp_i)*q_i*p_j )*factor; 
     dp_ds = ( p_i*(dq_i*p_j+q_i*dp_j) )*factor;
     %Derivatives with respect to (xi,eta)
     dp_dxi(ncount)  = dp_dr*dr_dxi;
     dp_deta(ncount) = dp_dr*dr_deta + dp_ds;
  end
end