function p = orthopoly2D_rst(x,n)
% p = orthopoly2D_rst(x,n)
% Computes the ortogonal base of 2D polynomials of degree less 
% or equal to n at the point x=(r,s) in [-1,1]^2
%

N = (n+1)*(n+2)/2 ;%number of nodes/polynomials
p = zeros(N,1);

r = x(1); s = x(2); 

%Ordering: 1st incresing the degree and 2nd lexicogafic order
ncount = 0; %counter for the polynomials order
%Loop on degree
for nDeg = 0:n
  %Loop increasing i
  for i = 0:nDeg
     if i==0
         p_i = 1;  q_i = 1;  
     else
         p_i = jacobiP(i,0,0,r);     q_i = q_i*(1-s)/2;
     end
     %Value for j
     j = nDeg-i;
     if j==0
        p_j = 1;
     else
        p_j = jacobiP(j,2*i+1,0,s);
     end
     ncount = ncount+1;
     factor = sqrt( (2*i+1)*(i+j+1)/2 );
     p(ncount) = ( p_i*q_i*p_j )*factor;
  end
end