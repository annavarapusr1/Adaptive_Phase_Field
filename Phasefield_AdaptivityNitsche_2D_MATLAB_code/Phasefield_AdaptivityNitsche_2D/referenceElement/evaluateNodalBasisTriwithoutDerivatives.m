function N=evaluateNodalBasisTriwithoutDerivatives(XI,XInodes,degree)
% N=evaluateNodalBasisTriwithoutDerivatives(XIs,XInodes,degree)
% Evaluates at XI the nodal basis of polynomials for the given degree
% and the nodes in XInodes.

% Vandermonde matrix
V=orthogonalPolynomialsTri(degree,XInodes);
invV = inv(V);

%Orthogonal basis at XIs and change of polynomial basis
[P,dPdxi,dPdeta]=orthogonalPolynomialsAndDerivativesTri(degree,XI);
N=P*invV; 


