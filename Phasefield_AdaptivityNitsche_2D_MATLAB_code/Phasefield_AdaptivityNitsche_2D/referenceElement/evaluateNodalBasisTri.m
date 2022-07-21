function [N,dNdxi,dNdeta]=evaluateNodalBasisTri(XI,XInodes,degree)
% [N,dNdxi,dNdeta]evaluateNodalBasisTri(XIs,XInodes,degree)
% Evaluates at XI the nodal basis of polynomials for the given degree
% and the nodes in XInodes.

% Vandermonde matrix
V=orthogonalPolynomialsTri(degree,XInodes);
%invV = inv(V);

%Orthogonal basis at XIs and change of polynomial basis
[P,dPdxi,dPdeta]=orthogonalPolynomialsAndDerivativesTri(degree,XI);
N=P/V; dNdxi = dPdxi/V; dNdeta = dPdeta/V;


