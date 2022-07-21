function [N,dNdxi,dNdeta]=evaluateNodalBasisQua(XI,nodesCoord,degree)
% [N,dNdxi,dNdeta]=evaluateNodalBasisQua(XI,nodesCoord,degree)
% Evaluates at XI the nodal basis of polynomials for the given degree
% with nodal coodinates at nodesCoord

% Vandermonde matrix
V=orthogonalPolynomialsQua(degree,nodesCoord);

%Orthogonal basis at XIs and change of polynomial basis
[P,dPdxi,dPdeta]=orthogonalPolynomialsAndDerivativesQua(degree,XI);
N=P/V; dNdxi = dPdxi/V; dNdeta = dPdeta/V; 


