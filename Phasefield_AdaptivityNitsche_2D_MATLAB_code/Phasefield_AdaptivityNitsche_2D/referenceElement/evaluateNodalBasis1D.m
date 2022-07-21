function [N,dNdxi]=evaluateNodalBasis1D(Xi,XiNodes,degree)

[V,kk]=orthogonalPolynomialsAndDerivatives1D(degree,XiNodes);

[P,dPdxi]=orthogonalPolynomialsAndDerivatives1D(degree,Xi);

N = P/V;
dNdxi = dPdxi/V;

