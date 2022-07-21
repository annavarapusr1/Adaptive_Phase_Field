function P=orthogonalPolynomialsQua(degree,XI)

nOfPoints=size(XI,1); nOfPolynomials=(degree+1)^2;
P=zeros(nOfPoints,nOfPolynomials);
for k=1:nOfPoints
    u=XI(k,1); v=XI(k,2);
    [Pu,dPudxi]=orthogonalPolynomialsAndDerivatives1D(degree,u);
    [Pv,dPvdxi]=orthogonalPolynomialsAndDerivatives1D(degree,v);
    P(k,:)=reshape(Pu'*Pv,1,nOfPolynomials);
end

