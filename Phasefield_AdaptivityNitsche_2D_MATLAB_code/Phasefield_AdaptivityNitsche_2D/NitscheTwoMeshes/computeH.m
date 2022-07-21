function H = computeH(X,T,Xref,Tref,refinedElements,referenceElementStd,referenceElementRef,ux,uy,H_previous,lambda,mu)

%H(:,iElem) = H evaluated in IP of element iElem

nOfElements = size(T,1);
nIPstd = length(referenceElementStd.IPweights);
nIPref = length(referenceElementRef.IPweights);

standardElements = setdiff(1:nOfElements,refinedElements);
nRefEls = length(refinedElements); nStdEls = nOfElements-nRefEls;

elasticEnergyDensity_positive = zeros(nStdEls*nIPstd+nRefEls*nIPref,1);

NxiStd = referenceElementStd.Nxi;
NetaStd = referenceElementStd.Neta;
NxiRef = referenceElementRef.Nxi;
NetaRef = referenceElementRef.Neta;

for i = 1:length(standardElements)
    iElem = standardElements(i);
    Te = T(iElem,:);
    Xe = X(Te,:);
    ux_e = ux(Te); uy_e = uy(Te);
    
    [J1,J2,J3,J4] = computeJ_IP(Xe,ux_e,uy_e,NxiStd,NetaStd);
    strain = [J1,J4,1/2*(J2+J3)]; %strain(i,:) = [eps_x, eps_y, eps_xy]
    trStrain = strain(:,1)+strain(:,2);
    detStrain = strain(:,1).*strain(:,2)-strain(:,3).^2;
    vap1 = trStrain/2 + sqrt(trStrain.^2/4-detStrain); vap2 = trStrain/2 - sqrt(trStrain.^2/4-detStrain);
    Apos = (1/2*(trStrain+abs(trStrain))).^2; Aneg = (1/2*(trStrain-abs(trStrain))).^2;
    Bpos = (1/2*(vap1+abs(vap1))).^2+(1/2*(vap2+abs(vap2))).^2; Bneg = (1/2*(vap1-abs(vap1))).^2+(1/2*(vap2-abs(vap2))).^2;
    elasticEnergyDensity_positive_e = 1/2*lambda*Apos+mu*Bpos;
    ind = (i-1)*nIPstd + (1:nIPstd);
    elasticEnergyDensity_positive(ind) = elasticEnergyDensity_positive_e;
end

for i = 1:length(refinedElements)
    Te = Tref(i,:); Xe = Xref(Te,:);
    ux_e = ux(Te+size(X,1)); uy_e = uy(Te+size(X,1));
    
    [J1,J2,J3,J4] = computeJ_IP(Xe,ux_e,uy_e,NxiRef,NetaRef);
    strain = [J1,J4,1/2*(J2+J3)]; %strain(i,:) = [eps_x, eps_y, eps_xy]
    trStrain = strain(:,1)+strain(:,2);
    detStrain = strain(:,1).*strain(:,2)-strain(:,3).^2;
    vap1 = trStrain/2 + sqrt(trStrain.^2/4-detStrain); vap2 = trStrain/2 - sqrt(trStrain.^2/4-detStrain);
    Apos = (1/2*(trStrain+abs(trStrain))).^2; Aneg = (1/2*(trStrain-abs(trStrain))).^2;
    Bpos = (1/2*(vap1+abs(vap1))).^2+(1/2*(vap2+abs(vap2))).^2; Bneg = (1/2*(vap1-abs(vap1))).^2+(1/2*(vap2-abs(vap2))).^2;
    elasticEnergyDensity_positive_e = 1/2*lambda*Apos+mu*Bpos;
    ind = nStdEls*nIPstd + (i-1)*nIPref + (1:nIPref);
    elasticEnergyDensity_positive(ind) = elasticEnergyDensity_positive_e;
end

H = max(H_previous,elasticEnergyDensity_positive);

end

function [J1,J2,J3,J4] = computeJ_IP(Xe,ux,uy,Nxi,Neta)
J11 = Nxi*Xe(:,1); J12 = Nxi*Xe(:,2);
J21 = Neta*Xe(:,1); J22 = Neta*Xe(:,2);
detJ = J11.*J22-J12.*J21;
%maybe we should use bsxfun instead of diagonal matrices...
invJ11 = diag(J22./detJ);
invJ12 = diag(-J12./detJ);
invJ21 = diag(-J21./detJ);
invJ22 = diag(J11./detJ);
% xy-derivatives
Nx = invJ11*Nxi + invJ12*Neta;
Ny = invJ21*Nxi + invJ22*Neta;

J1 = Nx*ux; J2 = Ny*ux;
J3 = Nx*uy; J4 = Ny*uy;
end