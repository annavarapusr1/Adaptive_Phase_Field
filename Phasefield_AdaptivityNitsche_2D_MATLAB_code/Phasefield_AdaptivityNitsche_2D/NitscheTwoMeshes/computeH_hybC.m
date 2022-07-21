function [H,hybridCondition,hybridCondition_NitscheFaces] = computeH(X,T,Xref,Tref,refinedElements,referenceElementStd,referenceElementRef,ux,uy,H_previous,lambda,mu,NitscheFaces)

%H(:,iElem) = H evaluated in IP of element iElem

nOfElements = size(T,1);
nIPstd = length(referenceElementStd.IPweights);
nIPref = length(referenceElementRef.IPweights);

standardElements = setdiff(1:nOfElements,refinedElements);
nRefEls = length(refinedElements); nStdEls = nOfElements-nRefEls;

elasticEnergyDensity_positive = zeros(nStdEls*nIPstd+nRefEls*nIPref,1);
hybridCondition = zeros(nStdEls*nIPstd+nRefEls*nIPref,1);
nIP1dref = length(referenceElementRef.IPweights1d);
hybridCondition_NitscheFaces = zeros(size(NitscheFaces,1)*nIP1dref,1);

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
    elasticEnergyDensity_negative_e = 1/2*lambda*Aneg+mu*Bneg;
    
    hybridCondition(ind) = (elasticEnergyDensity_positive_e < elasticEnergyDensity_negative_e);
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
    elasticEnergyDensity_negative_e = 1/2*lambda*Aneg+mu*Bneg;

    hybridCondition(ind) = (elasticEnergyDensity_positive_e < elasticEnergyDensity_negative_e);
end

H = max(H_previous,elasticEnergyDensity_positive);


N2dfaces=referenceElementRef.NfacesIP; N2dxifaces=referenceElementRef.NxifacesIP; N2detafaces=referenceElementRef.NetafacesIP;

for iFace = 1:size(NitscheFaces,1)
    %% Info de la cara
    infoiFace = NitscheFaces(iFace,:);
    refElem = infoiFace(2); refFaceNum = infoiFace(3);

    %% Funcions aprox
    refN2d = N2dfaces{refFaceNum};
    refN2dxi = N2dxifaces{refFaceNum};
    refN2deta = N2detafaces{refFaceNum};

    %% Ref element derivatives
    iElem = refElem; equivRef = find(refinedElements==iElem);
    Nxi = refN2dxi; Neta = refN2deta;
    Te = Tref(equivRef,:); Xe = Xref(Te,:);
    ux_ref_e = ux(Te+size(X,1)); uy_ref_e = uy(Te+size(X,1));
    
    [J1,J2,J3,J4] = computeJ_IP(Xe,ux_ref_e,uy_ref_e,Nxi,Neta);
    
    strain = [J1,J4,1/2*(J2+J3)]; %strain(i,:) = [eps_x, eps_y, eps_xy]
    trStrain = strain(:,1)+strain(:,2);
    detStrain = strain(:,1).*strain(:,2)-strain(:,3).^2;
    vap1 = trStrain/2 + sqrt(trStrain.^2/4-detStrain); vap2 = trStrain/2 - sqrt(trStrain.^2/4-detStrain);
    Apos = (1/2*(trStrain+abs(trStrain))).^2; Aneg = (1/2*(trStrain-abs(trStrain))).^2;
    Bpos = (1/2*(vap1+abs(vap1))).^2+(1/2*(vap2+abs(vap2))).^2; Bneg = (1/2*(vap1-abs(vap1))).^2+(1/2*(vap2-abs(vap2))).^2;
    elasticEnergyDensity_positive_e = 1/2*lambda*Apos+mu*Bpos;
    elasticEnergyDensity_negative_e = 1/2*lambda*Aneg+mu*Bneg;

    ind = (iFace-1)*nIP1dref + (1:nIP1dref);
    hybridCondition_NitscheFaces(ind) = (elasticEnergyDensity_positive_e < elasticEnergyDensity_negative_e);
end



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