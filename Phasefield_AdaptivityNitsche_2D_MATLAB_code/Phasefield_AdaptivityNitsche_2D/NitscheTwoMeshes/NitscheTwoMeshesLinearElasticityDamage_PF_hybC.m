function [K,f] = NitscheTwoMeshesLinearElasticityDamage_PF(X,T,Xref,Tref,NitscheFaces,refinedElements,referenceElementStd,referenceElementRef,E_elems,nu_elems,eta,beta,d,hybridCondition,hybridCondition_NitscheFaces)

nOfElements = size(T,1);
nOfNitscheFaces = size(NitscheFaces,1);
standardElements = setdiff(1:nOfElements,refinedElements);

nOfNodes = size(X,1)+size(Xref,1);
nDOF = 2*nOfNodes;
f = zeros(nDOF,1);

nIPstd = length(referenceElementStd.IPweights);
nIPref = length(referenceElementRef.IPweights);

%[mKe,mKe] = size(Ke)
mKe = 2*size(referenceElementStd.NodesCoord,1); nStdEls = length(standardElements);
ind_i_std  = zeros(mKe*mKe,nStdEls);
ind_j_std  = zeros(mKe*mKe,nStdEls);
coefK_std = zeros(mKe*mKe,nStdEls);

% Volume integrals (loop in elements)
for i = 1:length(standardElements)
    iElem = standardElements(i);
    E = E_elems(iElem); nu = nu_elems(iElem);
    Te = T(iElem,:); Xe = X(Te,:);
    de = d(Te);
    hC_e = hybridCondition((i-1)*nIPstd + (1:nIPstd));
    
    Ke = computeVolumeMatrices(Xe,referenceElementStd,E,nu,eta,de,hC_e);
    ind = [Te,Te+nOfNodes];
    [mi,mj] = meshgrid(ind,ind);
    ind_i_std(:,i) = mi(:); ind_j_std(:,i) = mj(:);
    coefK_std(:,i) = Ke(:);
end

mKe = 2*size(referenceElementRef.NodesCoord,1); nRefEls = length(refinedElements);
ind_i_ref  = zeros(mKe*mKe,nRefEls);
ind_j_ref  = zeros(mKe*mKe,nRefEls);
coefK_ref = zeros(mKe*mKe,nRefEls);

for i = 1:length(refinedElements)
    iElem = refinedElements(i);
    E = E_elems(iElem); nu = nu_elems(iElem);
    Te = Tref(i,:); Xe = Xref(Te,:);
    de = d(Te+size(X,1));
    hC_e = hybridCondition(nStdEls*nIPstd + (i-1)*nIPref + (1:nIPref));
    Ke = computeVolumeMatrices(Xe,referenceElementRef,E,nu,eta,de,hC_e);
    ind = [Te,Te+nOfNodes]+size(X,1);
    [mi,mj] = meshgrid(ind,ind);
    ind_i_ref(:,i) = mi(:); ind_j_ref(:,i) = mj(:);
    coefK_ref(:,i) = Ke(:);
end

K = sparse([ind_i_std(:);ind_i_ref(:)],[ind_j_std(:);ind_j_ref(:)],[coefK_std(:);coefK_ref(:)]);

%Integrals on Nitsche faces
%Calculo com si refinat per les dues bandes + projeccio a l'espai std on
%toqui
n = size(referenceElementRef.NfacesIP{1},2);
N2dfaces=referenceElementRef.NfacesIP; N2dxifaces=referenceElementRef.NxifacesIP; N2detafaces=referenceElementRef.NetafacesIP;
N2dxiGeofaces = referenceElementRef.NxiGeofacesIP; N2detaGeofaces = referenceElementRef.NetaGeofacesIP;
IPw_f = referenceElementRef.IPweights1d; ngf = length(IPw_f);
faceNodesRef = referenceElementRef.faceNodes;

for iFace = 1:nOfNitscheFaces
    %% Info de la cara
    infoiFace = NitscheFaces(iFace,:);
    refElem = infoiFace(2); refFaceNum = infoiFace(3);
    stdElem = infoiFace(4); stdFaceNum = infoiFace(5);
    
    hC_face = hybridCondition_NitscheFaces((iFace-1)*ngf + (1:ngf));
    
    %% Funcions aprox
    refN2d = N2dfaces{refFaceNum};
    stdN2d = N2dfaces{stdFaceNum}; stdN2d = stdN2d(end:-1:1,:); %flipping of integration points for 2nd element
    refN2dxi = N2dxifaces{refFaceNum};
    stdN2dxi = N2dxifaces{stdFaceNum}; stdN2dxi = stdN2dxi(end:-1:1,:);
    refN2deta = N2detafaces{refFaceNum};
    stdN2deta = N2detafaces{stdFaceNum}; stdN2deta = stdN2deta(end:-1:1,:);
    
    %% Funcions geo
    stdN2dxiGeo = N2dxiGeofaces{stdFaceNum}; stdN2dxiGeo = stdN2dxiGeo(end:-1:1,:);
    stdN2detaGeo = N2detaGeofaces{stdFaceNum}; stdN2detaGeo = stdN2detaGeo(end:-1:1,:);
    
    %% Std element derivatives
    iElem = stdElem;
    Nxi = stdN2dxi; Neta = stdN2deta;
    NxiGeo = stdN2dxiGeo; NetaGeo = stdN2detaGeo;
    Te = T(iElem,:); Xe = X(Te,:);  Testd = [Te,Te+nOfNodes];
    damage_std_e = d(Te);
    J11 = NxiGeo*Xe(:,1); J12 = NxiGeo*Xe(:,2);
    J21 = NetaGeo*Xe(:,1); J22 = NetaGeo*Xe(:,2);
    detJ = J11.*J22-J12.*J21;
    invJ11 = spdiags(J22./detJ,0,ngf,ngf); invJ12 = spdiags(-J12./detJ,0,ngf,ngf);
    invJ21 = spdiags(-J21./detJ,0,ngf,ngf); invJ22 = spdiags(J11./detJ,0,ngf,ngf);
    stdN2dx = invJ11*Nxi + invJ12*Neta; stdN2dy = invJ21*Nxi + invJ22*Neta;
    
    %% Ref element derivatives
    iElem = refElem; equivRef = find(refinedElements==iElem);
    Nxi = refN2dxi; Neta = refN2deta;
    Te = Tref(equivRef,:); Xe = Xref(Te,:);
    damage_ref_e = d(Te+size(X,1));
    Teref = [Te+size(X,1),Te+size(X,1)+nOfNodes];
    J11 = Nxi*Xe(:,1); J12 = Nxi*Xe(:,2); J21 = Neta*Xe(:,1); J22 = Neta*Xe(:,2); detJ = J11.*J22-J12.*J21;
    invJ11 = spdiags(J22./detJ,0,ngf,ngf); invJ12 = spdiags(-J12./detJ,0,ngf,ngf);
    invJ21 = spdiags(-J21./detJ,0,ngf,ngf); invJ22 = spdiags(J11./detJ,0,ngf,ngf);
    refN2dx = invJ11*Nxi + invJ12*Neta; refN2dy = invJ21*Nxi + invJ22*Neta;
    
    %% Integrals on the face (as seen from ref element)
    Xfref = Xref(Tref(equivRef,faceNodesRef(refFaceNum,:)),:); % Nodes on the face
    dxdxi = referenceElementRef.N1dxi*Xfref(:,1);
    dydxi = referenceElementRef.N1dxi*Xfref(:,2);
    
    % d als punts d'integració de la cara
    damage_ref_g = refN2d*damage_ref_e;
    NGeofacesIP = referenceElementRef.NGeofacesIP{stdFaceNum}; NGeofacesIP = NGeofacesIP(end:-1:1,:);
    damage_std_g = NGeofacesIP*damage_std_e;
    
    damage_ref_g(hC_face == 1) = 0;
    damage_std_g(hC_face == 1) = 0;
    
    dxdxiNorm = sqrt(dxdxi.^2+dydxi.^2);
    dline = dxdxiNorm.*IPw_f;
    nx = dydxi./dxdxiNorm; ny=-dxdxi./dxdxiNorm; %normal exterior to the refined element
    Ddline = spdiags(dline,0,ngf,ngf);
    
    Ddline_nx_dstd = spdiags(dline.*nx.*((1-damage_std_g).^2+eta),0,ngf,ngf);
    Ddline_ny_dstd = spdiags(dline.*ny.*((1-damage_std_g).^2+eta),0,ngf,ngf);
    Ddline_nx_dref = spdiags(dline.*nx.*((1-damage_ref_g).^2+eta),0,ngf,ngf);
    Ddline_ny_dref = spdiags(dline.*ny.*((1-damage_ref_g).^2+eta),0,ngf,ngf);
    
    Er = E_elems(refElem); Es = E_elems(stdElem);
    nur = nu_elems(refElem); nus = nu_elems(stdElem);
    Ar = Er/((1+nur)*(1-2*nur)); As = Es/((1+nus)*(1-2*nus));
    Z = zeros(n,n);
    KnxNx = refN2d'*Ddline_nx_dref*refN2dx; KnyNy = refN2d'*Ddline_ny_dref*refN2dy;
    KnxNy = refN2d'*Ddline_nx_dref*refN2dy; KnyNx = refN2d'*Ddline_ny_dref*refN2dx;
    KnxNxT = refN2dx'*Ddline_nx_dref*refN2d; KnyNyT = refN2dy'*Ddline_ny_dref*refN2d;
    KnxNyT = refN2dy'*Ddline_nx_dref*refN2d; KnyNxT = refN2dx'*Ddline_ny_dref*refN2d;
    KeRR = -0.5*Ar*[(1-nu)*KnxNx+(1-2*nu)/2*KnyNy,nu*KnxNy+(1-2*nu)/2*KnyNx; (1-2*nu)/2*KnxNy+nu*KnyNx, (1-2*nu)/2*KnxNx+(1-nu)*KnyNy]...
        -0.5*Ar*[(1-nu)*KnxNxT+(1-2*nu)/2*KnyNyT,nu*KnyNxT+(1-2*nu)/2*KnxNyT; (1-2*nu)/2*KnyNxT+nu*KnxNyT, (1-2*nu)/2*KnxNxT+(1-nu)*KnyNyT]...
        +beta*[refN2d'*(Ddline*refN2d),Z;Z,refN2d'*(Ddline*refN2d)];
    
    KnxNx = stdN2d'*Ddline_nx_dstd*stdN2dx; KnyNy = stdN2d'*Ddline_ny_dstd*stdN2dy;
    KnxNy = stdN2d'*Ddline_nx_dstd*stdN2dy; KnyNx = stdN2d'*Ddline_ny_dstd*stdN2dx;
    KnxNxT = stdN2dx'*Ddline_nx_dstd*stdN2d; KnyNyT = stdN2dy'*Ddline_ny_dstd*stdN2d;
    KnxNyT = stdN2dy'*Ddline_nx_dstd*stdN2d; KnyNxT = stdN2dx'*Ddline_ny_dstd*stdN2d;
    KeSS = +0.5*As*[(1-nu)*KnxNx+(1-2*nu)/2*KnyNy,nu*KnxNy+(1-2*nu)/2*KnyNx; (1-2*nu)/2*KnxNy+nu*KnyNx, (1-2*nu)/2*KnxNx+(1-nu)*KnyNy]...
        +0.5*As*[(1-nu)*KnxNxT+(1-2*nu)/2*KnyNyT,nu*KnyNxT+(1-2*nu)/2*KnxNyT; (1-2*nu)/2*KnyNxT+nu*KnxNyT, (1-2*nu)/2*KnxNxT+(1-nu)*KnyNyT]...
        +beta*[stdN2d'*(Ddline*stdN2d),Z;Z,stdN2d'*(Ddline*stdN2d)];
    
    KnxNx = refN2d'*Ddline_nx_dstd*stdN2dx; KnyNy = refN2d'*Ddline_ny_dstd*stdN2dy;
    KnxNy = refN2d'*Ddline_nx_dstd*stdN2dy; KnyNx = refN2d'*Ddline_ny_dstd*stdN2dx;
    KnxNxT = refN2dx'*Ddline_nx_dref*stdN2d; KnyNyT = refN2dy'*Ddline_ny_dref*stdN2d;
    KnxNyT = refN2dy'*Ddline_nx_dref*stdN2d; KnyNxT = refN2dx'*Ddline_ny_dref*stdN2d;
    KeRS = -0.5*As*[(1-nu)*KnxNx+(1-2*nu)/2*KnyNy,nu*KnxNy+(1-2*nu)/2*KnyNx; (1-2*nu)/2*KnxNy+nu*KnyNx, (1-2*nu)/2*KnxNx+(1-nu)*KnyNy]...
        +0.5*As*[(1-nu)*KnxNxT+(1-2*nu)/2*KnyNyT,nu*KnyNxT+(1-2*nu)/2*KnxNyT; (1-2*nu)/2*KnyNxT+nu*KnxNyT, (1-2*nu)/2*KnxNxT+(1-nu)*KnyNyT]...
        -beta*[refN2d'*(Ddline*stdN2d),Z;Z,refN2d'*(Ddline*stdN2d)];
    
    %     KnxNx = stdN2d'*Ddline_nx*refN2dx; KnyNy = stdN2d'*Ddline_ny*refN2dy;
    %     KnxNy = stdN2d'*Ddline_nx*refN2dy; KnyNx = stdN2d'*Ddline_ny*refN2dx;
    %     KnxNxT = stdN2dx'*Ddline_nx*refN2d; KnyNyT = stdN2dy'*Ddline_ny*refN2d;
    %     KnxNyT = stdN2dy'*Ddline_nx*refN2d; KnyNxT = stdN2dx'*Ddline_ny*refN2d;
    %     KeRL = +0.5*Ar*[(1-nu)*KnxNx+(1-2*nu)/2*KnyNy,nu*KnxNy+(1-2*nu)/2*KnyNx; (1-2*nu)/2*KnxNy+nu*KnyNx, (1-2*nu)/2*KnxNx+(1-nu)*KnyNy]...
    %         -0.5*Ar*[(1-nu)*KnxNxT+(1-2*nu)/2*KnyNyT,nu*KnyNxT+(1-2*nu)/2*KnxNyT; (1-2*nu)/2*KnyNxT+nu*KnxNyT, (1-2*nu)/2*KnxNxT+(1-nu)*KnyNyT]...
    %         -beta*[stdN2d'*(Ddline*refN2d),Z;Z,stdN2d'*(Ddline*refN2d)];
    
    %% Reduccio matrius en cada component
    P = referenceElementRef.P2d;
    Pxy = blkdiag(P,P);
    KeSS = Pxy*KeSS*Pxy';
    KeRS = KeRS*Pxy';
    %KeRL = P2*KeRL;
    
    %% Ensamblatge
    K(Teref,Teref)=K(Teref,Teref)+KeRR;
    K(Teref,Testd)=K(Teref,Testd)+KeRS;
    K(Testd,Teref)=K(Testd,Teref)+KeRS';
    K(Testd,Testd)=K(Testd,Testd)+KeSS;
end

function Ke = computeVolumeMatrices(Xe,referenceElement,E,nu,eta,de,hybridCondition)
N = referenceElement.N; Nxi = referenceElement.Nxi; Neta = referenceElement.Neta;
IPw = referenceElement.IPweights; ngauss = length(IPw);
J11 = Nxi*Xe(:,1); J12 = Nxi*Xe(:,2); J21 = Neta*Xe(:,1); J22 = Neta*Xe(:,2);
detJ = J11.*J22-J12.*J21;
%maybe we should use bsxfun instead of diagonal matrices...
invJ11 = spdiags(J22./detJ,0,ngauss,ngauss); invJ12 = spdiags(-J12./detJ,0,ngauss,ngauss);
invJ21 = spdiags(-J21./detJ,0,ngauss,ngauss); invJ22 = spdiags(J11./detJ,0,ngauss,ngauss);
Nx = invJ11*Nxi + invJ12*Neta; Ny = invJ21*Nxi + invJ22*Neta;

dg = N*de;
dg(hybridCondition == 1) = 0;

dvolu_d = spdiags(((1-dg).^2+eta).*(referenceElement.IPweights.*detJ),0,ngauss,ngauss);
Kxxe = Nx'*(dvolu_d*Nx);
Kyye = Ny'*(dvolu_d*Ny);
Kxye = Nx'*(dvolu_d*Ny);
Kyxe = Ny'*(dvolu_d*Nx);
Ke=E/((1+nu)*(1-2*nu))*[(1-nu)*Kxxe+1/2*(1-2*nu)*Kyye,nu*Kxye+(1-2*nu)/2*Kyxe;
    nu*Kyxe+(1-2*nu)/2*Kxye,(1-nu)*Kyye+1/2*(1-2*nu)*Kxxe];
