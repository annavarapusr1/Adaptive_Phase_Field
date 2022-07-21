function [K,f] = NitscheTwoMeshesSystemLinearElasticity(refinedElsGlobalNum,E_elements,nu_elements,X,T,Xref,Tref,refinedElements,referenceElementStd,referenceElementRef,infoFaces,beta,NitscheFaces,sourceTermFunction)

nOfElements = size(T,1);
nOfNitscheFaces = size(NitscheFaces,1);
nOfFaceNodes = size(referenceElementStd.NodesCoord1d,1);
faceNodes = referenceElementStd.faceNodes;

standardElements = setdiff(1:nOfElements,refinedElements);

nOfNodes = size(X,1)+size(Xref,1);
nDOF = 2*nOfNodes;
f = zeros(nDOF,1);
K = spalloc(nDOF,nDOF,nOfFaceNodes*nDOF);

% Volume integrals (loop in elements)
for i = 1:length(standardElements)
    iElem = standardElements(i);
    E = E_elements(iElem); nu = nu_elements(iElem);
    Te = T(iElem,:); Xe = X(Te,:);
    %figure(1),hold on, plot(Xe(:,1),Xe(:,2),'ko');
    [Ke,fe] = computeVolumeMatrices(Xe,referenceElementStd,E,nu,sourceTermFunction);
    %Assembly
    indRC = [Te,Te+nOfNodes];
    f(indRC) = f(indRC) + fe;
    K(indRC,indRC)=K(indRC,indRC)+ Ke;
end

for i = 1:length(refinedElements)
    iElem = refinedElements(i);
    E = E_elements(iElem); nu = nu_elements(iElem);
    equivRef = find(refinedElsGlobalNum==iElem);
    Te = Tref(equivRef,:); Xe = Xref(Te,:);
    %figure(1),hold on, plot(Xe(:,1),Xe(:,2),'ko','MarkerSize',8);
    
    [Ke,fe] = computeVolumeMatrices(Xe,referenceElementRef,E,nu,sourceTermFunction);
    %Assembly
    
    indRC = [Te,Te+nOfNodes]+size(X,1);
    f(indRC) = f(indRC) + fe;
    K(indRC,indRC)=K(indRC,indRC)+ Ke;
end

%Integrals on INTERIOR faces (loop in faces)
%Calculo com si refinat per les dues bandes + projeccio a l'espai std on
%toqui
n = size(referenceElementRef.NfacesIP,2);
N2dfaces=referenceElementRef.NfacesIP; N2dxifaces=referenceElementRef.NxifacesIP; N2detafaces=referenceElementRef.NetafacesIP;
N2dxiGeofaces = referenceElementRef.NxiGeofacesIP;
N2detaGeofaces = referenceElementRef.NetaGeofacesIP;
N1dxiGeo = referenceElementRef.N1dxiGeo;
IPw_f = referenceElementRef.IPweights1d; ngf = length(IPw_f);
for iFace=1:nOfNitscheFaces
    %% Info de la cara
    infoiFace = NitscheFaces(iFace,:);
    leftElem = infoiFace(2); leftFaceNum = infoiFace(3);
    rightElem = infoiFace(4); rightFaceNum = infoiFace(5);
    
    %% Funcions aprox
    refN2d = N2dfaces(:,:,leftFaceNum);
    stdN2d = N2dfaces(end:-1:1,:,rightFaceNum); %flipping of integration points for 2nd element
    refN2dxi = N2dxifaces(:,:,leftFaceNum); stdN2dxi = N2dxifaces(end:-1:1,:,rightFaceNum);
    refN2deta = N2detafaces(:,:,leftFaceNum); stdN2deta = N2detafaces(end:-1:1,:,rightFaceNum);
    
    %% Funcions geo
    refN2dxiGeo = N2dxiGeofaces(:,:,leftFaceNum);
    stdN2dxiGeo = N2dxiGeofaces(end:-1:1,:,rightFaceNum);
    refN2detaGeo = N2detaGeofaces(:,:,leftFaceNum);
    stdN2detaGeo = N2detaGeofaces(end:-1:1,:,rightFaceNum);
    
    %% Std element derivatives
    iElem = rightElem;
    Nxi = stdN2dxi; Neta = stdN2deta;
    NxiGeo = stdN2dxiGeo; NetaGeo = stdN2detaGeo;
    Te = T(iElem,:); Xe = X(Te,:);  Testd = [Te,Te+nOfNodes];
    %  figure(1), hold on, plot(Xe(:,1),Xe(:,2),'r*'), hold off
    J11 = NxiGeo*Xe(:,1); J12 = NxiGeo*Xe(:,2);
    J21 = NetaGeo*Xe(:,1); J22 = NetaGeo*Xe(:,2);
    detJ = J11.*J22-J12.*J21;
    invJ11 = spdiags(J22./detJ,0,ngf,ngf); invJ12 = spdiags(-J12./detJ,0,ngf,ngf);
    invJ21 = spdiags(-J21./detJ,0,ngf,ngf); invJ22 = spdiags(J11./detJ,0,ngf,ngf);
    stdN2dx = invJ11*Nxi + invJ12*Neta;    stdN2dy = invJ21*Nxi + invJ22*Neta;
    %% Ref element derivatives
    iElem = leftElem; equivRef = find(refinedElsGlobalNum==iElem);
    Nxi = refN2dxi; Neta = refN2deta;
    NxiGeo = refN2dxiGeo; NetaGeo = refN2detaGeo;
    Te = T(iElem,:); Xe = X(Te,:);    
    Teref = [Tref(equivRef,:)+size(X,1),Tref(equivRef,:)+size(X,1)+nOfNodes];
    J11 = NxiGeo*Xe(:,1); J12 = NxiGeo*Xe(:,2); J21 = NetaGeo*Xe(:,1); J22 = NetaGeo*Xe(:,2); detJ = J11.*J22-J12.*J21;
    invJ11 = spdiags(J22./detJ,0,ngf,ngf); invJ12 = spdiags(-J12./detJ,0,ngf,ngf);
    invJ21 = spdiags(-J21./detJ,0,ngf,ngf); invJ22 = spdiags(J11./detJ,0,ngf,ngf);
    refN2dx = invJ11*Nxi + invJ12*Neta;    refN2dy = invJ21*Nxi + invJ22*Neta;
    %Elemental matrices
    %% Integrals on the face (as seen from ref element)
    nodes = faceNodes(leftFaceNum,:); Xf = X(T(leftElem,nodes),:); % Nodes in the face
    %figure(1), hold on, plot(Xf(:,1),Xf(:,2),'gs','MarkerSize',15), hold off
    dxdxi = N1dxiGeo*Xf(:,1); dydxi = N1dxiGeo*Xf(:,2);
    dxdxiNorm = sqrt(dxdxi.^2+dydxi.^2);
    dline = dxdxiNorm.*IPw_f;
    nx = dydxi./dxdxiNorm; ny=-dxdxi./dxdxiNorm; %normal exterior to the left element
    Ddline = spdiags(dline,0,ngf,ngf); Ddline_nx = spdiags(dline.*nx,0,ngf,ngf); Ddline_ny = spdiags(dline.*ny,0,ngf,ngf);
    
    Er = E_elements(leftElem); Es = E_elements(rightElem);
    nur = nu_elements(leftElem); nus = nu_elements(rightElem);
    Ar = Er/((1+nur)*(1-2*nur)); As = Es/((1+nus)*(1-2*nus));
    Z = zeros(n,n);
    KnxNx = refN2d'*Ddline_nx*refN2dx; KnyNy = refN2d'*Ddline_ny*refN2dy;
    KnxNy = refN2d'*Ddline_nx*refN2dy; KnyNx = refN2d'*Ddline_ny*refN2dx;   
    KnxNxT = refN2dx'*Ddline_nx*refN2d; KnyNyT = refN2dy'*Ddline_ny*refN2d;
    KnxNyT = refN2dy'*Ddline_nx*refN2d; KnyNxT = refN2dx'*Ddline_ny*refN2d;    
    KeLL = -0.5*Ar*[(1-nu)*KnxNx+(1-2*nu)/2*KnyNy,nu*KnxNy+(1-2*nu)/2*KnyNx; (1-2*nu)/2*KnxNy+nu*KnyNx, (1-2*nu)/2*KnxNx+(1-nu)*KnyNy]...
        -0.5*Ar*[(1-nu)*KnxNxT+(1-2*nu)/2*KnyNyT,nu*KnyNxT+(1-2*nu)/2*KnxNyT; (1-2*nu)/2*KnyNxT+nu*KnxNyT, (1-2*nu)/2*KnxNxT+(1-nu)*KnyNyT]...
        +beta*[refN2d'*(Ddline*refN2d),Z;Z,refN2d'*(Ddline*refN2d)];
    
    KnxNx = stdN2d'*Ddline_nx*stdN2dx; KnyNy = stdN2d'*Ddline_ny*stdN2dy;
    KnxNy = stdN2d'*Ddline_nx*stdN2dy; KnyNx = stdN2d'*Ddline_ny*stdN2dx;   
    KnxNxT = stdN2dx'*Ddline_nx*stdN2d; KnyNyT = stdN2dy'*Ddline_ny*stdN2d;
    KnxNyT = stdN2dy'*Ddline_nx*stdN2d; KnyNxT = stdN2dx'*Ddline_ny*stdN2d;    
    KeRR = +0.5*As*[(1-nu)*KnxNx+(1-2*nu)/2*KnyNy,nu*KnxNy+(1-2*nu)/2*KnyNx; (1-2*nu)/2*KnxNy+nu*KnyNx, (1-2*nu)/2*KnxNx+(1-nu)*KnyNy]...
        +0.5*As*[(1-nu)*KnxNxT+(1-2*nu)/2*KnyNyT,nu*KnyNxT+(1-2*nu)/2*KnxNyT; (1-2*nu)/2*KnyNxT+nu*KnxNyT, (1-2*nu)/2*KnxNxT+(1-nu)*KnyNyT]...
        +beta*[stdN2d'*(Ddline*stdN2d),Z;Z,stdN2d'*(Ddline*stdN2d)];
    
    KnxNx = refN2d'*Ddline_nx*stdN2dx; KnyNy = refN2d'*Ddline_ny*stdN2dy;
    KnxNy = refN2d'*Ddline_nx*stdN2dy; KnyNx = refN2d'*Ddline_ny*stdN2dx;
    KnxNxT = refN2dx'*Ddline_nx*stdN2d; KnyNyT = refN2dy'*Ddline_ny*stdN2d;
    KnxNyT = refN2dy'*Ddline_nx*stdN2d; KnyNxT = refN2dx'*Ddline_ny*stdN2d;    
    KeLR = -0.5*As*[(1-nu)*KnxNx+(1-2*nu)/2*KnyNy,nu*KnxNy+(1-2*nu)/2*KnyNx; (1-2*nu)/2*KnxNy+nu*KnyNx, (1-2*nu)/2*KnxNx+(1-nu)*KnyNy]...
        +0.5*As*[(1-nu)*KnxNxT+(1-2*nu)/2*KnyNyT,nu*KnyNxT+(1-2*nu)/2*KnxNyT; (1-2*nu)/2*KnyNxT+nu*KnxNyT, (1-2*nu)/2*KnxNxT+(1-nu)*KnyNyT]...
        -beta*[refN2d'*(Ddline*stdN2d),Z;Z,refN2d'*(Ddline*stdN2d)];
    
%     KnxNx = stdN2d'*Ddline_nx*refN2dx; KnyNy = stdN2d'*Ddline_ny*refN2dy;
%     KnxNy = stdN2d'*Ddline_nx*refN2dy; KnyNx = stdN2d'*Ddline_ny*refN2dx;
%     KnxNxT = stdN2dx'*Ddline_nx*refN2d; KnyNyT = stdN2dy'*Ddline_ny*refN2d;
%     KnxNyT = stdN2dy'*Ddline_nx*refN2d; KnyNxT = stdN2dx'*Ddline_ny*refN2d;    
%     KeRL = +0.5*Ar*[(1-nu)*KnxNx+(1-2*nu)/2*KnyNy,nu*KnxNy+(1-2*nu)/2*KnyNx; (1-2*nu)/2*KnxNy+nu*KnyNx, (1-2*nu)/2*KnxNx+(1-nu)*KnyNy]...
%         -0.5*Ar*[(1-nu)*KnxNxT+(1-2*nu)/2*KnyNyT,nu*KnyNxT+(1-2*nu)/2*KnxNyT; (1-2*nu)/2*KnyNxT+nu*KnxNyT, (1-2*nu)/2*KnxNxT+(1-nu)*KnyNyT]...
%         -beta*[stdN2d'*(Ddline*refN2d),Z;Z,stdN2d'*(Ddline*refN2d)];
    
    
    %% Reduccio matrius (right = std)
    %% REDUCCIO EN CADA COMPONENT X I Y
    P = referenceElementRef.P2d;
    P2 = blkdiag(P,P);
    KeRR = P2*KeRR*P2';
    KeLR = KeLR*P2';
    %KeRL = P2*KeRL;
    
    %% Ensamblatge    
    K(Teref,Teref)=K(Teref,Teref)+KeLL;
    K(Teref,Testd)=K(Teref,Testd)+KeLR;
    K(Testd,Teref)=K(Testd,Teref)+KeLR';
    K(Testd,Testd)=K(Testd,Testd)+KeRR;
end

function [Ke,fe] = computeVolumeMatrices(Xe,referenceElement,E,nu,source)
N = referenceElement.N; Nxi = referenceElement.Nxi; Neta = referenceElement.Neta;
IPw = referenceElement.IPweights; ngauss = length(IPw);
% Jacobian
J11 = Nxi*Xe(:,1); J12 = Nxi*Xe(:,2); J21 = Neta*Xe(:,1); J22 = Neta*Xe(:,2);
detJ = J11.*J22-J12.*J21;
%maybe we should use bsxfun instead of diagonal matrices...
dvolu = spdiags(referenceElement.IPweights.*detJ,0,ngauss,ngauss);
% dvolu_d = spdiags(dvolu,0,ngauss,ngauss);
invJ11 = spdiags(J22./detJ,0,ngauss,ngauss); invJ12 = spdiags(-J12./detJ,0,ngauss,ngauss);
invJ21 = spdiags(-J21./detJ,0,ngauss,ngauss); invJ22 = spdiags(J11./detJ,0,ngauss,ngauss);
% xy-derivatives
Nx = invJ11*Nxi + invJ12*Neta;    Ny = invJ21*Nxi + invJ22*Neta;
%Computation of r.h.s. source term
Xg = N*Xe; %sourceTerm = -Gc*sourceTermFunction(Xg);

%Elemental matrices
Kxxe = Nx'*(dvolu*Nx);
Kyye = Ny'*(dvolu*Ny);
Kxye = Nx'*(dvolu*Ny);
Kyxe = Ny'*(dvolu*Nx);
Ke=E/((1+nu)*(1-2*nu))*[(1-nu)*Kxxe+1/2*(1-2*nu)*Kyye,nu*Kxye+(1-2*nu)/2*Kyxe;
    nu*Kyxe+(1-2*nu)/2*Kxye,(1-nu)*Kyye+1/2*(1-2*nu)*Kxxe];

test= 1;
b = source(Xg,E,nu,test);
fe = [N'*(dvolu*b(:,1)); N'*(dvolu*b(:,2))];