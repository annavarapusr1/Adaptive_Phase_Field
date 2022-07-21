function [K,f] = NitscheTwoMeshesSystemLaplace(EquivalenciaNumeracio,muElem,X,T,Xref,Tref,refinedElements,referenceElementStd,referenceElementRef,infoFaces,beta,NitscheFaces,sourceTermFunction)

nOfElements = size(T,1);
nOfNitscheFaces = size(NitscheFaces,1);
nOfInteriorFaces = size(infoFaces.intFaces,1);
nOfExteriorFaces = size(infoFaces.extFaces,1);
nOfFacesElem = size(referenceElementStd.faceNodes,1);
nOfFaceNodes = size(referenceElementStd.NodesCoord1d,1);
nOfElementNodes = size(referenceElementStd.NodesCoord,1);
faceNodes = referenceElementStd.faceNodes;

standardElements = setdiff(1:nOfElements,refinedElements);

nDOF = size(X,1)+size(Xref,1);
f = zeros(nDOF,1);
K = spalloc(nDOF,nDOF,nOfFaceNodes*nDOF);
%n = nOfElements*(nOfElementNodes)^2; indK = 1:nOfElementNodes; ind_i=zeros(1,n); ind_j=ind_i; coef_K=ind_i;

% Volume integrals (loop in elements)
N = referenceElementStd.N; Nxi = referenceElementStd.Nxi; Neta = referenceElementStd.Neta;
IPw = referenceElementStd.IPweights; ngauss = length(IPw);
for i = 1:length(standardElements)
    iElem = standardElements(i);
    mu=muElem(iElem);
    Te = T(iElem,:); Xe = X(Te,:);
    figure(1),hold on, plot(Xe(:,1),Xe(:,2),'k*');
    % Jacobian
    J11 = Nxi*Xe(:,1); J12 = Nxi*Xe(:,2); J21 = Neta*Xe(:,1); J22 = Neta*Xe(:,2);
    detJ = J11.*J22-J12.*J21;
    %maybe we should use bsxfun instead of diagonal matrices...
    dvolu = spdiags(referenceElementStd.IPweights.*detJ,0,ngauss,ngauss);
    invJ11 = spdiags(J22./detJ,0,ngauss,ngauss); invJ12 = spdiags(-J12./detJ,0,ngauss,ngauss);
    invJ21 = spdiags(-J21./detJ,0,ngauss,ngauss); invJ22 = spdiags(J11./detJ,0,ngauss,ngauss);
    % xy-derivatives
    Nx = invJ11*Nxi + invJ12*Neta;    Ny = invJ21*Nxi + invJ22*Neta;
    %Computation of r.h.s. source term
    Xg = N*Xe; sourceTerm = -mu*sourceTermFunction(Xg);
    fe = N'*(dvolu*sourceTerm);
    %Elemental matrix
    Ke = Nx'*(dvolu*Nx)+Ny'*(dvolu*Ny);
    %Assembly
    %indRC = (iElem-1)*nOfElementNodes+(1:nOfElementNodes);
    indRC = Te;
    f(indRC) = f(indRC) + fe;
    K(indRC,indRC)=K(indRC,indRC)+ mu*Ke;
end


%NGeo = referenceElementRef.N; NxiGeo = referenceElementRef.Nxi; NetaGeo = referenceElementRef.Neta;
N = referenceElementRef.N; Nxi = referenceElementRef.Nxi; Neta = referenceElementRef.Neta;
IPw = referenceElementRef.IPweights; ngauss = length(IPw);
for i = 1:length(refinedElements)
    iElem = refinedElements(i);
    mu=muElem(iElem);
    equivRef = find(EquivalenciaNumeracio==iElem);
    Te = Tref(equivRef,:); Xe = Xref(Te,:);
    %     Te = T(iElem,:); Xe = X(Te,:);
    figure(1),hold on, plot(Xe(:,1),Xe(:,2),'ko','MarkerSize',8);
    % Jacobian
    J11 = Nxi*Xe(:,1); J12 = Nxi*Xe(:,2); J21 = Neta*Xe(:,1); J22 = Neta*Xe(:,2);
    detJ = J11.*J22-J12.*J21;
    %maybe we should use bsxfun instead of diagonal matrices...
    dvolu = spdiags(referenceElementRef.IPweights.*detJ,0,ngauss,ngauss);
    invJ11 = spdiags(J22./detJ,0,ngauss,ngauss); invJ12 = spdiags(-J12./detJ,0,ngauss,ngauss);
    invJ21 = spdiags(-J21./detJ,0,ngauss,ngauss); invJ22 = spdiags(J11./detJ,0,ngauss,ngauss);
    % xy-derivatives
    Nx = invJ11*Nxi + invJ12*Neta;    Ny = invJ21*Nxi + invJ22*Neta;
    %Computation of r.h.s. source term
    Xg = N*Xe; sourceTerm = -mu*sourceTermFunction(Xg);
    fe = N'*(dvolu*sourceTerm);
    %Elemental matrix
    Ke = Nx'*(dvolu*Nx)+Ny'*(dvolu*Ny);
    %Assembly
    %indRC = (iElem-1)*nOfElementNodes+(1:nOfElementNodes);
    indRC = Tref(equivRef,:)+size(X,1);
    f(indRC) = f(indRC) + fe;
    K(indRC,indRC)=K(indRC,indRC)+ mu*Ke;
end


%Integrals on INTERIOR faces (loop in faces)
%Calculo com si refinat per les dues bandes + projeccio a l'espai std on
%toqui
N2dfaces=referenceElementRef.NfacesIP; N2dxifaces=referenceElementRef.NxifacesIP; N2detafaces=referenceElementRef.NetafacesIP;
N2dGeofaces = referenceElementRef.NGeofacesIP;
N2dxiGeofaces = referenceElementRef.NxiGeofacesIP;
N2detaGeofaces = referenceElementRef.NetaGeofacesIP;
%N1d=referenceElementRef.N1d;
N1dxiGeo = referenceElementRef.N1dxiGeo;
IPw_f = referenceElementRef.IPweights1d; IPz_f=referenceElementRef.IPcoordinates1d; ngf = length(IPw_f);
for iFace=1:nOfNitscheFaces
    
    %% Info de la cara
    infoiFace = NitscheFaces(iFace,:);
    refElem = infoiFace(2); refFaceNum = infoiFace(3);
    stdElem = infoiFace(4); stdFaceNum = infoiFace(5);
    
    %% Funcions aprox
    refN2d = N2dfaces(:,:,refFaceNum);
    stdN2d = N2dfaces(end:-1:1,:,stdFaceNum); %flipping of integration points for 2nd element
    refN2dxi = N2dxifaces(:,:,refFaceNum); stdN2dxi = N2dxifaces(end:-1:1,:,stdFaceNum);
    refN2deta = N2detafaces(:,:,refFaceNum); stdN2deta = N2detafaces(end:-1:1,:,stdFaceNum);
    
    %% Funcions geo
    refN2dGeo = N2dGeofaces(:,:,refFaceNum);
    stdN2dGeo = N2dGeofaces(end:-1:1,:,stdFaceNum); %flipping of integration points for 2nd element
    refN2dxiGeo = N2dxiGeofaces(:,:,refFaceNum);
    stdN2dxiGeo = N2dxiGeofaces(end:-1:1,:,stdFaceNum);
    refN2detaGeo = N2detaGeofaces(:,:,refFaceNum);
    stdN2detaGeo = N2detaGeofaces(end:-1:1,:,stdFaceNum);
    
    %% Std element derivatives
    iElem = stdElem;
    Nxi = stdN2dxi; Neta = stdN2deta;
    NxiGeo = stdN2dxiGeo; NetaGeo = stdN2detaGeo;
    Te = T(iElem,:); Xe = X(Te,:); Testd = Te;
    %  figure(1), hold on, plot(Xe(:,1),Xe(:,2),'r*'), hold off
    J11 = NxiGeo*Xe(:,1); J12 = NxiGeo*Xe(:,2);
    J21 = NetaGeo*Xe(:,1); J22 = NetaGeo*Xe(:,2);
    detJ = J11.*J22-J12.*J21;
    invJ11 = spdiags(J22./detJ,0,ngf,ngf); invJ12 = spdiags(-J12./detJ,0,ngf,ngf);
    invJ21 = spdiags(-J21./detJ,0,ngf,ngf); invJ22 = spdiags(J11./detJ,0,ngf,ngf);
    stdN2dx = invJ11*Nxi + invJ12*Neta;    stdN2dy = invJ21*Nxi + invJ22*Neta;
    %% Ref element derivatives
    iElem = refElem; equivRef = find(EquivalenciaNumeracio==iElem);
    Nxi = refN2dxi; Neta = refN2deta;
    NxiGeo = refN2dxiGeo; NetaGeo = refN2detaGeo;
    Te = T(iElem,:); Xe = X(Te,:);
    Teref = Tref(equivRef,:)+size(X,1);
    %  figure(1), hold on, plot(Xe(:,1),Xe(:,2),'bo'), hold off
    J11 = NxiGeo*Xe(:,1); J12 = NxiGeo*Xe(:,2); J21 = NetaGeo*Xe(:,1); J22 = NetaGeo*Xe(:,2); detJ = J11.*J22-J12.*J21;
    invJ11 = spdiags(J22./detJ,0,ngf,ngf); invJ12 = spdiags(-J12./detJ,0,ngf,ngf);
    invJ21 = spdiags(-J21./detJ,0,ngf,ngf); invJ22 = spdiags(J11./detJ,0,ngf,ngf);
    refN2dx = invJ11*Nxi + invJ12*Neta;    refN2dy = invJ21*Nxi + invJ22*Neta;
    %Elemental matrices
    muL = muElem(refElem); muR = muElem(stdElem);
    %% Integrals on the face (as seen from ref element)
    nodes = faceNodes(refFaceNum,:); Xf = X(T(refElem,nodes),:); % Nodes in the face
    figure(1), hold on, plot(Xf(:,1),Xf(:,2),'gs','MarkerSize',15), hold off
    dxdxi = N1dxiGeo*Xf(:,1); dydxi = N1dxiGeo*Xf(:,2);
    dxdxiNorm = sqrt(dxdxi.^2+dydxi.^2);
    dline = dxdxiNorm.*IPw_f;
    nx = dydxi./dxdxiNorm; ny=-dxdxi./dxdxiNorm; %normal exterior to the left element
    Ddline = spdiags(dline,0,ngf,ngf); Ddline_nx = spdiags(dline.*nx,0,ngf,ngf); Ddline_ny = spdiags(dline.*ny,0,ngf,ngf);
    ref_dline_n_muGradN=muL*(Ddline_nx*refN2dx+Ddline_ny*refN2dy);
    std_dline_n_muGradN=muR*(Ddline_nx*stdN2dx+Ddline_ny*stdN2dy);
    KeLL=-0.5*refN2d'*ref_dline_n_muGradN ...
        -0.5*ref_dline_n_muGradN'*refN2d ...
        +beta*refN2d'*(Ddline*refN2d);
    KeRR=0.5*stdN2d'*std_dline_n_muGradN ...
        +0.5*std_dline_n_muGradN'*stdN2d ...
        +beta*stdN2d'*(Ddline*stdN2d);
    KeLR=-0.5*refN2d'*std_dline_n_muGradN ...
        +0.5*ref_dline_n_muGradN'*stdN2d ...
        -beta*refN2d'*(Ddline*stdN2d);
    %     KeRL=+0.5*rightN2d'*left_dline_n_muGradN ...
    %          -0.5*right_dline_n_muGradN'*leftN2d ...
    %          -beta*rightN2d'*(Ddline*leftN2d);
    
    
    %% Reduccio matrius (right = std)
    P = referenceElementRef.P2d;
    KeRR = P*KeRR*P';
    KeLR = KeLR*P';
    
    
    %% Ensamblatge
    K(Teref,Teref)=K(Teref,Teref)+KeLL;
    K(Teref,Testd)=K(Teref,Testd)+KeLR;
    K(Testd,Teref)=K(Testd,Teref)+KeLR';
    K(Testd,Testd)=K(Testd,Testd)+KeRR;
    
end

% %Integrals on EXTERIOR faces (loop in faces)
% for iFace=1:nOfExteriorFaces
%     infoiFace=infoFaces.extFaces(iFace,:);
%     leftElem=infoiFace(1); leftFaceNum=infoiFace(2);
%     N2d=N2dfaces(:,:,leftFaceNum);
%     N2dxi=N2dxifaces(:,:,leftFaceNum);
%     N2deta=N2detafaces(:,:,leftFaceNum);
%     %left element derivatives
%     iElem=leftElem; Nxi=N2dxi; Neta=N2deta;
%     Te = T(iElem,:); Xe = X(Te,:);
%     J11 = Nxi*Xe(:,1); J12 = Nxi*Xe(:,2); J21 = Neta*Xe(:,1); J22 = Neta*Xe(:,2); detJ = J11.*J22-J12.*J21;
%     invJ11 = spdiags(J22./detJ,0,ngf,ngf); invJ12 = spdiags(-J12./detJ,0,ngf,ngf);
%     invJ21 = spdiags(-J21./detJ,0,ngf,ngf); invJ22 = spdiags(J11./detJ,0,ngf,ngf);
%     N2dx = invJ11*Nxi + invJ12*Neta;    N2dy = invJ21*Nxi + invJ22*Neta;
%     %Elemental matrices
%     muL=muElem(leftElem);
%     nodes = faceNodes(leftFaceNum,:); Xf = Xe(nodes,:); % Nodes in the face
%     dxdxi = N1dxi*Xf(:,1); dydxi = N1dxi*Xf(:,2);
%     dxdxiNorm = sqrt(dxdxi.^2+dydxi.^2);
%     dline = dxdxiNorm.*IPw_f;
%     nx = dydxi./dxdxiNorm; ny=-dxdxi./dxdxiNorm; %normal exterior to the left element
%     Ddline = spdiags(dline,0,ngf,ngf); Ddline_nx = spdiags(dline.*nx,0,ngf,ngf); Ddline_ny = spdiags(dline.*ny,0,ngf,ngf);
%     Ke=-muL*(Ddline_nx*N2dx+Ddline_ny*N2dy)'*N2d;
%     Ke = Ke + Ke' + beta*N2d'*(Ddline*N2d);
%     ind = (leftElem-1)*nOfElementNodes+(1:nOfElementNodes);
%     K(ind,ind)=K(ind,ind)+Ke;
%     Xg = N1d*Xf; uD = funcDirichletCondition(Xg);
%     fe = -(Ddline_nx*N2dx+Ddline_ny*N2dy)'*uD+ beta*N2d'*(dline.*uD);
%     f(ind)=f(ind)+fe;
% end