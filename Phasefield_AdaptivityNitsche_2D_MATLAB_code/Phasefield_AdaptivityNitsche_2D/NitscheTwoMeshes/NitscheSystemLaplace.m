function [K,f] = IPMsystemLaplace(muElem,X,T,referenceElement,infoFaces,beta,sourceTermFunction,funcDirichletCondition)

nOfElements = size(T,1);
nOfInteriorFaces = size(infoFaces.intFaces,1);
nOfExteriorFaces = size(infoFaces.extFaces,1);
nOfFacesElem = size(referenceElement.faceNodes,1);
nOfFaceNodes = size(referenceElement.NodesCoord1d,1);
nOfElementNodes = size(referenceElement.NodesCoord,1);
faceNodes = referenceElement.faceNodes;

%TEST (u=xy^2, v=x+y)
%   utest=zeros(nOfElements*nOfElementNodes,1); vtest=utest;
%   for iElem=1:nOfElements
%       ind=(iElem-1)*nOfElementNodes+[1:nOfElementNodes];
%       utest(ind)=X(T(iElem,:),1).*(X(T(iElem,:),2).^2);
%       vtest(ind)=X(T(iElem,:),1)+X(T(iElem,:),2);
%   end

nDOF = nOfElements*nOfElementNodes;
f = zeros(nDOF,1);
K = spalloc(nDOF,nDOF,nOfFaceNodes*nDOF);
%n = nOfElements*(nOfElementNodes)^2; indK = 1:nOfElementNodes; ind_i=zeros(1,n); ind_j=ind_i; coef_K=ind_i;

% Volume integrals (loop in elements)
N = referenceElement.N; Nxi = referenceElement.Nxi; Neta = referenceElement.Neta;
IPw = referenceElement.IPweights; ngauss = length(IPw);
for iElem = 1:nOfElements
    mu=muElem(iElem);
    Te = T(iElem,:); Xe = X(Te,:);
    % Jacobian
    J11 = Nxi*Xe(:,1); J12 = Nxi*Xe(:,2); J21 = Neta*Xe(:,1); J22 = Neta*Xe(:,2);
    detJ = J11.*J22-J12.*J21;
    %maybe we should use bsxfun instead of diagonal matrices...
    dvolu = spdiags(referenceElement.IPweights.*detJ,0,ngauss,ngauss);
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
    indRC = (iElem-1)*nOfElementNodes+(1:nOfElementNodes);
    f(indRC) = f(indRC) + fe;
    K(indRC,indRC)=K(indRC,indRC)+ mu*Ke;
%     for irow = 1:nOfElementNodes
%         ind_i(indK)=indRC(irow); ind_j(indK)=indRC; coef_K(indK) = mu*Ke(irow,:); indK = indK+nOfElementNodes;
%     end
%TEST
%     ue=utest(indRC); ve=vtest(indRC);
%     test1=ve'*Ke*ue; test1int=diag(dvolu)'*(Xg(:,2).*(Xg(:,2)+2*Xg(:,1)));
%     [test1,test1int,test1int-test1]
end
%K = sparse(ind_i,ind_j,coef_K);

%Integrals on INTERIOR faces (loop in faces)
N2dfaces=referenceElement.NfacesIP; N2dxifaces=referenceElement.NxifacesIP; N2detafaces=referenceElement.NetafacesIP;
N1d=referenceElement.N1d; N1dxi=referenceElement.N1dxi;
IPw_f = referenceElement.IPweights1d; IPz_f=referenceElement.IPcoordinates1d; ngf = length(IPw_f);
for iFace=1:nOfInteriorFaces
    infoiFace=infoFaces.intFaces(iFace,:);
    leftElem=infoiFace(1); leftFaceNum=infoiFace(2); rightElem=infoiFace(3); rightFaceNum=infoiFace(4); 
    leftN2d=N2dfaces(:,:,leftFaceNum); rightN2d=N2dfaces(end:-1:1,:,rightFaceNum); %flipping of integration points for 2nd element
    leftN2dxi=N2dxifaces(:,:,leftFaceNum); rightN2dxi=N2dxifaces(end:-1:1,:,rightFaceNum); 
    leftN2deta=N2detafaces(:,:,leftFaceNum); rightN2deta=N2detafaces(end:-1:1,:,rightFaceNum); 
    %right element derivatives
    iElem=rightElem; Nxi=rightN2dxi; Neta=rightN2deta;
    Te = T(iElem,:); Xe = X(Te,:);
    %  figure(1), hold on, plot(Xe(:,1),Xe(:,2),'r*'), hold off
    J11 = Nxi*Xe(:,1); J12 = Nxi*Xe(:,2); J21 = Neta*Xe(:,1); J22 = Neta*Xe(:,2); detJ = J11.*J22-J12.*J21;
    invJ11 = spdiags(J22./detJ,0,ngf,ngf); invJ12 = spdiags(-J12./detJ,0,ngf,ngf);
    invJ21 = spdiags(-J21./detJ,0,ngf,ngf); invJ22 = spdiags(J11./detJ,0,ngf,ngf);
    rightN2dx = invJ11*Nxi + invJ12*Neta;    rightN2dy = invJ21*Nxi + invJ22*Neta;  
    %left element derivatives
    iElem=leftElem; Nxi=leftN2dxi; Neta=leftN2deta;
    Te = T(iElem,:); Xe = X(Te,:);
    %  figure(1), hold on, plot(Xe(:,1),Xe(:,2),'bo'), hold off
    J11 = Nxi*Xe(:,1); J12 = Nxi*Xe(:,2); J21 = Neta*Xe(:,1); J22 = Neta*Xe(:,2); detJ = J11.*J22-J12.*J21;
    invJ11 = spdiags(J22./detJ,0,ngf,ngf); invJ12 = spdiags(-J12./detJ,0,ngf,ngf);
    invJ21 = spdiags(-J21./detJ,0,ngf,ngf); invJ22 = spdiags(J11./detJ,0,ngf,ngf);
    leftN2dx = invJ11*Nxi + invJ12*Neta;    leftN2dy = invJ21*Nxi + invJ22*Neta;  
    %Elemental matrices
    muL=muElem(leftElem); muR=muElem(rightElem);
    %Integrals on the face (as seen from left element)
    nodes = faceNodes(leftFaceNum,:); Xf = X(T(leftElem,nodes),:); % Nodes in the face
    %  figure(1), hold on, plot(Xf(:,1),Xf(:,2),'gs'), hold off
    dxdxi = N1dxi*Xf(:,1); dydxi = N1dxi*Xf(:,2);
    dxdxiNorm = sqrt(dxdxi.^2+dydxi.^2);
    dline = dxdxiNorm.*IPw_f;
    nx = dydxi./dxdxiNorm; ny=-dxdxi./dxdxiNorm; %normal exterior to the left element
    Ddline = spdiags(dline,0,ngf,ngf); Ddline_nx = spdiags(dline.*nx,0,ngf,ngf); Ddline_ny = spdiags(dline.*ny,0,ngf,ngf);
    left_dline_n_muGradN=muL*(Ddline_nx*leftN2dx+Ddline_ny*leftN2dy);
    right_dline_n_muGradN=muR*(Ddline_nx*rightN2dx+Ddline_ny*rightN2dy);
    KeLL=-0.5*leftN2d'*left_dline_n_muGradN ...
         -0.5*left_dline_n_muGradN'*leftN2d ...
         +beta*leftN2d'*(Ddline*leftN2d);
    KeRR=0.5*rightN2d'*right_dline_n_muGradN ...
        +0.5*right_dline_n_muGradN'*rightN2d ...
        +beta*rightN2d'*(Ddline*rightN2d);
    KeLR=-0.5*leftN2d'*right_dline_n_muGradN ...
         +0.5*left_dline_n_muGradN'*rightN2d ...
         -beta*leftN2d'*(Ddline*rightN2d);
%     KeRL=+0.5*rightN2d'*left_dline_n_muGradN ...
%          -0.5*right_dline_n_muGradN'*leftN2d ...
%          -beta*rightN2d'*(Ddline*leftN2d);     
    Lind = (leftElem-1)*nOfElementNodes+(1:nOfElementNodes);
    Rind = (rightElem-1)*nOfElementNodes+(1:nOfElementNodes);
    K(Lind,Lind)=K(Lind,Lind)+KeLL;
    K(Lind,Rind)=K(Lind,Rind)+KeLR;
    K(Rind,Lind)=K(Rind,Lind)+KeLR';
    K(Rind,Rind)=K(Rind,Rind)+KeRR;
    %TEST integration points in face (need for flipping :-))
    %XfL=leftN2d*X(T(leftElem,:),:);    XfR=rightN2d*X(T(rightElem,:),:);
    
%     %TEST
%     uL=utest(Lind); uR=utest(Rind);
%     vL=vtest(Lind); vR=vtest(Rind); 
%     Xg=N1d*Xf; xf=Xg(:,1);yf=Xg(:,2);
% 
%     disp(sprintf('\nCara %d, elements %d i %d',iFace,leftElem,rightElem))
%     testa=dline'*((nx+ny).*(xf.*(yf.^2)));
%     Ketest=right_dline_n_muGradN'*rightN2d;
%     testb=vR'*Ketest*uR; 
%     [testa,testb],error=testa-testb, if(abs(error)>1.e-6) disp('AINxxxxxxxxxxxxxxxxxxx!'), end
% 
%     testa=dline'*((xf+yf).*(nx.*(yf.^2)+ny.*(2*xf.*yf)));
%     Ketest=leftN2d'*left_dline_n_muGradN;
%     testb=vL'*Ketest*uL; 
%     [testa,testb],error=testa-testb, if(abs(error)>1.e-6) disp('AINxxxxxxxxxxxxxxxxxxx!'), end
% 
%     testa=dline'*((xf+yf).*(nx.*(yf.^2)+ny.*(2*xf.*yf)));
%     Ketest=rightN2d'*right_dline_n_muGradN;
%     testb=vR'*Ketest*uR; 
%     [testa,testb],error=testa-testb, if(abs(error)>1.e-6) disp('AINxxxxxxxxxxxxxxxxxxx!'), end

end

%Integrals on EXTERIOR faces (loop in faces)
for iFace=1:nOfExteriorFaces
    infoiFace=infoFaces.extFaces(iFace,:);
    leftElem=infoiFace(1); leftFaceNum=infoiFace(2); 
    N2d=N2dfaces(:,:,leftFaceNum); 
    N2dxi=N2dxifaces(:,:,leftFaceNum); 
    N2deta=N2detafaces(:,:,leftFaceNum); 
    %left element derivatives
    iElem=leftElem; Nxi=N2dxi; Neta=N2deta;
    Te = T(iElem,:); Xe = X(Te,:);
    J11 = Nxi*Xe(:,1); J12 = Nxi*Xe(:,2); J21 = Neta*Xe(:,1); J22 = Neta*Xe(:,2); detJ = J11.*J22-J12.*J21;
    invJ11 = spdiags(J22./detJ,0,ngf,ngf); invJ12 = spdiags(-J12./detJ,0,ngf,ngf);
    invJ21 = spdiags(-J21./detJ,0,ngf,ngf); invJ22 = spdiags(J11./detJ,0,ngf,ngf);
    N2dx = invJ11*Nxi + invJ12*Neta;    N2dy = invJ21*Nxi + invJ22*Neta;  
    %Elemental matrices
    muL=muElem(leftElem); 
    nodes = faceNodes(leftFaceNum,:); Xf = Xe(nodes,:); % Nodes in the face
    dxdxi = N1dxi*Xf(:,1); dydxi = N1dxi*Xf(:,2);
    dxdxiNorm = sqrt(dxdxi.^2+dydxi.^2);
    dline = dxdxiNorm.*IPw_f;
    nx = dydxi./dxdxiNorm; ny=-dxdxi./dxdxiNorm; %normal exterior to the left element
    Ddline = spdiags(dline,0,ngf,ngf); Ddline_nx = spdiags(dline.*nx,0,ngf,ngf); Ddline_ny = spdiags(dline.*ny,0,ngf,ngf);
    Ke=-muL*(Ddline_nx*N2dx+Ddline_ny*N2dy)'*N2d; 
    Ke = Ke + Ke' + beta*N2d'*(Ddline*N2d);
    ind = (leftElem-1)*nOfElementNodes+(1:nOfElementNodes);
    K(ind,ind)=K(ind,ind)+Ke;
    Xg = N1d*Xf; uD = funcDirichletCondition(Xg);
    fe = -(Ddline_nx*N2dx+Ddline_ny*N2dy)'*uD+ beta*N2d'*(dline.*uD);
    f(ind)=f(ind)+fe;
end

%TEST
%   utest=zeros(nOfElements*nOfElementNodes,1); 
%   for iElem=1:nOfElements
%       ind=(iElem-1)*nOfElementNodes+[1:nOfElementNodes];
%       utest(ind)=X(T(iElem,:),1);
%   end
% residu=K*utest-f;
% residuElem=reshape(residu,nOfElementNodes,nOfElements);
% residuElemMaxAbs=sum(abs(residuElem),1);
% wrongElements=find(residuElemMaxAbs>1.e-6);
% XwrongElements=[];
% for iElem=wrongElements, XwrongElements=[XwrongElements;sum(X(T(iElem,1:3),:))/3]; end
% if length(XwrongElements)>0, figure(1), hold on, plot(XwrongElements(:,1),XwrongElements(:,2),'ro','markerSize',10), hold off, end
