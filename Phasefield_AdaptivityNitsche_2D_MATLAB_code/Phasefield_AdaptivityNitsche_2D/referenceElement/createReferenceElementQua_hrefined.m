function referenceElementRef = createReferenceElementQua_hrefined(degree,refinementFactor)

referenceElementGeo = createReferenceElement(0,degree);
nOfFaces = 4;

% Nodes coordinates refined referenceElement
href = 2/refinementFactor;
nOfElements = refinementFactor^2;
xs = -1:href:1; 
QUA = zeros(nOfElements,4);
[X,Y] = meshgrid(xs); j = 1;
nodes= [X(:),Y(:)];
for i = 1:length(nodes)
    x = nodes(i,1); y = nodes(i,2); 
    if(x>(-1+1.e-6) && y>(-1+1.e-6))
        QUA(j,:) = [i,i-1,(i-1)-length(xs),i-length(xs)];
        j = j+1;
    end 
end

NodesCoordGeo = referenceElementGeo.NodesCoord; nOfNodesGeo = size(NodesCoordGeo,1);
T = zeros(size(QUA,1),nOfNodesGeo);
NodesCoord = zeros(size(QUA,1)*nOfNodesGeo,2);
for i = 1:size(QUA,1)
    Te = QUA(i,:);
    nodesTe = nodes(Te,:);
    Xe = (NodesCoordGeo+1)/2*href + nodesTe(3,:);
    NodesCoord((i-1)*nOfNodesGeo+(1:nOfNodesGeo),:) = Xe;
    T(i,:) = (i-1)*size(Xe,1)+(1:size(Xe,1));
end
[NodesCoord,~,IC] = uniquetol(NodesCoord,'ByRows',true);
T = IC(T); %matriu de connectivitats amb numeracio actual de NodesCoord
if(size(T,2)==1)
    T = T';
end
%figure(2),clf,kk=plotMesh(NodesCoord,T,1,'plotNodes');

nOfNodes = size(NodesCoord,1);

nodes1d = xs'; 
NodesCoord1dGeo = referenceElementGeo.NodesCoord1d; nOfNodes1dGeo = size(NodesCoord1dGeo,1);
T1d = zeros(refinementFactor,nOfNodes1dGeo);
NodesCoord1d = zeros(refinementFactor*nOfNodes1dGeo,1);
for i = 1:refinementFactor
    Tf = i:i+1; nodesTf = nodes1d(Tf);
    Xf = (NodesCoord1dGeo+1)/2*href + nodesTf(1);
    NodesCoord1d((i-1)*nOfNodes1dGeo+(1:nOfNodes1dGeo),:) = Xf;
    T1d(i,:) = (i-1)*size(Xf,1)+(1:size(Xf,1));
end
[NodesCoord1d,~,IC] = uniquetol(NodesCoord1d);
T1d = IC(T1d);
if(size(T1d,2)==1) T1d = T1d'; end

referenceElementRef.NodesCoord = NodesCoord;
referenceElementRef.TNodesCoord = T;
referenceElementRef.NodesCoord1d = NodesCoord1d;
referenceElementRef.NodesCoordGeo = referenceElementGeo.NodesCoord;
referenceElementRef.NodesCoord1dGeo = referenceElementGeo.NodesCoord1d;
referenceElementRef.refinementFactor = refinementFactor;

% nodes vertices quadrilateral
vertexNodes = zeros(4,1);
vertexNodes(1) = find(abs(NodesCoord(:,1)+1)<1.e-6 & abs(NodesCoord(:,2)+1)<1.e-6);
vertexNodes(2) = find(abs(NodesCoord(:,1)-1)<1.e-6 & abs(NodesCoord(:,2)+1)<1.e-6);
vertexNodes(3) = find(abs(NodesCoord(:,1)-1)<1.e-6 & abs(NodesCoord(:,2)-1)<1.e-6);
vertexNodes(4) = find(abs(NodesCoord(:,1)+1)<1.e-6 & abs(NodesCoord(:,2)-1)<1.e-6);
referenceElementRef.vertexNodes = vertexNodes;

% faceNodesGeo
referenceElementRef.faceNodesGeo = referenceElementGeo.faceNodes;

% faceNodes
% face1: (-1,-1) to (1,-1)
nodesInFace = find(abs(NodesCoord(:,2)+1)<1.e-6);
coordNodesFace = NodesCoord(nodesInFace,:);
[~,indexs] = sort(coordNodesFace(:,1),'ascend');
nodesInFace = nodesInFace(indexs);
faceNodes(1,:) = nodesInFace;
% face2: (1,-1) to (1,1)
nodesInFace = find(abs(NodesCoord(:,1)-1)<1.e-6);
coordNodesFace = NodesCoord(nodesInFace,:);
[~,indexs] = sort(coordNodesFace(:,2),'ascend');
nodesInFace = nodesInFace(indexs);
faceNodes(2,:) = nodesInFace;
% face3: (1,1) to (-1,1)
nodesInFace = find(abs(NodesCoord(:,2)-1)<1.e-6);
coordNodesFace = NodesCoord(nodesInFace,:);
[~,indexs] = sort(coordNodesFace(:,1),'descend');
nodesInFace = nodesInFace(indexs);
faceNodes(3,:) = nodesInFace;
% face4: (-1,1) to (-1,-1)
nodesInFace = find(abs(NodesCoord(:,1)+1)<1.e-6);
coordNodesFace = NodesCoord(nodesInFace,:);
[~,indexs] = sort(coordNodesFace(:,2),'descend');
nodesInFace = nodesInFace(indexs);
faceNodes(4,:) = nodesInFace;

referenceElementRef.faceNodes = faceNodes;

% innerNodes
innerNodes = setdiff(1:size(NodesCoord,1),faceNodes);
referenceElementRef.innerNodes = innerNodes;

nIP1d = length(referenceElementGeo.IPweights1d); nofnodes1d = length(NodesCoord1d);
IPweights1d = zeros(refinementFactor*nIP1d,1); 
N1d = zeros(refinementFactor*nIP1d,nofnodes1d); N1dxi = N1d;
IPcoordinates1d = zeros(refinementFactor*nIP1d,1);
NGeo1d = referenceElementGeo.N1d;
Nx1dGeo = referenceElementGeo.N1dxi;

for i = 1:refinementFactor
    Tf = T1d(i,:);
    Xf = NodesCoord1d(Tf);
    Xg = NGeo1d*Xf;
    IPcoordinates1d((i-1)*nIP1d+(1:nIP1d)) = Xg;
    
    dxdxi = Nx1dGeo*Xf;
    dxdxiNorm = sqrt(dxdxi.^2);
    IPweights1d((i-1)*nIP1d+(1:nIP1d)) = dxdxiNorm.*(referenceElementGeo.IPweights1d);
    
    Nx1d = diag(1./dxdxiNorm)*Nx1dGeo;   
    N1d((i-1)*nIP1d+(1:nIP1d),Tf) = NGeo1d;
    N1dxi((i-1)*nIP1d+(1:nIP1d),Tf) = Nx1d;
end

referenceElementRef.N1d = sparse(N1d);
referenceElementRef.N1dxi = sparse(N1dxi);
referenceElementRef.faceNodes1d = 1:nofnodes1d;
referenceElementRef.IPcoordinates1d = IPcoordinates1d;
referenceElementRef.IPweights1d = IPweights1d;

NGeo = referenceElementGeo.N;
NxiGeo = referenceElementGeo.Nxi; NetaGeo = referenceElementGeo.Neta;

% IPCoordinates = []; IPweights = [];
nIP = length(referenceElementGeo.IPweights);
IPCoordinates = zeros(nOfElements*nIP,2); 
IPweights = zeros(nOfElements*nIP,1);
N = zeros(size(T,1)*nIP,size(NodesCoord,1));
Nxi = N; Neta = N;

% NfacesIP = N2d evaluades als IP de cadascuna de les cares
NfacesIP = zeros(size(IPcoordinates1d,1),nOfNodes,nOfFaces);
NxifacesIP = NfacesIP; NetafacesIP = NfacesIP;

faceNodesGeo = referenceElementGeo.faceNodes;
NodesCoordGeo = referenceElementGeo.NodesCoord;
IPcoordinates1dGeo = referenceElementGeo.IPcoordinates1d;
nIP1dgeo = length(IPcoordinates1dGeo);
facesIPcoordsGeo = zeros(size(IPcoordinates1dGeo,1),2,nOfFaces);
a = (1-IPcoordinates1dGeo)/2; b = (1+IPcoordinates1dGeo)/2;
for i=1:nOfFaces
    A = NodesCoordGeo(faceNodesGeo(i,1),:); B = NodesCoordGeo(faceNodesGeo(i,end),:);
    facesIPcoordsGeo(:,:,i) = a*A+b*B;
end

cont = zeros(4,1);
for i = 1:size(T,1)
    Te = T(i,:);
    Xe = NodesCoord(Te,:);
    Xg = NGeo*Xe;
    ind = (i-1)*nIP+(1:nIP);
    
    IPCoordinates(ind,:) = Xg;
    
    J11 = NxiGeo*Xe(:,1); J12 = NxiGeo*Xe(:,2);
    J21 = NetaGeo*Xe(:,1); J22 = NetaGeo*Xe(:,2);
    detJ = J11.*J22-J12.*J21;
    IPweights(ind) = referenceElementGeo.IPweights.*detJ;
    invJ11 = spdiags(J22./detJ,0,nIP,nIP);
    invJ12 = spdiags(-J12./detJ,0,nIP,nIP);
    invJ21 = spdiags(-J21./detJ,0,nIP,nIP);
    invJ22 = spdiags(J11./detJ,0,nIP,nIP);
    % xy-derivatives
    Nx = invJ11*NxiGeo + invJ12*NetaGeo;
    Ny = invJ21*NxiGeo + invJ22*NetaGeo;

    N(ind,Te) = NGeo;
    Nxi(ind,Te) = Nx;
    Neta(ind,Te) = Ny;
    
    for iface = 1:nOfFaces
        nodes_iface = Te(faceNodesGeo(iface,:));
        if (sum(ismember(nodes_iface,faceNodes(iface,:)))==length(nodes_iface)) %la cara del subelement es exterior
            [N_iface,Nxi_iface,Neta_iface] = evaluateNodalBasisQua(facesIPcoordsGeo(:,:,iface),NodesCoordGeo,degree);
            J11 = Nxi_iface*Xe(:,1); J12 = Nxi_iface*Xe(:,2);
            J21 = Neta_iface*Xe(:,1); J22 = Neta_iface*Xe(:,2);
            detJ = J11.*J22-J12.*J21;
            invJ11 = spdiags(J22./detJ,0,nIP1dgeo,nIP1dgeo);
            invJ12 = spdiags(-J12./detJ,0,nIP1dgeo,nIP1dgeo);
            invJ21 = spdiags(-J21./detJ,0,nIP1dgeo,nIP1dgeo);
            invJ22 = spdiags(J11./detJ,0,nIP1dgeo,nIP1dgeo);
            Nx_iface = invJ11*Nxi_iface+invJ12*Neta_iface;
            Ny_iface = invJ21*Nxi_iface+invJ22*Neta_iface;
            if(iface < 3)
                ind = (cont(iface))*nIP1dgeo+(1:nIP1dgeo);
            else
                ind = ((refinementFactor-1)*nIP1dgeo)+(1:nIP1dgeo)-(cont(iface)*nIP1dgeo);
            end
            NfacesIP(ind,Te,iface) = N_iface;
            NxifacesIP(ind,Te,iface) = Nx_iface;
            NetafacesIP(ind,Te,iface) = Ny_iface;
            cont(iface) = cont(iface)+1;
        end
    end
end

%hold on, plot(IPCoordinates(:,1),IPCoordinates(:,2),'k*','MarkerSize',5);

cellNfacesIP = cell(nOfFaces,1);
cellNxifacesIP = cell(nOfFaces,1); cellNetafacesIP = cell(nOfFaces,1);
for i = 1:nOfFaces
    cellNfacesIP{i} = sparse(NfacesIP(:,:,i));
    cellNxifacesIP{i} = sparse(NxifacesIP(:,:,i));
    cellNetafacesIP{i} = sparse(NetafacesIP(:,:,i));
end

referenceElementRef.IPcoordinates = IPCoordinates;
referenceElementRef.IPweights = IPweights;
referenceElementRef.N = sparse(N);
referenceElementRef.Nxi = sparse(Nxi);
referenceElementRef.Neta = sparse(Neta);
referenceElementRef.NfacesIP = cellNfacesIP;
referenceElementRef.NxifacesIP = cellNxifacesIP;
referenceElementRef.NetafacesIP = cellNetafacesIP;

[N,Nxi,Neta]=evaluateNodalBasisQua(referenceElementRef.IPcoordinates,referenceElementRef.NodesCoordGeo,degree);
[N1,Nxi1]=evaluateNodalBasis1D(referenceElementRef.IPcoordinates1d,referenceElementRef.NodesCoord1dGeo,degree);

referenceElementRef.NGeo=sparse(N);
referenceElementRef.NxiGeo=sparse(Nxi);
referenceElementRef.NetaGeo=sparse(Neta);
referenceElementRef.N1dGeo=sparse(N1);
referenceElementRef.N1dxiGeo=sparse(Nxi1);
referenceElementRef.degree = degree;
referenceElementRef.degreeGeo = degree;

referenceElementRef.IPcoordinates1dGeo = referenceElementGeo.IPcoordinates1d;

% Matrix P st N1dGeo=N1d*P (equiv. N1dgeo'=P*N1d')
[Ptransposed,~]=evaluateNodalBasis1D(referenceElementRef.NodesCoord1d,referenceElementGeo.NodesCoord1d,degree);
referenceElementRef.P = Ptransposed';

%IPM faces (N2dGeo = N2d*P2d)
facesIPcoords=zeros(size(IPcoordinates1d,1),2,nOfFaces);
a=(1-IPcoordinates1d)/2; b=(1+IPcoordinates1d)/2;
for i=1:nOfFaces
    A=NodesCoord(faceNodes(i,1),:); B=NodesCoord(faceNodes(i,end),:);
    facesIPcoords(:,:,i)=a*A+b*B;
end

[P2dtransposed,~]=evaluateNodalBasisQua(NodesCoord,NodesCoordGeo,degree);
referenceElementRef.P2d = P2dtransposed';

cellNGeofacesIP = cell(nOfFaces,1);
cellNxiGeofacesIP = cell(nOfFaces,1); cellNetaGeofacesIP = cell(nOfFaces,1);
for i=1:nOfFaces
    [N,Nxi,Neta]=evaluateNodalBasisQua(facesIPcoords(:,:,i),referenceElementRef.NodesCoordGeo,degree);
    N(abs(N)<1.e-6)=0;
    cellNGeofacesIP{i} = sparse(N);
    cellNxiGeofacesIP{i} = sparse(Nxi);
    cellNetaGeofacesIP{i} = sparse(Neta);
end

%referenceElementRef.facesIPcoords=facesIPcoords;
referenceElementRef.NGeofacesIP = cellNGeofacesIP;
referenceElementRef.NxiGeofacesIP = cellNxiGeofacesIP;
referenceElementRef.NetaGeofacesIP = cellNetaGeofacesIP;

[NGeoNodesCoord,~,~] = evaluateNodalBasisQua(NodesCoord,referenceElementGeo.NodesCoord,degree);
referenceElementRef.NGeoNodesCoord = NGeoNodesCoord;
