function referenceElementRef = addReferenceElementIPMfaces(referenceElementRef)
%Adds to theReferenceElement the basis functions and their derivatives at the integration points on the faces

nOfFaces=size(referenceElementRef.faceNodes,1);
nOfNodes=size(referenceElementRef.N,2);
nOfNodesStd = size(referenceElementRef.NGeo,2);
nDeg=referenceElementRef.degree;
xiIP1d=referenceElementRef.IPcoordinates1d;
faceNodes=referenceElementRef.faceNodes;
faceNodesGeo=referenceElementRef.faceNodesGeo;

nodesCoord=referenceElementRef.NodesCoord;
nodesCoordGeo = referenceElementRef.NodesCoordGeo;
a=(1-xiIP1d)/2; b=(1+xiIP1d)/2;
%Coordinates of the integration points on the faces
facesIPcoords=zeros(size(xiIP1d,1),2,nOfFaces);
for i=1:nOfFaces
    A=nodesCoord(faceNodes(i,1),:); B=nodesCoord(faceNodes(i,end),:);
    facesIPcoords(:,:,i)=a*A+b*B;
end
%Basis functions and derivatives at these points
NfacesIP=zeros(size(xiIP1d,1),nOfNodes,nOfFaces);
NxifacesIP=zeros(size(xiIP1d,1),nOfNodes,nOfFaces);
NetafacesIP=zeros(size(xiIP1d,1),nOfNodes,nOfFaces);

%%%%%%
T = referenceElementRef.TNodesCoord;
xiIP1dgeo = referenceElementRef.IPcoordinates1dGeo;
facesIPcoordsGeo=zeros(size(xiIP1dgeo,1),2,nOfFaces);
NGeofacesIP=zeros(size(xiIP1dgeo,1),nOfNodesStd,nOfNodesStd);
NxiGeofacesIP=zeros(size(xiIP1dgeo,1),nOfNodesStd,nOfNodesStd);
NetaGeofacesIP=zeros(size(xiIP1dgeo,1),nOfNodesStd,nOfNodesStd);
a=(1-xiIP1dgeo)/2; b=(1+xiIP1dgeo)/2;
for i=1:nOfFaces
    A=nodesCoordGeo(faceNodesGeo(i,1),:); B=nodesCoordGeo(faceNodesGeo(i,end),:);
    facesIPcoordsGeo(:,:,i)=a*A+b*B;
end
for i=1:nOfFaces
    [N,Nxi,Neta]=evaluateNodalBasisQua(facesIPcoordsGeo(:,:,i),nodesCoordGeo,nDeg);
    N(abs(N)<1.e-6)=0;
    NGeofacesIP(:,:,i)=N; NxiGeofacesIP(:,:,i)=Nxi; NetaGeofacesIP(:,:,i)=Neta;
end
% nofip1d = length(xiIP1dgeo);
% for m = 1:size(T,1)
%     Xe = nodesCoord(T(m,:),:);
%     figure(2),hold on, plot(Xe(:,1),Xe(:,2),'bo');
%     
%     for i = 1:nOfFaces
%         Ni = NGeofacesIP(:,:,i); Nxii = NxiGeofacesIP(:,:,i); Netai = NetaGeofacesIP(:,:,i);
%         
%         NfacesIP((m-1)*nofip1d+(1:nofip1d),Te,i)=Ni; 
%         NxifacesIP((m-1)*nofip1d+(1:nofip1d),Te,i)=Nx_i; 
%         NetafacesIP((m-1)*nofip1d+(1:nofip1d),Te,i)=Ny_i;
%     end
% end

%%%%%

NGeofacesIP=zeros(size(xiIP1d,1),nOfNodes,nOfFaces);
NxiGeofacesIP=zeros(size(xiIP1d,1),nOfNodes,nOfFaces);
NetaGeofacesIP=zeros(size(xiIP1d,1),nOfNodes,nOfFaces);

for i=1:nOfFaces
    [N,Nxi,Neta]=evaluateNodalBasisQua(facesIPcoords(:,:,i),nodesCoordGeo,nDeg);
    N(abs(N)<1.e-6)=0;
    NGeofacesIP(:,:,i)=N; NxiGeofacesIP(:,:,i)=Nxi; NetaGeofacesIP(:,:,i)=Neta;
end

referenceElementRef.facesIPcoords=facesIPcoords;

referenceElementRef.NfacesIP=NfacesIP;
referenceElementRef.NxifacesIP=NxifacesIP;
referenceElementRef.NetafacesIP=NetafacesIP;

referenceElementRef.NGeofacesIP = NGeofacesIP;
referenceElementRef.NxiGeofacesIP = NxiGeofacesIP;
referenceElementRef.NetaGeofacesIP = NetaGeofacesIP;


