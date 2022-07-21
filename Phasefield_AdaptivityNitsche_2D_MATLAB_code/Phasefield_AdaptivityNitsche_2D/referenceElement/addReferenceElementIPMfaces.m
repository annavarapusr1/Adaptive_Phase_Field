function theReferenceElement = addReferenceElementIPMfaces(theReferenceElement)
%Adds to theReferenceElement the basis functions and their derivatives at the integration points on the faces

nOfFaces=size(theReferenceElement.faceNodes,1);
nOfNodes=size(theReferenceElement.N,2);
nDeg=theReferenceElement.degree;
xiIP1d=theReferenceElement.IPcoordinates1d;
faceNodes=theReferenceElement.faceNodes;
nodesCoord=theReferenceElement.NodesCoord;
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
if nOfFaces==3
    for i=1:nOfFaces
        [N,Nxi,Neta]=evaluateNodalBasisTri(facesIPcoords(:,:,i),nodesCoord,nDeg);
        N(abs(N)<1.e-6)=0;
        NfacesIP(:,:,i)=N; NxifacesIP(:,:,i)=Nxi; NetafacesIP(:,:,i)=Neta;
    end
else %Quadrilaterals are not tested...
    for i=1:nOfFaces
        [N,Nxi,Neta]=evaluateNodalBasisQua(facesIPcoords(:,:,i),nDeg);
        N(abs(N)<1.e-6)=0;
        NfacesIP(:,:,i)=N; NxifacesIP(:,:,i)=Nxi; NetafacesIP(:,:,i)=Neta;
    end   
end

theReferenceElement.NfacesIP=NfacesIP;
theReferenceElement.NxifacesIP=NxifacesIP;
theReferenceElement.NetafacesIP=NetafacesIP;
theReferenceElement.facesIPcoords=facesIPcoords;




