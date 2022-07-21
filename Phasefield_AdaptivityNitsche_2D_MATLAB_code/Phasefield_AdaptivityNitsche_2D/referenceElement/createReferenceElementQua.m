function theReferenceElement=createReferenceElementQua(degree)

% theReferenceElement=createReferenceElementTri(degree)
% Output:
%  theReferenceElement: struct containing
%     .IPcoordinates: coordinates of the integration points for 2D elemens
%     .IPweights: weights of the integration points for 2D elements
%     .N: shape functions at the IP
%     .Nxi,.Neta: derivatives of the shape functions at the IP
%     .IPcoordinates1d: coordinates of the integration points for 1D boundary elemens
%     .IPweights1d: weights of the integration points for 1D boundary elements
%     .N1d: 1D shape functions at the IP
%     .N1dxi: derivatives of the 1D shape functions at the IP
%     .faceNodes: matrix [nOfFaces nOfNodesPerFace] with the edge nodes numbering
%     .innerNodes: vector [1 nOfInnerNodes] with the inner nodes numbering
%     .faceNodes1d: vector [1 nOfNodesPerElement] with the 1D nodes numbering
%     .NodesCoord: spatial coordinates of the element nodes
%     .NodesCoord1d: spatial coordinates of the 1D element nodes

switch degree
    case 1 %Q1
        faceNodes = [1 2; 2 3; 3 4; 4 1];
        innerNodes = [];
        faceNodes1d = 1:2;
        coord2d = [-1 -1; 1 -1; 1 1; -1 1];
        coord1d = [-1; 1];
        vertexNodes=1:4;
%     case 2 %Q2
%         faceNodes = [1 5 2; 2 6 3; 3 7 4; 4 8 1];
%         innerNodes = 9;
%         faceNodes1d = 1:3;
%         coord2d = [-1 -1; 1 -1; 1 1; -1 1; 0 -1; 1 0; 0 1; -1 0; 0 0];
%         coord1d = [-1; 0; 1];
%         vertexNodes=1:4;
    case {2,3,4,5,6,7,8,9,10}
        faceNodes1d = 1:degree+1;
        coord1d= feketeNodes1D(degree,[1:degree+1]);
        coord2d=ez4u_2Dqua_NodalCoordinatesReferenceElement(degree);
        %Identificatin of vertexes
        vX=[-1 -1; 1 -1; 1 1; -1 1];
        vertexNodes=[];
        for i=1:4
            vertexNodes=[vertexNodes find(abs(coord2d(:,1)-vX(i,1))+abs(coord2d(:,2)-vX(i,2))<1.e-10 )];
        end
        %Identification of faces nodes (faceNodes)
        nodes=find( abs(coord2d(:,2)+1)<1.e-10 ); [kk,perm]=sort(coord2d(nodes,1)); nodes=nodes(perm); %{eta=-1}
        faceNodes=nodes';
        nodes=find( abs(coord2d(:,1)-1)<1.e-10 ); [kk,perm]=sort(coord2d(nodes,2)); nodes=nodes(perm); %{xi=1}
        faceNodes=[faceNodes; nodes'];
        nodes=find( abs(coord2d(:,2)-1)<1.e-10 ); [kk,perm]=sort(-coord2d(nodes,1)); nodes=nodes(perm); %{eta=1}
        faceNodes=[faceNodes; nodes'];
        nodes=find( abs(coord2d(:,1)+1)<1.e-10 ); [kk,perm]=sort(-coord2d(nodes,2)); nodes=nodes(perm); %{xi=-1}
        faceNodes=[faceNodes; nodes'];
        %InteriorNodes
        innerNodes=setdiff([1:(degree+1)^2],reshape(faceNodes,1,(degree+1)*4));
    otherwise
        error('??? Reference element not implemented yet')
end

nOfGaussPoints1D = 2*degree+1; %number of integration points 1D
[gp1d,gw1d]=gaussLegendre(nOfGaussPoints1D);
%2D numerical quadrature (tensor product of 1D quadrature)
nOfGaussPoints2D=nOfGaussPoints1D^2;
[a1,a2]=meshgrid(gp1d,gp1d);
gp2d=[reshape(a1,nOfGaussPoints2D,1),reshape(a2,nOfGaussPoints2D,1)];
[a1,a2]=meshgrid(gw1d,gw1d);
gw2d=reshape(a1.*a2,nOfGaussPoints2D,1);

[N,Nxi,Neta]=evaluateNodalBasisQua(gp2d,coord2d,degree);
[N1,Nxi1]=evaluateNodalBasis1D(gp1d,coord1d,degree);

%Creating reference element structure
theReferenceElement = struct('IPcoordinates',gp2d,...
    'IPweights',gw2d,'N',N,'Nxi',Nxi,'Neta',Neta,...
    'IPcoordinates1d',gp1d,'IPweights1d',gw1d',...
    'N1d',N1,'N1dxi',Nxi1,'faceNodes',faceNodes,'vertexNodes',vertexNodes,...
    'innerNodes',innerNodes,'faceNodes1d',faceNodes1d,...
    'NodesCoord',coord2d,'NodesCoord1d',coord1d,'degree',degree);

