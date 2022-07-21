%Identification of nodes at corners, faces, etc 
%for a given list of nodal coordinates in a quadrilateral [-1,1]^2

%Example with equally spaced nodes
degree=6;
x=linspace(-1,1,degree+1);
[x,y]=meshgrid(x,x);
coord2d=[reshape(x,(degree+1)^2,1),reshape(y,(degree+1)^2,1)];

%Identification of vertexes
Xi=[-1,-1]; [d,v1]=min((coord2d(:,1)-Xi(1)).^2+(coord2d(:,2)-Xi(2)).^2);
Xi=[ 1,-1]; [d,v2]=min((coord2d(:,1)-Xi(1)).^2+(coord2d(:,2)-Xi(2)).^2);
Xi=[ 1, 1]; [d,v3]=min((coord2d(:,1)-Xi(1)).^2+(coord2d(:,2)-Xi(2)).^2);
Xi=[-1, 1]; [d,v4]=min((coord2d(:,1)-Xi(1)).^2+(coord2d(:,2)-Xi(2)).^2);
vertexNodes=[v1,v2,v3,v4];

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

% %Plot for checking
% aux=reshape(faceNodes,1,(degree+1)*4); 
% figure(1), plot(coord2d(innerNodes,1),coord2d(innerNodes,2),'bo',coord2d(aux,1),coord2d(aux,2),'r*')





