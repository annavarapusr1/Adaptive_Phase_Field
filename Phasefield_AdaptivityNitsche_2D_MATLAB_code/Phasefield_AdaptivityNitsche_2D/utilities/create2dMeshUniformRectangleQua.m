function [X,T]=create2dMeshUniformRectangleQua(refNodalCoord,intervals,nxy)
%[X,T]=create2dMeshUniformRectangleQua(referenceElement.NodalCoordinates,intervals,nxy)
% Input:
%  refNodalCoord: coordinates of the nodes in the referenceElement [-1,1]^3
%  intervals = [a,b,c,d] for domain [a,b]x[c,d]
%  nxy=[nx,ny] for number of divisions in each direction

a=intervals(1); b=intervals(2); hx=(b-a)/nxy(1);
c=intervals(3); d=intervals(4); hy=(d-c)/nxy(2);

x=linspace(a,b,nxy(1)+1);
y=linspace(c,d,nxy(2)+1);

X1=(refNodalCoord+1)/2;
x1=X1(:,1)*hx; y1=X1(:,2)*hy;  %coordinates in [0,hx]x[0,hy]

X=[]; T=[]; aux=[1:length(x1)]; kk=0;
for i=1:nxy(1)
    for j=1:nxy(2)
        X=[X; [x1+x(i) y1+y(j)]];
        T=[T ; aux+kk];
        kk=kk+length(x1);
    end
end

%Identification of repeated coordinates
inode=1;
while(inode<size(X,1))
    Xi=X(inode,:);
    ii=find( abs(Xi(1)-X(inode+1:end,1))+ abs(Xi(2)-X(inode+1:end,2))<1.e-10);
    for k=1:length(ii)
        i=ii(k)+inode;
        X(i,:)=[];
        T(T==i)=inode; T(T>i)=T(T>i)-1;
        ii=ii-1;
    end
    inode=inode+1;
end

