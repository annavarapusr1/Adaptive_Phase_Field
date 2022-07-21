function [X,T]=glueMeshesRef(X1,X2,T1,T2)

X=[X1;X2]; T=[T1;T2+size(X1,1)];

%Identification of repeated coordinates
inode=size(X1,1)+1;
while(inode<size(X,1))
    Xi=X(inode,:);
    ii=find( abs(Xi(1)-X(inode+1:end,1))+ abs(Xi(2)-X(inode+1:end,2))<1.e-9);
    
       
    for k=1:length(ii)
        i=ii(k)+inode;
        X(i,:)=[];
        T(T==i)=inode; T(T>i)=T(T>i)-1;
        ii=ii-1;
    end
    
    inode=inode+1;
end
