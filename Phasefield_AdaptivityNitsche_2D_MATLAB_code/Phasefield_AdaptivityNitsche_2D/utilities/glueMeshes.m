function [X,T]=glueMeshes3d(X1,X2,T1,T2)

% X=[X1;X2]; T=[T1;T2+size(X1,1)];

% %Identification of repeated coordinates
% inode=1;
% while(inode<size(X,1))
%     Xi=X(inode,:);
%     ii=find((abs(Xi(1)-X(inode+1:end,1))+abs(Xi(2)-X(inode+1:end,2))+abs(Xi(3)-X(inode+1:end,3)))<1.e-10);
%     for k=1:length(ii)
%         i=ii(k)+inode;
%         X(i,:)=[];
%         T(T==i)=inode; T(T>i)=T(T>i)-1;
%         ii=ii-1;
%     end
%     inode=inode+1;
% end

% if (~isempty(X1))
%     [i_node_X1,i_node_X2]=find(squeeze(sum(abs(bsxfun(@minus,X1,permute(X2,[3 2 1]))),2))<1e-8);
%     T2new = T2+size(X1,1);
%     
%     for i = 1:length(i_node_X2)
%         inode = i_node_X2(i);
%         T2new(T2==inode) = i_node_X1(i);
%     end
%     
%     X = [X1;X2(setdiff(1:size(X2,1),i_node_X2),:)];
%     T = [T1;T2new];
%     new_num_T = sparse(unique(T),1,1:numel(unique(T)));
%     T = full(new_num_T(T));
% else
%     X = X2; T = T2; 
% end

if(isempty(X1))
   X1 = X2; X2 = [];
   T1 = T2; T2 = [];
end

X = [X1;zeros(size(X2,1),size(X2,2))]; %T = [T1;T2+size(X1,1)];
T2new = T2+size(X1,1);
cont = 0;
for inode = 1:size(X2,1)
    Xi = X2(inode,:);
    X1minusXi =  sum(abs(X1-Xi),2);
    node_duplicated = find((X1minusXi)<1.e-8);

    
%     %% only for shear test!!
%     if((Xi(:,1)<-1.e-4) && abs(Xi(:,2))<1.e-4 && ~isempty(node_duplicated))
%        %saber si estem a dalt o a baix
%        [els,cols] = find(T1 == node_duplicated(1));
%        Xaux = X1(T1(els(1),:),:);
%        
%        if numel(node_duplicated)==1
%            if ~(any(X2(:,2)>1.e-5)&&any(Xaux(:,2)>1.e-5)||any(X2(:,2)<-1.e-5)&&any(Xaux(:,2)<-1.e-5))
%                node_duplicated = [];
%            end
%        else
%            if (any(X2(:,2)>1.e-5)&&any(Xaux(:,2)>1.e-5)||any(X2(:,2)<-1.e-5)&&any(Xaux(:,2)<-1.e-5))
%                node_duplicated = node_duplicated(1);
%            else
%                node_duplicated = node_duplicated(2);
%            end   
%        end
%     elseif abs(Xi(:,1))<1.e-4 && abs(Xi(:,2))<1.e-4        
%         a = 3;
%     end
    
    if node_duplicated
        
        T2new(T2==inode) = node_duplicated;
        T2new(T2>inode)=T2new(T2>inode)-1;
       
    else
        cont = cont+1;
        X(cont+size(X1,1),:) = Xi;
        T2new(T2==inode) = cont+size(X1,1);
    end
end
X = X(1:(size(X1,1)+cont),:);
T = [T1;T2new];
