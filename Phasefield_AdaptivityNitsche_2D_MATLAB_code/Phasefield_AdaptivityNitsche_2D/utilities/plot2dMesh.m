function plot2dMesh(X,T,refElemVertexNodes)

Tvertexes=T(:,refElemVertexNodes);

[nOfElements,nOfElementVertexes]=size(Tvertexes);
if nOfElementVertexes==4
    aux=[1:4,1];
else
    aux=[1:3,1];
end

clf, hold on
for i=1:nOfElements
    Xe=X(Tvertexes(i,:),:);
    plot(Xe(aux,1),Xe(aux,2),'k-')
end
plot(X(:,1),X(:,2),'ko')
hold off