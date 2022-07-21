function plotMeshTwoMeshes(NitscheFaces,referenceElementStd,referenceElementRef,refinedElements,X,T,Xref,Tref,u,nDegRef)

% Check input
if nargin == 4
    nDegRef = 20;
end

nOfVertexes=length(referenceElementStd.vertexNodes);

% Plotting element (equal spaced points)
if nOfVertexes==3 %triangle
    nodes = [];
    h = 1/nDegRef;
    for j = 0:nDegRef
        i = (0:nDegRef-j)';
        aux = j*ones(size(i));
        nodes = [nodes; [i aux]*h];
    end
    nodes = 2*nodes - 1;
else
    xs = -1:1/nDegRef:1;
    [nodesx,nodesy] = meshgrid(xs);
    nodes = [nodesx(:),nodesy(:)];
end

elemTriRef = delaunayn(nodes);
nDeg = referenceElementStd.degree;
nOfElements = size(T,1);
standardElements = setdiff(1:nOfElements,refinedElements);
nOfRefEls = length(refinedElements); nOfStdEls = length(standardElements);

% Standard elements
Tstd = T(standardElements,:);
hold on, kk = plotMesh(X+u(1:size(X,1),:),Tstd,1);

% Refined elements
elemTriRef = delaunayn([referenceElementRef.IPcoordinates;referenceElementRef.NodesCoord]);
TNodesCoord=referenceElementRef.TNodesCoord;
% elemTriRef = delaunayn(referenceElementRef.IPcoordinates);
for i = 1:nOfRefEls
    Te = Tref(i,:);
    Xe = Xref(Te,:);
    
    ue = u(Te+size(X,1),:);
    hold on, kk = plotMesh(Xe+ue,TNodesCoord,1);
end


hold off
axis equal

