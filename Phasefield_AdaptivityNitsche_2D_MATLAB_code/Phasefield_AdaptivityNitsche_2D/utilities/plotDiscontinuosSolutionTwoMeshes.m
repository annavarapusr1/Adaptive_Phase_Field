function plotDiscontinuosSolutionTwoMeshes(NitscheFaces,referenceElementStd,referenceElementRef,refinedElements,X,T,Xref,Tref,u,nDegRef)

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
N = evaluateNodalBasisQua(nodes,referenceElementStd.NodesCoord,nDeg);
for i = 1:nOfStdEls
    iElem = standardElements(i);
    Te = T(iElem,:);
    Xplot = N*X(Te,:);
    ue = u(Te,:);
    uplot = N*ue;
    hold on
    patch('Faces',elemTriRef,'Vertices',Xplot,'FaceVertexCData',uplot,...
        'FaceColor','interp','EdgeAlpha',0);
end

% Refined elements
Nref = referenceElementRef.N;
elemTriRef = delaunayn([referenceElementRef.IPcoordinates;referenceElementRef.NodesCoord]);
% elemTriRef = delaunayn(referenceElementRef.IPcoordinates);
for i = 1:nOfRefEls
    Te = Tref(i,:);
    Xe = Xref(Te,:);
    Xplot = [Nref*Xe;Xe];
    ue = u(Te+size(X,1),:);
    uplot = [Nref*ue;ue];
    hold on
    patch('Faces',elemTriRef,'Vertices',Xplot,'FaceVertexCData',uplot,...
        'FaceColor','interp','EdgeAlpha',0);
    plot(Xe(referenceElementRef.vertexNodes,1),Xe(referenceElementRef.vertexNodes,2),'w-','LineWidth',0.05);
end

% Nitsche faces
faceNodesStd = referenceElementStd.faceNodes;
for i  = 1:size(NitscheFaces,1)
    el = NitscheFaces(i,4);
    f = NitscheFaces(i,5);
    Xe = X(T(el,:),:); face = faceNodesStd(f,:);
    Xf = Xe(face,:);
    plot(Xf(:,1),Xf(:,2),'w-','LineWidth',1);
end

hold off
axis equal
%shading interp
colormap jet
