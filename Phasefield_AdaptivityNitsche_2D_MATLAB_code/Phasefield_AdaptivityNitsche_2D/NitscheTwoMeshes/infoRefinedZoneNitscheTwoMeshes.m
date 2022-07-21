function [Xref,Tref,NitscheFaces,nodesCCDStd,nodesCCDRef,nodesCCDStdCorrected,nonStdNodes,refinedElements,RefinedDirichletFaces,nodesTopBottomStd,nodesTopBottomRef] = infoRefinedZoneNitscheTwoMeshes(X,T,F,infoFaces,elementsToRefine,referenceElementRef,referenceElementStd,test)

nOfElements = size(T,1);
nOfExteriorFaces = size(infoFaces.extFaces,1);
nOfElementFaces = size(F,2);
nOfInteriorFaces = size(infoFaces.intFaces,1);
nOfTopBottomFaces = size(infoFaces.topBottomFaces,1);

%% Xref, Tref
standardElements = setdiff(1:nOfElements,elementsToRefine);
stdNodes = T(standardElements,:); stdNodes = unique(stdNodes(:));
refinedElements = elementsToRefine';
Tref = []; Xref = [];
gluingFaces = [];

for iref=refinedElements
    XrefNew = referenceElementRef.NGeoNodesCoord*X(T(iref,:),:); %transformation
    
    [Xref,Tref] = glueMeshes(Xref,XrefNew,Tref,1:size(XrefNew,1));
    
    for f = 1:nOfElementFaces
        iFace = F(iref,f);
        aux = find(gluingFaces==iFace);
        if isempty(aux)
            gluingFaces = [gluingFaces,iFace];
        else
            gluingFaces(aux) = [];
        end
    end
end
%refinedExteriorFaces = gluingFaces(gluingFaces>nOfInteriorFaces & gluingFaces<(nOfInteriorFaces+nOfExteriorFaces+1));
refinedDirichletFaces = gluingFaces(gluingFaces > (nOfInteriorFaces+nOfExteriorFaces+nOfTopBottomFaces));
gluingFaces = gluingFaces(gluingFaces<(nOfInteriorFaces+1));
refinedTopBottomFaces = gluingFaces(gluingFaces > (nOfInteriorFaces+nOfExteriorFaces) & gluingFaces < (nOfInteriorFaces+nOfExteriorFaces+nOfTopBottomFaces+1));

NitscheFaces = zeros(length(gluingFaces),5);
faceNodesStd = referenceElementStd.faceNodes;
for j = 1:length(gluingFaces)
    f = gluingFaces(j);
    [elements,faces] = find(F == f);
    %Xe1 = X(T(elements(1),:),:);
    %Xf = Xe1(faceNodesStd(faces(1),:),:);
    if (ismember(elements(1),refinedElements)) %primer element ref
        NitscheFaces(j,:) = [f,elements(1),faces(1),elements(2),faces(2)];
    else
        NitscheFaces(j,:) = [f,elements(2),faces(2),elements(1),faces(1)];
    end
end

StandardDirichletFaces = infoFaces.dirichletFaces;
RefinedDirichletFaces = zeros(length(refinedDirichletFaces),2);
faceNodesRef = referenceElementRef.faceNodes;
for j = 1:length(refinedDirichletFaces)
    f = refinedDirichletFaces(j);
    [element,face] = find(F == f);
    equivRef = find(refinedElements==element);
    Xe = Xref(Tref(equivRef,:),:);
    %Xf = Xe(faceNodesRef(face,:),:);
    RefinedDirichletFaces(j,:) = [equivRef,face];
    %Remove refined face from StandardExteriorFaces
    aux = find(StandardDirichletFaces(:,1)==element);
    aux2 = find(StandardDirichletFaces(aux,2)==face);
    StandardDirichletFaces(aux(aux2),:) = [];
end

nNodes1dstd = size(referenceElementStd.faceNodes,2);
nodesCCDStd = zeros(1,size(StandardDirichletFaces,1)*nNodes1dstd);
for f = 1:size(StandardDirichletFaces,1)
    infoFace = StandardDirichletFaces(f,:);
    Te = T(infoFace(1),:); Xe = X(Te,:);
    %Xf = Xe(faceNodesStd(infoFace(2),:),:);
    nodesFace = Te(faceNodesStd(infoFace(2),:));
    nodesCCDStd((f-1)*nNodes1dstd+(1:nNodes1dstd)) = nodesFace;
end
nodesCCDStd = unique(nodesCCDStd);

nNodes1dref = size(referenceElementRef.faceNodes,2);
nodesCCDRef = zeros(1,size(RefinedDirichletFaces,1)*nNodes1dref);
for f = 1:size(RefinedDirichletFaces,1)
    infoFace = RefinedDirichletFaces(f,:);
    Te = Tref(infoFace(1),:); 
    nodesFace = Te(faceNodesRef(infoFace(2),:));
    nodesCCDRef((f-1)*nNodes1dref+(1:nNodes1dref)) = nodesFace;
end
nodesCCDRef = unique(nodesCCDRef);

%% Rows/columns to be removed from K,f
nonStdNodes = setdiff(1:size(X,1),stdNodes);

nodesCCDStdCorrected = nodesCCDStd;
for j = length(nonStdNodes):-1:1
    row = nonStdNodes(j);
    fn = find(nodesCCDStdCorrected>row);
    nodesCCDStdCorrected(fn) = nodesCCDStdCorrected(fn)-1;
end


