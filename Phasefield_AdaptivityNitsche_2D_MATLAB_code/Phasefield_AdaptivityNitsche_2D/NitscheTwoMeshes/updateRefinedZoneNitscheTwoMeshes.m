function [Xref,Tref,NitscheFaces,nodesCCDStd,nodesCCDRef,nodesCCDStdCorrected,nonStdNodes,refinedElements,d,ux,uy,H_previous,nodesTopBottomStd,nodesTopBottomRef] = updateRefinedZoneNitscheTwoMeshes(X,T,Xref,Tref,F,infoFaces,NitscheFaces,refinedElementsOLD,elementsToRefine,referenceElementRef,referenceElementStd,nodesCCDStd,nodesCCDRef,d,ux,uy,H_previous,test,nodesTopBottomStd,nodesTopBottomRef)

nOfElements = size(T,1);
nOfExteriorFaces = size(infoFaces.extFaces,1);
nOfElementFaces = size(F,2);
nOfInteriorFaces = size(infoFaces.intFaces,1);
nOfTopBottomFaces = size(infoFaces.topBottomFaces,1);

% Xref, Tref
gluingFaces = NitscheFaces(:,1);
newGluingFaces = [];
newRefDirichletFaces = [];
newRefExteriorFaces = [];
newRefTopBottomFaces = [];
for iref = elementsToRefine'
    XrefNew = referenceElementRef.NGeoNodesCoord*X(T(iref,:),:); %transformation
    
    [Xref,Tref] = glueMeshes(Xref,XrefNew,Tref,1:size(XrefNew,1));
    
    for f = 1:nOfElementFaces
        iFace = F(iref,f);
        aux = find(gluingFaces==iFace);
        if isempty(aux)
            if(iFace < (nOfInteriorFaces+1))
                auxNew = find(newGluingFaces==iFace);
                if isempty(auxNew)
                    newGluingFaces = [newGluingFaces,iFace];
                else
                    newGluingFaces(auxNew) = [];
                end
            elseif(iFace>(nOfInteriorFaces+nOfExteriorFaces+nOfTopBottomFaces))
                newRefDirichletFaces = [newRefDirichletFaces,iFace];
            elseif(iFace < (nOfInteriorFaces+nOfExteriorFaces+1))
                newRefExteriorFaces = [newRefExteriorFaces,iFace];
            else
                newRefTopBottomFaces = [newRefTopBottomFaces,iFace];
            end
        else
            gluingFaces(aux) = [];
            NitscheFaces(aux,:) = [];
        end
    end
end
newGluingFaces = unique(newGluingFaces);
refinedElements = [refinedElementsOLD,elementsToRefine'];
standardElements = setdiff(1:nOfElements,refinedElements);

% NitscheFaces
faceNodesStd = referenceElementStd.faceNodes;
for j = 1:length(newGluingFaces)
    f = newGluingFaces(j);
    [elements,faces] = find(F == f);
    Xe1 = X(T(elements(1),:),:);
    Xf = Xe1(faceNodesStd(faces(1),:),:);
    if (ismember(elements(1),refinedElements)) %first element ref
        info = [f,elements(1),faces(1),elements(2),faces(2)];
    else
        info = [f,elements(2),faces(2),elements(1),faces(1)];
    end
    NitscheFaces = [NitscheFaces; info];
end

% nodesCCD
stdNodes = T(standardElements,:); stdNodes = unique(stdNodes(:));
nonStdNodes = setdiff(1:size(X,1),stdNodes);
nodesCCDStd = setdiff(nodesCCDStd,nonStdNodes); 
faceNodesRef = referenceElementRef.faceNodes;
for j = 1:length(newRefDirichletFaces)
    f = newRefDirichletFaces(j);
    [element,face] = find(F == f);
    equivRef = find(refinedElements==element);
    Te = Tref(equivRef,:);
    Xe = Xref(Te,:);
    Xf = Xe(faceNodesRef(face,:),:);
    nodesFace = Te(faceNodesRef(face,:));
    nodesCCDRef = [nodesCCDRef,nodesFace];
end
nodesCCDRef = unique(nodesCCDRef);

nodesCCDStdCorrected = nodesCCDStd;
for j = length(nonStdNodes):-1:1
    row = nonStdNodes(j);
    fn = find(nodesCCDStdCorrected>row);
    nodesCCDStdCorrected(fn) = nodesCCDStdCorrected(fn)-1;
end

% nodesTopBottom
if (test==2)
    nodesTopBottomStd = setdiff(nodesTopBottomStd,nonStdNodes); %remove non standard
    faceNodesRef = referenceElementRef.faceNodes;
    for j = 1:length(newRefTopBottomFaces)
        f = newRefTopBottomFaces(j);
        [element,face] = find(F == f);
        equivRef = find(refinedElements==element);
        Te = Tref(equivRef,:);
        Xe = Xref(Te,:);
        %Xf = Xe(faceNodesRef(face,:),:);
        nodesFace = Te(faceNodesRef(face,:));
        nodesTopBottomRef = [nodesTopBottomRef,nodesFace];
        %figure(1),hold on, plot(Xf(:,1),Xf(:,2),'k-','LineWidth',3);
    end
    nodesTopBottomRef = unique(nodesTopBottomRef);
end

%% Projection of the solution: d, u, H, Hprevious
dstd = d(1:size(X,1)); dref = d(size(X,1)+1:end);
ux_std = ux(1:size(X,1)); ux_ref = ux(size(X,1)+1:end);
uy_std = uy(1:size(X,1)); uy_ref = uy(size(X,1)+1:end);

Nstd_nodesRef = referenceElementRef.NGeoNodesCoord;
nOfNewDofs = size(Xref,1) - length(dref);

d = [dstd;dref;zeros(nOfNewDofs,1)];
ux = [ux_std;ux_ref;zeros(nOfNewDofs,1)];
uy = [uy_std;uy_ref;zeros(nOfNewDofs,1)];

%toRemoveFromStd = [];
for iElem = elementsToRefine'
    Te = T(iElem,:); Xe = X(Te,:);
    de = dstd(Te); ux_e = ux_std(Te); uy_e = uy_std(Te);
    
    d_projection = Nstd_nodesRef*de;
    ux_projection = Nstd_nodesRef*ux_e;
    uy_projection = Nstd_nodesRef*uy_e;
    
    %toRemoveFromStd = [toRemoveFromStd,setdiff(Te,stdNodes)];
    
    Tref_e = Tref(find(refinedElements==iElem),:);
    Xref_e = Xref(Tref_e,:);
    ind = size(X,1) + Tref_e;
    d(ind) = d_projection;
    ux(ind) = ux_projection;
    uy(ind) = uy_projection;
end

% Projection H_previous
nIPstd = length(referenceElementStd.IPweights);
nStdElsOLD = length(setdiff(1:nOfElements,refinedElementsOLD));
Hprev_std = H_previous(1:(nStdElsOLD*nIPstd)); Hprev_ref = H_previous((nStdElsOLD*nIPstd+1):end);

Hprev_refNEW = [];

standardElementsOLD = setdiff(1:nOfElements,refinedElementsOLD);
[Nstd_ipref,~,~] = evaluateNodalBasisQua(referenceElementRef.IPcoordinates,referenceElementStd.NodesCoord,referenceElementStd.degree);

toremove = [];

for iElem = elementsToRefine'
    i  = find(standardElementsOLD == iElem);
    ind = (i-1)*nIPstd + (1:nIPstd);
    Hprev_e = Hprev_std(ind);
    
    toremove = [toremove,ind];
    
    N = referenceElementStd.N;
    NxiGeo = referenceElementStd.Nxi;
    NetaGeo = referenceElementStd.Neta;
    IPwstd = referenceElementStd.IPweights;
    
    Xe = X(T(iElem,:),:);
    J11 = NxiGeo*Xe(:,1); J12 = NxiGeo*Xe(:,2);
    J21 = NetaGeo*Xe(:,1); J22 = NetaGeo*Xe(:,2);
    detJ = J11.*J22-J12.*J21;
    dvolu = diag(IPwstd.*detJ);
    M = N'*(dvolu*N);
    b_Hprev = N'*(dvolu*Hprev_e);
    
    Hprev_nodesstd = M\b_Hprev;
    Hprev_projection = Nstd_ipref*Hprev_nodesstd;
    Hprev_refNEW = [Hprev_refNEW; Hprev_projection];
end

Hprev_std(toremove) = [];
H_previous = [Hprev_std;Hprev_ref;Hprev_refNEW];
