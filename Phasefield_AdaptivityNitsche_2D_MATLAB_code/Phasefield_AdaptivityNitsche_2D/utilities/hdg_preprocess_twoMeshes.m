function [F, infoFaces] = hdg_preprocess_twoMeshes(X,T,referenceElement,test)
% For every face i:
% infoFaces.intFaces(i,:)=[element1 nface1 element2 nface2 node1] for interior faces
% infoFaces.extFaces(i,:)=[element1 nface1] for exterior faces (NO
% dirichlet)
% infoFaces.dirichletFaces(i,:)=[element1 nface1] for dirichlet faces

[F, infoFaces] = hdg_preprocess(T(:,referenceElement.vertexNodes));
nOfInteriorFaces = size(infoFaces.intFaces,1);
nOfExteriorFaces = size(infoFaces.extFaces,1);

extFaces = []; dirichletFaces = []; topBottomFaces = [];
topFaces = []; bottomFaces = [];
topCircleFaces = []; bottomCircleFaces = [];
for iFace = 1:nOfExteriorFaces
    face = infoFaces.extFaces(iFace,:); %[#element, #face]
    
    EFaces = referenceElement.faceNodes;
    
    nodesFace = T(face(1),EFaces(face(2),:)); %els dels vertexs
    XnodesFace = X(nodesFace,:);
    
    if (test==1) %shear
        if(((abs(XnodesFace(1,2)-0.5)<1.e-6)&&(abs(XnodesFace(2,2)-0.5)<1.e-6)) ...
                || ((abs(XnodesFace(1,2)+0.5)<1.e-6)&&(abs(XnodesFace(2,2)+0.5)<1.e-6))) %face on dirichlet boundary
            dirichletFaces = [dirichletFaces;face];
        else
            extFaces = [extFaces;face];
        end
    elseif(test==2) %multi
        if((abs(XnodesFace(1,2))<1.e-5)&&(abs(XnodesFace(2,2))<1.e-5)) %face on dirichlet boundary
            bottomFaces = [bottomFaces;face];
        elseif((abs(XnodesFace(1,2)-2)<1.e-5)&&(abs(XnodesFace(2,2)-2)<1.e-5))
            topFaces = [topFaces;face];
        else
            dirichletFaces = [dirichletFaces;face];
        end
    elseif(test==3) %L-shaped
        if((abs(XnodesFace(1,2)+250)<1.e-6)&&(abs(XnodesFace(2,2)+250)<1.e-6)) %face on dirichlet boundary
            dirichletFaces = [dirichletFaces;face];
        else
            extFaces = [extFaces;face];
        end
    elseif(test==4) %branching
        if((abs(XnodesFace(1,2)+1)<1.e-5)&&(abs(XnodesFace(2,2)+1)<1.e-5)) %face on dirichlet boundary
            bottomFaces = [bottomFaces;face];
        elseif((abs(XnodesFace(1,2)-1)<1.e-5)&&(abs(XnodesFace(2,2)-1)<1.e-5))
            topFaces = [topFaces;face];
        elseif((abs(XnodesFace(1,1)-1)<1.e-6)&&(abs(XnodesFace(2,1)-1)<1.e-6))
            dirichletFaces = [dirichletFaces;face];
        else
            extFaces = [extFaces;face];
        end
    elseif(test==5) %plate with a hole
        m = size(XnodesFace,1);
        if((length(find(XnodesFace(:,1)==0))==m) || (length(find(XnodesFace(:,1)==65))==m) ...
                || (length(find(XnodesFace(:,2)==0))==m) || (length(find(XnodesFace(:,2)==120))==m))
            extFaces = [extFaces;face];
        elseif((length(find(XnodesFace(:,2)>80))==m))
            topCircleFaces = [topCircleFaces;face];
        elseif(length(find(XnodesFace(:,2)<30))==m)
            bottomCircleFaces = [bottomCircleFaces;face];
        else
            extFaces = [extFaces;face];
        end
    end
end

topBottomFaces = [topFaces;bottomFaces];

if (test == 5)
    dirichletFaces = [bottomCircleFaces;topCircleFaces];
    infoFaces.bottomCircleFaces = bottomCircleFaces;
    infoFaces.topCircleFaces = topCircleFaces;
end


infoFaces.extFaces = extFaces;
infoFaces.dirichletFaces = dirichletFaces;
infoFaces.topBottomFaces = topBottomFaces;

% Rearrangement F
nOfExteriorFaces = size(extFaces,1); nOfDirichletFaces = size(dirichletFaces,1);
for iFace = 1:nOfExteriorFaces
    infoFace = extFaces(iFace,:);
    F(infoFace(1),infoFace(2)) = iFace + nOfInteriorFaces;
end

for iFace = 1:size(topBottomFaces,1)
    infoFace = topBottomFaces(iFace,:);
    F(infoFace(1),infoFace(2)) = iFace + nOfInteriorFaces+nOfExteriorFaces;
end

for iFace = 1:nOfDirichletFaces
    infoFace = dirichletFaces(iFace,:);
    F(infoFace(1),infoFace(2)) = iFace + nOfInteriorFaces+nOfExteriorFaces+size(topBottomFaces,1);
end

