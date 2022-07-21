clear, close all
setpath

%% PARAMETERS
% Resolution parameters
degree = 1;
test = 1; %1: shear, 3: Lshaped, 5: hole
tol = 1.e-2;
h = 1/24; refinementFactor = 10; refinementValue = 0.2;
dt = 1.e-4; tf = 0.02;
increments = dt:dt:tf;
eta = 0;
Gc = 2.7e-3; l = 0.015;
E = 210; nu = 0.3;
saveDisplacements = []; %1.e-3:1.e-3:tf;

resultsPath = 'Resultats2D/ShearTest/p1_m24r10/';

% Nitsche parameters
alphaLE = 100; alphaD = 100;
betaLE = E*alphaLE/((h/refinementFactor));
betaD = Gc*l*alphaD/((h/refinementFactor));

% Phase-field parameters
lambda = E*nu/((1+nu)*(1-2*nu));
mu = E/(2*(1+nu));

%% REFERENCE ELEMENTS
referenceElementStd = createReferenceElement(0,degree);
referenceElementRef = createReferenceElementQua_hrefined(degree,refinementFactor);

%% MESH
[X,T] = create2dMeshUniformRectangleQua(referenceElementStd.NodesCoord,[-0.5 0.5 -0.5 0.5],[1/h, 1/h]);
%figure(1),clf,aux=plotMesh(X,T,1,'plotNodes');

%% Crack in the domain -> modification of the mesh
nodesCrack = find((X(:,1)<-1.e-8)&(abs(X(:,2))<1.e-6));
inode = 1;
for node = nodesCrack'
    [rows,cols] = find(T == node); %rows = elements containing the node
    for element = rows'
        elNodes = T(element,:);
        if (find(X(elNodes,2) > 1.e-4)) %new node (duplicated) for the element above the crack
            aux = find(T(element,:) == node);
            T(element,aux) = size(X,1)+inode;
        end
    end
    X = [X; X(node,:)+[0,1.e-6]];
    ++inode;
end

nOfElements = size(T,1);
nOfElementNodes = size(T,2);

% Material parameters: constant in every element
E_elems = E*ones(nOfElements,1); nu_elems = nu*ones(nOfElements,1);
Gc_elems = Gc*ones(nOfElements,1); l_elems = l*ones(nOfElements,1);

%% INITIAL REFINED ZONE
% DG preprocess (identification of faces)
[F,infoFaces] = hdg_preprocess_twoMeshes(X,T(:,:),referenceElementStd,test);
nOfFaces = max(max(F));
nOfExteriorFaces = size(infoFaces.extFaces,1);
nOfDirichletFaces = size(infoFaces.dirichletFaces,1);
nOfElementFaces = size(F,2);
nOfInteriorFaces = nOfFaces - (nOfExteriorFaces-nOfDirichletFaces);

% Elements to refine
vN = referenceElementStd.vertexNodes;
X1 = X(T(:,vN(1)),:); X2 = X(T(:,vN(2)),:); X3 = X(T(:,vN(3)),:); X4 = X(T(:,vN(4)),:);
Xmig = (X1+X2+X3+X4)/4;
elementsToRefine = find(abs(Xmig(:,1))<h & abs(Xmig(:,2))<h);

[Xref,Tref,NitscheFaces,nodesCCDStd,nodesCCDRef,nodesCCDStdCorrected,nonStdNodes,refinedElements] = infoRefinedZoneNitscheTwoMeshes(X,T,F,infoFaces,elementsToRefine,referenceElementRef,referenceElementStd);

nonStdNodes_xy = [nonStdNodes,nonStdNodes+size(X,1)+size(Xref,1)];
stdNodes = setdiff(1:size(X,1),nonStdNodes);
nDOFStd = length(stdNodes); nDOFRef = size(Xref,1);
nDOF = (nDOFStd+nDOFRef);
nodesCCD = [nodesCCDStdCorrected,nodesCCDStdCorrected+nDOF];
notCCD = setdiff(1:(2*nDOF),nodesCCD); %actual degrees of freedom (not boundary nodes)

d = zeros(size(X,1)+size(Xref,1),1);

nOfRefEls = length(refinedElements); nOfStdEls = nOfElements-nOfRefEls;
nIPref = length(referenceElementRef.IPweights); nIPstd = length(referenceElementStd.IPweights);
H_previous = zeros(nIPref*nOfRefEls + nIPstd*nOfStdEls,1);

loads = zeros(length(increments),1); 

%% LOOP IN TIME STEPS
for ind = 1:length(increments)
    
    k = increments(ind);
    fprintf('\nCurrent applied displacement: %f \n',k)

    s = 0; stoppingCriterion = 0;
    
    while (stoppingCriterion == 0)
        s = s + 1;
        
        %% EQUILIBRIUM EQUATION
        uCCD = dirichletCondition(X,Xref,nodesCCDStd,nodesCCDRef,k,test);
        
        [K,f] = NitscheTwoMeshesLinearElasticityDamage_PF(X,T,Xref,Tref,NitscheFaces,refinedElements,referenceElementStd,referenceElementRef,E_elems,nu_elems,eta,betaLE,d);
        
        K(nonStdNodes_xy,:) = [];
        K(:,nonStdNodes_xy) = [];
        f(nonStdNodes_xy) = [];
        
        Kdn = K(nodesCCD,notCCD); Kdd = K(nodesCCD,nodesCCD);
        
        f = f(notCCD)-K(notCCD,nodesCCD)*uCCD;
        K = K(notCCD,notCCD);
        sol = K\f;
        u = zeros(2*nDOF,1);
        u(notCCD) = sol; u(nodesCCD) = uCCD;
        ux = u(1:nDOF); uy = u(nDOF+1:end);
        
        % reconstruccio ux,uy
        ux_std = ux(1:nDOFStd); ux_ref = ux(nDOFStd+1:end);
        uy_std = uy(1:nDOFStd); uy_ref = uy(nDOFStd+1:end);
        ux = zeros(size(X,1) + size(Xref,1),1); uy = ux;
        ux(setdiff(1:size(X,1),nonStdNodes)) = ux_std; ux(size(X,1)+1:end) = ux_ref;
        uy(setdiff(1:size(X,1),nonStdNodes)) = uy_std; uy(size(X,1)+1:end) = uy_ref;
        
        %% HISTORY FIELD H (at integration points)
        H = computeH(X,T,Xref,Tref,refinedElements,referenceElementStd,referenceElementRef,ux,uy,H_previous,lambda,mu);
        
        %% DAMAGE EQUATION
        
        %% Computation
        [K,f] = NitscheTwoMeshesDamage_PF(X,T,Xref,Tref,NitscheFaces,refinedElements,referenceElementStd,referenceElementRef,Gc_elems,l_elems,betaD,H);
        K(nonStdNodes,:) = [];
        K(:,nonStdNodes) = [];
        f(nonStdNodes) = [];
        d = K\f;
        
        % reconstruction d
        dStd = d(1:nDOFStd); dRef = d(nDOFStd+1:end);
        d = zeros(size(X,1) + size(Xref,1),1);
        d(setdiff(1:size(X,1),nonStdNodes)) = dStd;
        d(size(X,1)+1:end) = dRef;
        ploT = [0:0.01:tf] ;
        if(ismembertol(k,ploT))
          figure(2),clf,plotDiscontinuosSolutionTwoMeshes(NitscheFaces,referenceElementStd,referenceElementRef,refinedElements,X,T,Xref,Tref,d,10)
        end 
        
        %% STOPPING CRITERION: segons damage field
        if (s == 1)
            d_previous = d;
        else
            errorDamage = computeEuclideanNormRelative(d,d_previous);
            if (errorDamage < tol)
                stoppingCriterion = 1;
            else
                d_previous = d;
            end
            fprintf('      Damage residual at iteration %d, |R|',s-1)
            fprintf(' = %.4e \n', errorDamage)
        end
        
        %% UPDATE REFINEMENT
        nodalValues = find(d(1:size(X,1)) > refinementValue);
        if (nodalValues)
            elementsToRefine = [];
            for i = 1:length(nodalValues)
                a = nodalValues(i);
                [elements,aux] = find(T == a);
                elementsToRefine = [elementsToRefine;elements];
            end
            elementsToRefine = unique(elementsToRefine);
            elementsToRefine = setdiff(elementsToRefine,refinedElements);
            
            if (elementsToRefine)
               fprintf('\nMesh adaptivity in progress... \n')
                %% Update refined zone
                [Xref,Tref,NitscheFaces,nodesCCDStd,nodesCCDRef,nodesCCDStdCorrected,nonStdNodes,refinedElements,d,ux,uy,H_previous,~,~] = updateRefinedZoneNitscheTwoMeshes(X,T,Xref,Tref,F,infoFaces,NitscheFaces,refinedElements,elementsToRefine,referenceElementRef,referenceElementStd,nodesCCDStd,nodesCCDRef,d,ux,uy,H_previous,test,[],[]);
                d_previous = d;
                H = computeH(X,T,Xref,Tref,refinedElements,referenceElementStd,referenceElementRef,ux,uy,H_previous,lambda,mu);
                
                nonStdNodes_xy = [nonStdNodes,nonStdNodes+size(X,1)+size(Xref,1)];
                stdNodes = setdiff(1:size(X,1),nonStdNodes);
                nDOFStd = length(stdNodes); nDOFRef = size(Xref,1);
                nDOF = (nDOFStd+nDOFRef);
                nodesCCD = [nodesCCDStdCorrected,nodesCCDRef+nDOFStd,nodesCCDStdCorrected+nDOF,nodesCCDRef+nDOF+nDOFStd];
                notCCD = setdiff(1:(2*nDOF),nodesCCD); %actual degrees of freedom (not boundary nodes)
                
                nOfRefEls = length(refinedElements); nOfStdEls = nOfElements-nOfRefEls;
            end
        end
    end
    
    fprintf('Solve Converged!\n') 
    stoppingCriterion = 0;
    H_previous = H;
    
    %% LOAD
    lagrangemultipliers = Kdn*sol+Kdd*uCCD;
    caresdaltx = find(uCCD>1.e-5);
    load = sum(lagrangemultipliers(caresdaltx));
    loads(ind) = load; 
    
    %% SAVE RESULTS
    if(ismembertol(k,saveDisplacements))
        structName =  [resultsPath,'k',num2str(k),'.mat'];
        struct.increments = increments;
        struct.ind = ind;
        struct.d = d;
        struct.u = [ux,uy] ;
        struct.H = H;
        struct.loads = loads;
        % infoRefinement
        struct.Xref = Xref; struct.Tref = Tref;
        struct.NitscheFaces = NitscheFaces;
        struct.refinedElements = refinedElements;
        save(structName,'struct');
    end
end