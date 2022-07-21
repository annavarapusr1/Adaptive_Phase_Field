function L2Norm = computeL2NormTwoMeshes(nonStdNodes,referenceElementStd,referenceElementRef,refinedElements,X,T,Xref,Tref,u,u0,varargin)

% The function computeL2Norm allows to compute the L2 Norm using an
% analytical function handle. This function has to depend on x-y
% coordinates.
%
% Input:
%  referenceElement: information of the reference element
%  X: nodal coordinates
%  T: connectivity matrix
%  u: FEM solution (nodal)
%  u0: analytical function handle. This function has to have the x-y
%      coordinates matrix [xi yi] as a first input argument.
%  varargin: ordered list of input arguments needed to execute u0 but x-y
%            coordinates.
% Output:
%  L2Norm: computed L2 norm = sqrt(integral[(u-u0)^2])/sqrt(integral[u0^2]) over the domain defined by
%          T,X

nOfElements = size(T,1);

standardElements = setdiff(1:nOfElements,refinedElements);
nOfRefEls = length(refinedElements); nOfStdEls = length(standardElements);
nDOFStd = size(X,1)-length(nonStdNodes);
uStd = u(1:nDOFStd,:); uRef = u(nDOFStd+1:end,:);

sol = zeros(size(X,1)+size(Xref,1),size(u,2));
sol(setdiff(1:size(X,1),nonStdNodes),:) = uStd;
sol(size(X,1)+1:end,:) = uRef;

%Loop in 2D elements
L2Norm = 0; L2NormSolution = 0;
for i = 1:nOfStdEls
    Te = T(standardElements(i),:);
    Xe = X(Te,:);
    ue = sol(Te,:);
    [elemError,elemNormSol] = ElementalL2Norm(Xe,referenceElementStd,ue,u0,varargin{:});
    L2Norm = L2Norm + elemError;
    L2NormSolution = L2NormSolution + elemNormSol;
end
for i = 1:nOfRefEls
    Te = Tref(i,:);
    Xe = Xref(Te,:);
    ue = sol(Te+size(X,1),:);
    [elemError,elemNormSol] = ElementalL2Norm(Xe,referenceElementRef,ue,u0,varargin{:});
    L2Norm = L2Norm + elemError;
    L2NormSolution = L2NormSolution + elemNormSol;
end
L2Norm = sqrt(L2Norm)/sqrt(L2NormSolution);


%_______________________________________________________________________
function [elemL2NormError,elemL2NormSol] = ElementalL2Norm(Xe,referenceElement,ue,u0,varargin)

if (size(Xe,2) == 2) %2D
    N = referenceElement.N;
    Nxi = referenceElement.Nxi;
    Neta = referenceElement.Neta;
    J11 = Nxi*Xe(:,1); J12 = Nxi*Xe(:,2); 
    J21 = Neta*Xe(:,1); J22 = Neta*Xe(:,2);
    detJ = J11.*J22-J12.*J21;
    dvolu = referenceElement.IPweights.*detJ;
elseif (size(Xe,2) == 3) %3D
    N = referenceElement.N;
    Nxi = referenceElement.Nxi;
    Neta = referenceElement.Neta;
    Nzeta = referenceElement.Nzeta;
    J11 = Nxi*Xe(:,1); J12 = Nxi*Xe(:,2); J13=Nxi*Xe(:,3);
    J21 = Neta*Xe(:,1); J22 = Neta*Xe(:,2); J23 = Neta*Xe(:,3);
    J31 = Nzeta*Xe(:,1); J32 = Nzeta*Xe(:,2); J33 = Nzeta*Xe(:,3);
    detJ = J11.*J22.*J33-J11.*J23.*J32-J12.*J21.*J33+J12.*J23.*J31+J13.*J21.*J32-J13.*J22.*J31;
    dvolu = referenceElement.IPweights.*detJ;
end

u0_g = u0(N*Xe,varargin{:});
ue_g = N*ue;
elemL2NormError = sum(sum((ue_g-u0_g).^2,2).*dvolu);
elemL2NormSol = sum(sum((u0_g).^2,2).*dvolu);

