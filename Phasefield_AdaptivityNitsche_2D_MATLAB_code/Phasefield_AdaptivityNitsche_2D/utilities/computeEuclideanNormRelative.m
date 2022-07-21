function EuclideanNormRelative = computeEuclideanNormRelative(u,uprevious)

% L2NormRelative = |u-u_previous|/|u|
%
% Input:
%  u: FEM solution (nodal)
%  u_previous: FEM solution from previous iteration
%
% Output:
%  EuclideanNormRelative 

normDifference = sqrt(sum((u-uprevious).^2));
normU = sqrt(sum(u.^2));
EuclideanNormRelative = normDifference/normU;



