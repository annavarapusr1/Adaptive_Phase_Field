function [coord2D coord1D] = equispacedNodesReferenceElement(nDeg)
coord2D = [];
h = 1/nDeg;
for j = 0:nDeg
    i = (0:nDeg-j)';
    aux = j*ones(size(i));
    coord2D = [coord2D; [i aux]*h];
end
coord2D = 2*coord2D - 1;

coord1D = (-1 : 2/nDeg : 1 )';

% renumber
switch nDeg
    case 3
        renumber = [1 4 10 2 3 7 9 8 5 6];
    case 4
        renumber = [1 5 15 2 3 4 9 12 14 13 10 6 7 8 11];
    case 5
        renumber = [1 6 21 2 3 4 5 11 15 18 20 19 16 12 7 8 9 10 13 14 17];
    case 8
        renumber = [1 9 45 2 3 4 5 6 7 8 17 24 30 35 39 42 44 43 40 36 31 25 18 10 11 12 13 14 15 16 19 20 21 22 23 ...
            26 27 28 29 32 33 34 37 38 41];
    case 9
        renumber = [1 10 55 2:9 19 27 34 40 45 49 52 54 53 50 46 41 35 28 20 11 12:18 21:26 29:33 ...
            36:39 42:44 47 48 51];
    case 11
        renumber = [1 12 78 2:11 23 33 42 50 57 63 68 72 75 77 76 73 69 64 58 51 43 34 24 13:22 25:32 35:41 44:49 52:56 ...
            59:62 65:67 70 71 74];
        
    otherwise
        error('Equispaced reference element not implemented')
end

coord2D = coord2D(renumber,:);
