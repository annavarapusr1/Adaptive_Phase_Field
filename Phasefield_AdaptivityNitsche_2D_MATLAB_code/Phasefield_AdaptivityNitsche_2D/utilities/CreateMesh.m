function [X,T] = CreateMesh(elem,nen,dom,npx,npy)
% 
% [X,T] = CreateMesh(elem,nen,dom,npx,npy)
% Creates the topology of an structured and uniform mesh
% over a rectangular domain [x1,x2]x[y1,y2]
%
% Input:    
%   elem: 0 for quadrilatera, 1 for triangles
%   dom = [x1,x2,y1,y2]:  vertices' coordinates
%   npx,npy: number of nodes on each direction 
% Output:   
%   X:  nodal coordinates
%   T:  connectivities


x1 = dom(1); x2 = dom(2); 
y1 = dom(3); y2 = dom(4); 

% Allocate space for the nodal coordinates matrix
X = zeros((npx)*(npy),2);
xs = linspace(x1,x2,npx)'; 
unos = ones(npx,1);
% Nodes' coordinates
yys = linspace(y1,y2,npy);
for i=1:npy
    ys = yys(i)*unos; 
    posi = (i-1)*(npx)+1:i*(npx); 
    X(posi,:)=[xs,ys];
end

% Number of divisions in each direction
mx = npx-1; my = npy-1;


% Connectivities
if elem == 0            % Quadrilaterals
    if nen == 4         % Q1
        nx = mx; ny = my; 
        T = zeros(nx*ny,4);
        for i=1:ny
            for j=1:nx
                ielem = (i-1)*nx+j;
                inode = (i-1)*(npx)+j;
                T(ielem,:) = [inode   inode+1   inode+npx+1   inode+npx];
            end   
        end
    elseif nen == 9     % Q2
        nx = mx/2; ny = my/2; 
        if abs(round(nx)-nx) > 1e-5 || abs(round(ny)-ny) > 1e-5
            error('Number of nodes in X or Y direction is not odd')
        end
        
        T = zeros(nx*ny,9);
        for i=1:ny
            for j=1:nx
                ielem = (i-1)*nx + j;
                inode = (i-1)*2*npx + 2*(j-1) + 1;
                nodes_aux = [inode+(0:2)  inode+npx+(0:2)  inode+2*npx+(0:2)];
                T(ielem,:) = nodes_aux([1  3  9  7  2  6  8  4  5]); 
            end
        end
    elseif nen == 16     % Q3
        nx = mx/3; ny = my/3; 
        if abs(round(nx)-nx) > 1e-5 || abs(round(ny)-ny) > 1e-5
            error('Number of nodes in X or Y direction is not odd')
        end
        
        T = zeros(nx*ny,16);
        for i=1:ny
            for j=1:nx
                ielem = (i-1)*nx + j;
                inode = (i-1)*3*npx + 3*(j-1) + 1;
                nodes_aux = [inode+(0:3)  inode+npx+(0:3)  ...
                    inode+2*npx+(0:3)  inode+3*npx+(0:3)];
                T(ielem,:) = nodes_aux([1  4  16  13  2  3  8  12  ...
                    15  14  9  5  6  7  11  10]); 
            end
        end
    elseif nen == 25     % Q4
        nx = mx/4; ny = my/4; 
        if abs(round(nx)-nx) > 1e-5 || abs(round(ny)-ny) > 1e-5
            error('Number of nodes in X or Y direction is not odd')
        end
        
        T = zeros(nx*ny,25);
        for i=1:ny
            for j=1:nx
                ielem = (i-1)*nx + j;
                inode = (i-1)*4*npx + 4*(j-1) + 1;
                nodes_aux = [inode+(0:4)  inode+npx+(0:4)  ...
                    inode+2*npx+(0:4)  inode+3*npx+(0:4)  inode+4*npx+(0:4)]; 
                T(ielem,:) = nodes_aux([1  5  25  21  2  3  4  10  15  20 ...
                    24  23  22  16  11  6  7  9  19  17  8  14  18  12  13]); 
            end
        end
    else 
        error('Unavailable quadrilatera');
    end
elseif elem == 1        % Triangles
    if nen == 3             % P1
        nx = mx; ny = my; 
        T = zeros(nx*ny,3);
        for i=1:ny
            for j=1:nx
                ielem = 2*((i-1)*nx+j)-1;
                inode = (i-1)*(npx)+j;
                nodes_qua = [inode, inode+1, inode+npx+1, inode+npx]; 
                if rem(i+j,2) == 0
                    T(ielem,:) = nodes_qua([1,2,3]); 
                    T(ielem+1,:) = nodes_qua([1,3,4]); 
                else
                    T(ielem,:) = nodes_qua([1,2,4]); 
                    T(ielem+1,:) = nodes_qua([2,3,4]);
                end
            end   
        end
        % Modification of left lower and right upper corner elements to avoid them 
        % having all their nodes on the boundary
        T(1,:) = [1  npx+2   npx+1];
        T(2,:) = [1    2     npx+2];
        aux = size(T,1);
        T(aux,:) = [npx*ny-1    npx*npy   npx*npy-1];
        T(aux-1,:)   = [npx*ny-1    npx*ny    npx*npy];
    elseif nen == 4         % P1+  (Note that "bubble" coordinates are not considered)
        nx = mx; ny = my; 
        T = zeros(nx*ny,4);
        for i=1:ny
            for j=1:nx
                ielem = 2*((i-1)*nx+j)-1;
                inode = (i-1)*(npx)+j;
                n_ad = npx*npy + 2*((i-1)*nx+j)-1;
                nodes_qua = [inode, inode+1, inode+npx+1, inode+npx]; 
                if rem(i+j,2) == 0
                    T(ielem,:) = [nodes_qua([1,2,3]), n_ad]; 
                    T(ielem+1,:) = [nodes_qua([1,3,4]), n_ad+1]; 
                else
                    T(ielem,:) = [nodes_qua([1,2,4]), n_ad]; 
                    T(ielem+1,:) = [nodes_qua([2,3,4]), n_ad+1];
                end
            end   
        end
        % Modification of left lower and right upper corner elements to avoid them 
        % having all their nodes on the boundary
        aux = size(T,1);
        T(1,:) = [1  npx+2   npx+1  npx*npy+1];
        T(2,:) = [1    2     npx+2  npx*npy+2];
        T(aux-1,:) = [npx*my-1    npx*npy   npx*npy-1   npx*npy+2*mx*my-1];
        T(aux,:)   = [npx*my-1    npx*my    npx*npy   npx*npy+2*mx*my];
    elseif nen == 6         % P2
        nx = mx/2; ny = my/2; 
        if abs(round(nx)-nx) > 1e-5 || abs(round(ny)-ny) > 1e-5
            error('Number of nodes in X or Y direction is not odd')
        end
        
        T = zeros(2*nx*ny,6); 
        for i=1:ny
            for j=1:nx
                ielem=2*((i-1)*nx+j)-1;
                inode=(i-1)*2*(npx)+2*(j-1)+1;
                nodes_aux = [inode+(0:2)  inode+npx+(0:2)  inode+2*npx+(0:2)];
                if rem(i+j,2) == 0
                    T(ielem,:) = nodes_aux([1  3  7  2  5  4]);
                    T(ielem+1,:) = nodes_aux([3  9  7  6  8  5]);
                else
                    T(ielem,:) = nodes_aux([1  3  9  2  6  5]);
                    T(ielem+1,:) = nodes_aux([1  9  7  5  8  4]);
                end
            end    
        end
        % Modification of left lower and right upper corner elements to avoid them 
        % having all their nodes on the boundary
        inode = 1; 
        nodes_aux = [inode+(0:2)  inode+npx+(0:2)  inode+2*npx+(0:2)];
        T(1,:) = nodes_aux([1  9  7  5  8  4]);
        T(2,:) = nodes_aux([1  3  9  2  6  5]);
        
        ielem  = size(T,1)-1;
        inode = npx*(npy-2)-2; 
        nodes_aux = [inode+(0:2)  inode+npx+(0:2)  inode+2*npx+(0:2)];
        T(ielem,:) = nodes_aux([1  9  7  5  8  4]);
        T(ielem+1,:) = nodes_aux([1  3  9  2  6  5]);
    elseif nen == 7
        nx = mx/2; ny = my/2; 
        if abs(round(nx)-nx) > 1e-5 || abs(round(ny)-ny) > 1e-5
            error('Number of nodes in X or Y direction is not odd')
        end
        
        T = zeros(2*nx*ny,7); 
        for i=1:ny
            for j=1:nx
                ielem=2*((i-1)*nx+j)-1;
                inode=(i-1)*2*(npx)+2*(j-1)+1;
                n_ad = npx*npy + 2*((i-1)*nx+j)-1;
                nodes_aux = [inode+(0:2)  inode+npx+(0:2)  inode+2*npx+(0:2)];
                if rem(i+j,2) == 0
                T(ielem,:) = [nodes_aux([1  3  7  2  5  4])  n_ad];
                T(ielem+1,:) = [nodes_aux([3  9  7  6  8  5])  n_ad+1];
                else
                    T(ielem,:) = [nodes_aux([1  3  9  2  6  5]), n_ad];
                    T(ielem+1,:) = [nodes_aux([1  9  7  5  8  4]), n_ad+1];
                end
            end    
        end
        % Modification of left lower and right upper corner elements to avoid them 
        % having all their nodes on the boundary
        inode = 1; 
        nodes_aux = [inode+(0:2)  inode+npx+(0:2)  inode+2*npx+(0:2)];
        T(1,:) = [nodes_aux([1  9  7  5  8  4]), npx*npy+1];
        T(2,:) = [nodes_aux([1  3  9  2  6  5]), npx*npy+2];

        ielem  = size(T,1)-1;
        inode = npx*(npy-2)-2; 
        nodes_aux = [inode+(0:2)  inode+npx+(0:2)  inode+2*npx+(0:2)];
        T(ielem,:) = [nodes_aux([1  9  7  5  8  4]), npx*npy+ielem];
        T(ielem+1,:) = [nodes_aux([1  3  9  2  6  5]), npx*npy+ielem+1];
    elseif nen == 10
        nx = mx/3; ny = my/3; 
        if abs(round(nx)-nx) > 1e-5 || abs(round(ny)-ny) > 1e-5
            error('Number of nodes in X or Y direction is not odd')
        end
        
        T = zeros(2*nx*ny,10); 
        for i=1:ny
            for j=1:nx
                ielem = 2*((i-1)*nx+j)-1;
                inode = (i-1)*3*npx + 3*(j-1) + 1;
                nodes_aux = [inode+(0:3)  inode+npx+(0:3)  ...
                    inode+2*npx+(0:3)  inode+3*npx+(0:3)];
                if rem(i+j,2) == 0
                    T(ielem,:)=nodes_aux([1  4  13  2  3  7  10  9  5  6]); 
                    T(ielem+1,:)=nodes_aux([4  16  13  8  12  15  14  10  7  11]); 
                else
                    T(ielem,:)=nodes_aux([1  4  16  2  3  8  12  11  6  7]); 
                    T(ielem+1,:)=nodes_aux([1  16  13  6  11  15  14  9  5  10]); 
                end
            end
        end
        % Modification of left lower and right upper corner elements to avoid them 
        % having all their nodes on the boundary
        inode = 1; 
        nodes_aux = [inode+(0:3)  inode+npx+(0:3)  ...
            inode+2*npx+(0:3)  inode+3*npx+(0:3)];
        T(1,:) = nodes_aux([1  16  13  6  11  15  14  9  5  10]); 
        T(2,:) = nodes_aux([1  4  16  2  3  8  12  11  6  7]); 

        inode = npx*(npy-3)-3; 
        ielem = 2*nx*ny-1;
        nodes_aux = [inode+(0:3)  inode+npx+(0:3)  ...
            inode+2*npx+(0:3)  inode+3*npx+(0:3)];
        T(ielem,:) = nodes_aux([1  16  13  6  11  15  14  9  5  10]); 
        T(ielem+1,:) = nodes_aux([1  4  16  2  3  8  12  11  6  7]); 
    elseif nen == 15
        nx = mx/4; ny = my/4; 
        if abs(round(nx)-nx) > 1e-5 || abs(round(ny)-ny) > 1e-5
            error('Number of nodes in X or Y direction is not odd')
        end
        
        T = zeros(2*nx*ny,15); 
        for i=1:ny
            for j=1:nx
                ielem = 2*((i-1)*nx+j)-1;
                inode = (i-1)*4*npx + 4*(j-1) + 1;
                nodes_aux = [inode+(0:4)  inode+npx+(0:4)  ...
                    inode+2*npx+(0:4)  inode+3*npx+(0:4)  inode+4*npx+(0:4)]; 
                if rem(i+j,2) == 0
                    T(ielem,:) = nodes_aux([1  5  21  2  3  4  9  13 ...
                        17  16  11  6  7  8  12]); 
                    T(ielem+1,:) = nodes_aux([5  25  21  10  15  20  24  23 ...
                        22  17  13  9  14  19  18]); 
                else
                    T(ielem,:) = nodes_aux([1  5  25  2  3  4  10  15  20 ...
                        19  13  7  8  9  14]); 
                    T(ielem+1,:) = nodes_aux([1  25  21  7  13  19  ...
                        24  23  22  16  11  6   12  18  17]); 
                end
            end
        end
        % Modification of left lower and right upper corner elements to avoid them 
        % having all their nodes on the boundary
        inode = 1; 
        nodes_aux = [inode+(0:4)  inode+npx+(0:4)  ...
            inode+2*npx+(0:4)  inode+3*npx+(0:4)  inode+4*npx+(0:4)]; 
        T(1,:) = nodes_aux([1  25  21  7  13  19  24  23  22  16  11  6  12  18  17]); 
        T(2,:) = nodes_aux([1  5  25  2  3  4  10  15  20  19  13  7  8  9  14]);     
   
        inode = npx*(npy-4)-4; 
        ielem = 2*nx*ny-1; 
        nodes_aux = [inode+(0:4)  inode+npx+(0:4)  ...
            inode+2*npx+(0:4)  inode+3*npx+(0:4)  inode+4*npx+(0:4)]; 
        T(ielem,:) = nodes_aux([1  25  21  7  13  19  24  23  22  16  11  6  12  18  17]); 
        T(ielem+1,:) = nodes_aux([1  5  25  2  3  4  10  15  20  19  13  7  8  9  14]);     
    else
        error('Unavailable triangle')
    end
end   

