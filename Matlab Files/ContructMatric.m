% Generate implicit matrix for 2D heat equation with Dirichlet in bottom and right and Neumann in top and left
% Assumes dx = dy
% Parameters:
%    ----------
%    nx   : int
%        number of discretization points in x
%    ny   : int
%        number of discretization points in y
%    sigma: float
%        alpha*dt/dx
        
%    Returns:
%    -------
%    matrixA: 2D array of floats
%        Matrix of implicit 2D heat equation
function matrixA = ContructMatric(nx,ny,sigma)
    x = round((nx-2)*(ny-2));
    matrixA = zeros(x,x);
    row_number = 1;
    for j=1:ny                                                    % looping
        for i=1:nx
            % construct the corners
            if i == 1 && j== 1 % Bottom left corner (Dirichlet down and left)
                matrixA(row_number,row_number) = 1/sigma +4;
                matrixA(row_number,row_number+1) = -1;
                matrixA(row_number,row_number+nx-2) = -1;
            elseif i == nx-2 && j == 1 % Bottom right corner (Dirichlet down, Neumann right)
                matrixA(row_number,row_number) = 1/sigma +3;
                matrixA(row_number,row_number+1) = -1;
                matrixA(row_number,row_number+nx-2) = -1;
            elseif i ==1 && j == ny-2 % Bottom right corner (Dirichlet down, Neumann right)
                matrixA(row_number,row_number) = 1/sigma + 3;
                matrixA(row_number,row_number+1) = -1;
                matrixA(row_number,row_number+nx-2) = -1;
            elseif i == nx-2 && j == ny-2 % Bottom right corner (Dirichlet down, Neumann right)
                matrixA(row_number,row_number) = 1/sigma + 2;
                matrixA(row_number,row_number+1) = -1;
                matrixA(row_number,row_number+nx-2) = -1;
            %construc the sides 
            elseif i == 1 % Left boundary (Dirichlet)
                matrixA(row_number,row_number) = 1/sigma+4; % Set diagonal
                matrixA(row_number,row_number+1) = -1;      % fetch i+1
                matrixA(row_number,row_number+nx-2) = -1;   % fetch j+1
                if row_number-(nx-2)>0
                    matrixA(row_number,row_number-(nx-2)) = -1; % fetch j-1
                end
            elseif i ==nx-2 % Right boundary (Neumann)
                matrixA(row_number,row_number) = 1/sigma+3; % Set diagonal
                matrixA(row_number,row_number+1) = -1;      % fetch i+1
                matrixA(row_number,row_number+nx-2) = -1;   % fetch j+1
                if row_number-(nx-2)>0
                    matrixA(row_number,row_number-(nx-2)) = -1; % fetch j-1
                end
            elseif j == 1 %Bottom boundary (Dirichlet)
                matrixA(row_number,row_number) = 1/sigma+4; % Set diagonal
                matrixA(row_number,row_number+1) = -1;      % fetch i+1
                matrixA(row_number,row_number+nx-2) = -1;   % fetch j+1
                if row_number-(nx-2)>0
                    matrixA(row_number,row_number-(nx-2)) = -1; % fetch j-1
                end
            elseif j == ny-2
                matrixA(row_number,row_number) = 1/sigma+3; % Set diagonal
                matrixA(row_number,row_number+1) = -1;      % fetch i+1
                matrixA(row_number,row_number+nx-2) = -1;   % fetch j+1
                matrixA(row_number,row_number-(nx-2)) = -1; % fetch j-1
            else % interior points
                matrixA(row_number,row_number) = 1/sigma+4; % Set diagonal
                matrixA(row_number,row_number+1) = -1;      % fetch i+1
                matrixA(row_number,row_number-1) = -1;      % fetch i+1
                matrixA(row_number,row_number+nx-2) = -1;   % fetch j+1
                matrixA(row_number,row_number-(nx-2)) = -1; % fetch j-1
            end
            row_number = row_number +1;
        end
    end