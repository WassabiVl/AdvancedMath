% Generates right-hand side for 2D implicit heat equation with Dirichlet in bottom and left and Neumann in top and right
%        Assumes dx=dy, Neumann BCs = 0, and constant Dirichlet BCs
% Assumes dx = dy
% Parameters:
%    ----------
%    nx   : int
%        number of discretization points in x
%    ny   : int
%        number of discretization points in y
%    sigma: float
%        alpha*dt/dx
%  T    : array of float
%            Temperature in current time step
%        T_bc : float
%            Temperature in Dirichlet BC
%        
%        Returns:
%        -------
%        RHS  : array of float
%            Right hand side of 2D implicit heat equation
function RHS = ContructRHS(nx,ny,sigma1,sigma2, T, T_bc)
    x = round((nx-2)*(ny-2));
    RHS = zeros(x,1);
    row_number = 1;
    for j=1:ny-1                                                    % looping
        for i=1:nx-1
            if i<= nx/2
                sigma = sigma1; % use plate 1 sigma
            else
                sigma = sigma2; % use plate 2 sigma
            end
            % construct the corners
            if i == 1 && j== 1 % Bottom left corner (Dirichlet down and left)
                RHS(row_number) = T(j,i)*1/sigma + 2*T_bc;
            elseif i == nx-2 && j == 1 % Bottom right corner (Dirichlet down, Neumann right)
                RHS(row_number) = T(j,i)*1/sigma + T_bc;
            elseif i ==1 && j == ny-2 % Bottom right corner (Dirichlet down, Neumann right)
                 RHS(row_number) = T(j,i)*1/sigma + T_bc;
            elseif i == nx-2 && j == ny-2 % Bottom right corner (Dirichlet down, Neumann right)
                 RHS(row_number) = T(j,i)*1/sigma;
            %construc the sides 
            elseif i == 1 % Left boundary (Dirichlet)
                 RHS(row_number) = T(j,i)*1/sigma + T_bc;
            elseif i ==nx-2 % Right boundary (Neumann)
                 RHS(row_number) = T(j,i)*1/sigma+ T_bc;
            elseif j == 1 %Bottom boundary (Dirichlet)
                 RHS(row_number) = T(j,i)*1/sigma + T_bc;
            elseif j == ny-2
                 RHS(row_number) = T(j,i)*1/sigma;
            else % interior points
                RHS(row_number) = T(j,i)*1/sigma;
            end
            row_number = row_number + 1;
        end
    end
end
    