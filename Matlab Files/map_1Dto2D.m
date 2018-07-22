% Takes temperatures of solution of linear system, stored in 1D, 
 %   and puts them in a 2D array with the BCs
  %  Valid for constant Dirichlet bottom, left, right, up, and Neumann with zero 
   % flux top and right
        
%    Parameters:
 %   ----------
 %       nx  : int
 %           number of nodes in x direction
 %       ny  : int
 %           number of nodes in y direction
 %       T_1D: array of floats
 %           solution of linear system
 %       T_bc: float
 %           Dirichlet BC
            
 %   Returns:
 %   -------
 %       T: 2D array of float
 %           Temperature stored in 2D array with BCs
 %
function MatrixTA = map_1Dto2D(nx, ny, T_1D, T_bc)
    
    MatrixTA = zeros(ny,nx);
    
    row_number = 1;
    for j=1:ny-1                                                    % looping
        for i=1:nx-1
            MatrixTA(j,i) = T_1D(row_number);
            row_number = row_number+1;
        end
    end
    % Dirichlet BC
    MatrixTA(1,:) = T_bc;
    MatrixTA(end,:) = T_bc;
    MatrixTA(:,1) = T_bc;
    MatrixTA(:,end) = T_bc;
    % Neumann BC
    MatrixTA(end,:) = MatrixTA(end-1,:);
    MatrixTA(1,:) = MatrixTA(2,:);
    MatrixTA(:,end) = MatrixTA(:,end-1);
    MatrixTA(:,1) = MatrixTA(:,2);
end