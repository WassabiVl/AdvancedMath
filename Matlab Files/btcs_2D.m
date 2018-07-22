function MatrixT = btcs_2D(T, A, nt, sigma, T_bc, nx, ny, dt, j_mid,i_mid,T_avg)
    for j=1:nt
       Tn = T;
       b = ContructRHS(nx, ny, sigma, Tn, T_bc);
       T_interior = gauss2(A,b);
       MatrixT = map_1Dto2D(nx, ny, T_interior, T_bc);
       if T(j_mid,i_mid)>=T_avg
           fprintf('Center of plate reached at time %d ', dt)
           break
       end
    end
end






