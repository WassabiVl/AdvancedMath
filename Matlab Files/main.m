% 2D Heat Equation with 2 Plates Using Finite Difference Method
%  Wael Al Atrash 118011
%   Advanced Numerical Mathmatics
%    Bauhaus University




% assumptions
% plate 1 and plate 2 are insulated
% plate 1 and plate 2 have the same dimensions

% Properties of Plate 1
name=('Plate 1'); % Aluminum example
Hc1 = 910; % heat capacity of Plate 1
p1 = 2700.0; % density of Plate 1
K1 = 204.3; % thermal conductivity of Plate 1
Alpha1 = K1/(p1*Hc1); % thermal diffusion of Plate 1


% Properties of Plate 2
name2=('Plate 2'); %copper example
Hc2 = 390; % heat capacity of Plate 2
p2 = 8940; % density of Plate 2
K2 = 401; % thermal conductivity of Plate 2
Alpha2 = K2/(p2*Hc2);

% Common properties between both
NodesL = 20; % nodes in x direction
NodesH = 20; % nodes in y direction
L = 1; H = 0.5; % Length and height of each plate in meters
T_initial= 0 ;                   % Initial temperature in all nodes ( the whole plate )
T_east   =  1000;                   % temperature on the upper side ( at y=0  "Dirichlet Conditions" )
T_west   =  1000;                   % temperature on the lower side ( at y=H "Dirichlet Conditions" )
T_north  = 1000 ;                   % temperature on the left  side ( at x=0  "Dirichlet Conditions" )
T_south  = 1000 ;                   % temperature on the right side ( at x=L "Dirichlet Conditions" ) 
T_bc = 1000;                        % Dirichlet Boundary Condidtions for implicit

t_end=1000 ;                        % final time for visual simulation (sec)
dt=0.6 ;                           % time step (1 sec)


tolerance = 0.5;                   % tolerance for numerical simulation (0.5 deg Celsius)
tolerance_ss=0.001;                % tolerance for steady state section (0.1 deg Celsius)


k=1;                                                                       % iteration counter
err_SS_max(k)=1;                                                                  % initial error
err_SS_min(k)=1;                                                                  % initial error
dx=L/NodesL;                                                                  % delta x
dy=H/NodesH; 
dt1 = (0.25 * min([dx/2 dy])^2) / Alpha1;                            % time step for plate 1
dt2 = (0.25 * min([dx/2 dy])^2) / Alpha2;  
sigma = Alpha1*dt1/(dx/2) + Alpha2*dt2/(dx/2);
sigma1 = Alpha1*dt1/(dx/2);
sigma2 = Alpha2*dt2/(dx/2);% time step for plate 2% delta y
n_time=round(t_end/(dt1+dt2));                                                           % number of iterations for time
T_max=max([T_east T_west T_north T_south T_initial]);                      % Max T to set axes limits in plotting
T_min=min([T_east T_west T_north T_south T_initial]);                      % Min T to set axes limits in plotting
A = ContructMatric(NodesL+2, NodesH+2, sigma1,sigma2);
mid_i = round(NodesL/2);
mid_j = round(NodesH/2);

Solution_type=questdlg('Which method you want to solve the time derivative with ?','Question','Euler','2nd order Runge-Kutte','Implicit','Euler');                     % solve with 2nd order Runge Kutte in time or 2 to solve with Euler

if dt1 <= 1/(2*Alpha1*((1/dx^2)+(1/dy^2))) || dt2 <= 1/(2*Alpha2*((1/dx^2)+(1/dy^2)))                              % test the stability condition 
else 
    fprintf('Error, the stability condition is not met\nPlease return to "Inputs Section" and choose a "dt" smaller than %f \n',1/(2*Alpha1*((1/dx^2)+(1/dy^2))))
    return
end

message=msgbox('Your computer is now solving the problem, Please wait..... ');    % Busy message 
% ----------------- Initial Conditions for finite difference section ---------------
T=zeros(NodesL+2,NodesH+2,75000);                       % set max iterations 75,000 due to memory limitations (T variable takes maximum 1GB in memory)
T(:,:,1)=T_initial;
T(:,1,1)=T_south;
T(:,2,:)=T_south;
T(:,end-1,1)=T_north;
T(:,end,1)=T_north;                            
T(end,:,1)=T_east;
T(end-1,:,1)=T_east;                             
T(1,:,1)=T_west;
T(2,:,1)=T_west;
T_avg = min(mean2(T(:,:,1)))/2;

%% 3-Test matrix for the solution

% rule of thermodynamics, if no external factors is affecting the system at time >= 0, the end state of the plate be an mean of all the tempretures
% and not a steady state solution
Tss(1:NodesL+2,1:NodesH+2) = T_avg;
x=zeros(1,NodesL+2);y=zeros(1,NodesH+2);            %Generate the plate
for i = 1:NodesL+2                 
x(i) =(i-1)*dx; 
end
for i = 1:NodesH+2                 
y(i) =(i-1)*dy; 
end
[x,y] = meshgrid(x,y);

%% 4- Finite difference section (Using 2nd order Runge Kutte or Euler in time or forward method)

k=1;
switch Solution_type
    case '2nd order Runge-Kutte'
        err_R_k_max(k)=100;                            % initial error
        err_R_k_min(k)=100;                            % initial error
        while err_R_k_max(k)>=tolerance || err_R_k_min(k)>=tolerance
          for i=2:NodesL+1
            for j=2:NodesH+1
                if i<=NodesL/2
                    k1=Alpha1*(((T(i-1,j,k)-2*T(i,j,k)+T(i+1,j,k))/dx^2)+((T(i,j-1,k)-2*T(i,j,k)+T(i,j+1,k))/dy^2));
                    Tk=T(:,:,k)+k1*dt1;
                    k2=Alpha1*(((Tk(i-1,j)-2*Tk(i,j)+Tk(i+1,j))/dx^2)+((Tk(i,j-1)-2*Tk(i,j)+Tk(i,j+1))/dy^2));
                    T(i,j,k+1) =T(i,j,k)+(dt1/2)*(k1+k2);
                else
                    k1=Alpha2*(((T(i-1,j,k)-2*T(i,j,k)+T(i+1,j,k))/dx^2)+((T(i,j-1,k)-2*T(i,j,k)+T(i,j+1,k))/dy^2));
                    Tk=T(:,:,k)+k1*dt2;
                    k2=Alpha2*(((Tk(i-1,j)-2*Tk(i,j)+Tk(i+1,j))/dx^2)+((Tk(i,j-1)-2*Tk(i,j)+Tk(i,j+1))/dy^2));
                    T(i,j,k+1) =T(i,j,k)+(dt2/2)*(k1+k2);
                end
            end
          end
          k=k+1;
          % reinforce Neumann boundary conditions insulation
          disp(T(:,:,k));
          T(:,end,k) = T(:,end-1,k); T(:,end,k+1) = T(:,end,k);
          T(:,1,k) =T(:,2,k); T(:,1,k+1) = T(:,1,k);
          T(end,:,k) = T(end-1,:,k); T(end,:,k+1) = T(end,:,k);
          T(1,:,k) = T(2,:,k); T(1,:,k+1) = T(1,:,k);
          disp(T(:,:,k));
          err_R_k_max(k)=abs(max(max(T(:,:,k)-Tss)));        %calculate error
          err_R_k_min(k)=abs(min(min(T(:,:,k)-Tss)));        %calculate error
          subplot(1,1,1)
           surf(x,y,T(:,:,k));
           title(sprintf('Temperature at time : %i seconds ',round(k*((dt1+dt2)/2))));
           cb=colorbar;
            caxis([T_min T_max]);
            view(90,-90);
            xlim([0 L+dx]); xlabel('Length');
            ylim([0 H+dy]); ylabel('Width');
            zlim([T_min T_max]); zlabel('Temprature');
           drawnow;
          
         if  T(mid_i,mid_j,k)>= T_avg || isequal(T(:,:,k),T(:,:,k-1))   % break out of loop
            break
          end
          if  T(mid_i,mid_j,k)>= T_avg || isequal(T(:,:,k),T(:,:,k-1))   % break out of loop
            break
          end
         end

    case'Euler' % with central second order derivative in space, forward first order derivative in time
        err_E_max(k)=100;                            % initial error
        err_E_min(k)=100;                            % initial error
        while err_E_max(k)>=tolerance || err_E_min(k)>=tolerance 
          
            for j=2:NodesH+1
                for i=2:NodesL+1
                    if i<=NodesL/2
                        T(i,j,k+1) = T(i,j,k)+dt1*Alpha1*(((T(i-1,j,k)-2*T(i,j,k)+T(i+1,j,k))/dx^2)+((T(i,j-1,k)-2*T(i,j,k)+T(i,j+1,k))/dy^2));
                    else
                        T(i,j,k+1) = T(i,j,k)+dt2*Alpha2*(((T(i-1,j,k)-2*T(i,j,k)+T(i+1,j,k))/dx^2)+((T(i,j-1,k)-2*T(i,j,k)+T(i,j+1,k))/dy^2));
                    end
                end
            end
          
          k=k+1;
          % reinforce Neumann boundary conditions insulation
           disp(T(:,:,k));
          T(:,end,k) = T(:,end-1,k); T(:,end,k+1) = T(:,end,k);
          T(:,1,k) =T(:,2,k); T(:,1,k+1) = T(:,1,k);
          T(end,:,k) = T(end-1,:,k); T(end,:,k+1) = T(end,:,k);
          T(1,:,k) = T(2,:,k); T(1,:,k+1) = T(1,:,k);
           disp(T(:,:,k));
          err_E_max(k)=abs(max(max(T(:,:,k)-Tss)));        %calculate error
          err_E_min(k)=abs(min(min(T(:,:,k)-Tss)));        %calculate error
          subplot(1,1,1)
          
           surf(x,y,T(:,:,k));
           title(sprintf('Temperature at time : %i seconds ',round(k*((dt1+dt2)/2))));
           cb=colorbar;
            caxis([T_min T_max]);
            view(90,-90);
            xlim([0 L+dx]); xlabel('Length');
            ylim([0 H+dy]); ylabel('Width');
            zlim([T_min T_max]); zlabel('Temprature');
           
           

           drawnow;
          if  T(mid_i,mid_j,k)>= T_avg || isequal(T(:,:,k),T(:,:,k-1))   % break out of loop
            break
          end
          if  T(mid_i,mid_j,k)>= T_avg || isequal(T(:,:,k),T(:,:,k-1))   % break out of loop
            break
          end
        end
        
     
    case'Implicit'
        T(:,:,2) = T(:,:,1);
        T(:,:,3) = T(:,:,1);
        err_E_max(k)=100;                            % initial error
        err_E_min(k)=100;                            % initial error
        while T(mid_i,mid_j,k)<= T_avg  
          b = ContructRHS(NodesL+2, NodesH+2, sigma1,sigma2, T(:,:,k), T_bc);
           T_interior = gauss2(A,b);
           %T_interior = linsolve(A,b);
           MatrixT = map_1Dto2D(NodesL+2, NodesH+2, T_interior, T_bc);
           k=k+1;
           T(:,:,k) = MatrixT; % add the matrix to the 2d matrix collection
           disp(T(:,:,k));
          if T(mid_i,mid_j,k)>= T_avg
                    break
          end
         surf(x,y,T(:,:,k));
           title(sprintf('Temperature at time : %i seconds ',round(k*((dt1+dt2)/2))));
           cb=colorbar;
            caxis([T_min T_max]);
            view(90,-90);
            xlim([0 L+dx]); xlabel('Length');
            ylim([0 H+dy]); ylabel('Width');
            zlim([T_min T_max]); zlabel('Temprature');
           
           

           drawnow;
          if  T(mid_i,mid_j,k)>= T_avg || isequal(T(:,:,k),T(:,:,k-1))   % break out of loop
            break
          end
          if  T(mid_i,mid_j,k)>= T_avg || isequal(T(:,:,k),T(:,:,k-1))   % break out of loop
            break
          end
        end
         
    

    case []
    close(message)
    msgbox('Error, Please re-run the code and choose Euler or 2nd order Runge-Kutte to continue the solution')
    return
end
T=T(:,:,1:k);                                            % delete the unused assigned zero layers
SStime=k*(dt1+dt2)/2;                                             % steady state time
close(message)                                               % close the busy message



