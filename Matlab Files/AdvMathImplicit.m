% 2D Heat Equation with 2 Plates Using Finite Difference Method
%  Wael Al Atrash 118011
%   Advanced Numerical Mathematics
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
name2=('Plate 2');     % copper example
Hc2 = 390;             % heat capacity of Plate 2
p2 = 8940;             % density of Plate 2
K2 = 401;              % thermal conductivity of Plate 2
Alpha2 = K2/(p2*Hc2);  % thermal difussion of Plate 2

% Common properties between both
NodesL = 20;                     % nodes in x direction
NodesH = 20;                     % nodes in y direction
L = 1; H = 0.5;                  % Length and height of each plate in meters
T_initial= 0 ;                   % Initial temperature in all nodes ( the whole plate )
T_bc = 100;                      % Dirichlet Boundary Condidtions
t_end=10 ;                       % final time for visual simulation (sec)
dt=0.6 ;                         % time step (1 sec)


tolerance = 0.5;                   % tolerance for numerical simulation (0.5 deg Celsius)
tolerance_ss=0.001;                % tolerance for steady state section (0.1 deg Celsius)

k=1;                               % iteration counter
err_SS_max(k)=1;                   % initial error
err_SS_min(k)=1;                   % initial error
dx=L/NodesL;                       % delta x
dy=H/NodesH;                       % delta y
dt1 = (0.25 * min([dx dy])^2) / Alpha1; % time step for plate 1
dt2 = (0.25 * min([dx dy])^2) / Alpha2; % time step for plate 2% delta y
sigma = Alpha1*dt1/(dx/2) + Alpha2*dt2/(dx/2);
sigma1 = Alpha1*dt1/(dx/2);
sigma2 = Alpha2*dt2/(dx/2)

n_time=round(t_end/dt1);                                                           % number of iterations for time
T_max=max([T_bc T_initial]);                      % Max T to set axes limits in plotting
T_min=min([T_bc T_initial]);                      % Min T to set axes limits in plotting
mid_i = round(NodesL/2);
mid_j = round(NodesH/2);

if dt1 <= 1/(2*Alpha1*((1/dx^2)+(1/dy^2))) || dt2 <= 1/(2*Alpha2*((1/dx^2)+(1/dy^2)))                            % test the stability condition
else
    fprintf('Error, the stability condition is not met\nPlease return to "Inputs Section" and fix the data so that "dt1" and "dt2" are smaller than %f, %f \n',1/(2*Alpha1*((1/dx^2)+(1/dy^2)),1/(2*Alpha2*((1/dx^2)+(1/dy^2)))
    return
end
message=msgbox('Your computer is now solving the problem, Please wait..... ');    % Busy message

% ----------------- Initial Conditions for Implected difference section ---------------
T=zeros(NodesL+2,NodesH+2,75000);                       % set max iterations 75,000 due to memory limitations (T variable takes maximum 1GB in memory)
T(:,:,1)=T_initial;
T(:,1,1)=T_bc;
T(:,2,1)=T_bc;
T(:,NodesH+1,1)=T_bc;
T(:,NodesH+2,1)=T_bc;
T(NodesL+1,:,1)=T_bc;
T(NodesL+2,:,1)=T_bc;
T(1,:,1)=T_bc;
T(2,:,1)=T_bc;
T_avg = mean2(T(:,:,1));


% ------------------- Initial Conditions for steady state section -------------------
Tss=zeros(NodesL+2,NodesH+2); Tss2=zeros(NodesL+2,NodesH+2);
Tss(:,1)=T_bc;                Tss2(:,1)=T_bc;
Tss(:,2)=T_bc;                Tss2(:,2)=T_bc;
Tss(:,NodesH+1)=T_bc;         Tss2(:,NodesH+1)=T_bc;
Tss(:,NodesH+2)=T_bc;         Tss2(:,NodesH+2)=T_bc;
Tss(NodesL+1,:)=T_bc;         Tss2(NodesL+1,:)=T_bc;
Tss(NodesL+2,:)=T_bc;         Tss2(NodesL+2,:)=T_bc;
Tss(1,:)=T_bc;                Tss2(1,:)=T_bc;
Tss(2,:)=T_bc;                Tss2(2,:)=T_bc;



%% 3- Steady-State section


while err_SS_max(k)>=tolerance_ss || err_SS_min(k)>=tolerance_ss
    for i=2:NodesL                                                    % looping
        for j=2:NodesH
            Tss2(i,j)=0.25*(Tss(i+1,j)+Tss(i,j+1)+Tss(i-1,j)+Tss(i,j-1));
        end
    end
    k=k+1;                                                        % update k
    err_SS_max(k)=abs(max(max(Tss2-Tss)));                        % calculate error
    err_SS_min(k)=abs(min(min(Tss2-Tss)));                        % calculate error
    Tss=Tss2;                                                     % update T
end


%% 4- Implected difference section (Using Advances diffusion equation in time with Guass elimination or matlab solve linear equation function)

k=1;
A = ContructMatric(NodesL+2, NodesH+2, sigma);  % Generate implicit matrix for 2D heat equation
err_R_k_max(k)=1000;                            % initial error
err_R_k_min(k)=1000;                            % initial error
while T(mid_i,mid_j,k)<= T_avg                  % target the central node to see if its temp is equal to the equalized tempreture of the plates
    %for j=1:t_end
       b = ContructRHS(NodesL+2, NodesH+2, sigma1,sigma2, T(:,:,k), T_bc);
       T_interior = gauss2(A,b);
       %T_interior = linsolve(A,b);
       MatrixT = map_1Dto2D(NodesL+2, NodesH+2, T_interior, T_bc);
       k=k+1;
       T(:,:,k) = MatrixT; % add the matrix to the 2d matrix collection
       disp(T(:,:,k));
       if T(mid_i,mid_j,k)>= T_avg || k > 72.000 || n_time >k
           break
       end
       err_R_k_max(k)=abs(max(max(T(:,:,k)-Tss)));        %calculate error
       err_R_k_min(k)=abs(min(min(T(:,:,k)-Tss)));        %calculate error
       if round(err_R_k_max(k),5)==round(err_R_k_max(k-1),5) && err_R_k_max(k)~= 0      % Test solution convergence
           errordlg('The solution is not converging, Please choose a larger tolerance','tolerance Error');
           close(message)
           return
       end
       if round(err_R_k_min(k),5)==round(err_R_k_min(k-1),5) && err_R_k_min(k)~= 0      % Test solution convergence
           errordlg('The solution is not converging, Please choose a larger tolerance','tolerance Error');
           close(message)
           return
       end
    %end
end
T=T(:,:,1:k);                                            % delete the unused assigned zero layers
SStime=k*dt1;                                             % steady state time
close(message)                                               % close the busy message


%% 5- Printed results section

fprintf('This is the solution of the heat equation through out a plate of dimensions %i X %i \n',L,H);

%% 6- Plotting section

x=zeros(1,NodesL+2);y=zeros(1,NodesH+2);            %Generate the plate
for i = 1:NodesL+2
    x(i) =(i-1)*dx;
end
for i = 1:NodesH+2
    y(i) =(i-1)*dy;
end
[x,y] = meshgrid(x,y);

% %%%            -------------- Constant plot ----------------

subplot(2,2,3)
hold on
title(sprintf('Temperature at iteration number : %i  ',round(k)))
surf(x,y,T(:,:,k))
plot3(  0,  0,T_max,'ko','markerfacecolor','r') % plot red point
plot3(  L/2,  H/2,T_max,'ko','markerfacecolor','g') % plot green point
plot3(L,H,T_max,'ko','markerfacecolor','b') % plot blue point
plot3(  0,  0,T_min,'ko','markerfacecolor','r') % plot red point
plot3(  L/2,  H/2,T_min,'ko','markerfacecolor','g') % plot green point
plot3(L,H,T_min,'ko','markerfacecolor','b') % plot blue point
cb=colorbar;
caxis('auto');
view(90,-90);
xlim([0 L+dx]); xlabel('Length');
ylim([0 H+dy]); ylabel('Width');
zlim([T_min T_max]); zlabel('Temprature');
drawnow
hold off

subplot(2,2,4)
hold on
title(sprintf('Temperature at iteration number : %i  ',round(k)))
scatter(k,T(1,1,k),'ko','markerfacecolor','r');
val=(sprintf('  T =  %0.2f   ',T(1,1,k)));
text(k,T(1,1,k),val,'HorizontalAlignment','Left');
scatter(k,T(floor(NodesL/2),floor(NodesH/2),k),'ko','markerfacecolor','g');
val=(sprintf('  T =  %0.2f   ',T(floor(NodesL/2),floor(NodesH/2),k)));
text(k,T(floor(NodesL/2),floor(NodesH/2),k),val,'HorizontalAlignment','right');
scatter(k,T(NodesL+2,NodesH+2,k),'ko','markerfacecolor','b');
val=(sprintf('  T =  %0.2f   ',T(NodesL+2,NodesH+2,k)));
text(k,T(NodesL+2,NodesH+2,k),val,'HorizontalAlignment','Left');
axis tight; xlabel('Time Iterations');
ylim([T_min T_max]); ylabel('Temperature');
legend('Red Point','Green Point ','Blue Point ','Location','northwest')
drawnow
hold off

%%%             ------------ Animated plot ----------

for j=1:k
    subplot(2,2,1)
    surf(x,y,T(:,:,j))
    hold on
    title(sprintf('Temperature at iteration number : %i ',round(k)))
    plot3(  0,  0,T_max,'ko','markerfacecolor','r') % plot red point
    plot3(  L/2,  H/2,T_max,'ko','markerfacecolor','g') % plot green point
    plot3(L,H,T_max,'ko','markerfacecolor','b') % plot blue point
    plot3(  0,  0,T_min,'ko','markerfacecolor','r') % plot red point
    plot3(  L/2,  H/2,T_min,'ko','markerfacecolor','g') % plot green point
    plot3(L,H,T_min,'ko','markerfacecolor','b') % plot blue point
    cb=colorbar;
    caxis([T_min T_max]);
    view(90,-90);
    xlim([0 L+dx]); xlabel('Length');
    ylim([0 H+dy]); ylabel('Width');
    zlim([T_min T_max]); zlabel('Temprature');
    drawnow
    hold off

    subplot(2,2,2)
    hold on
    title(sprintf('Temperature at time : %i seconds ',round(k)))
    scatter(j,T(1,1,j),'r.');
    scatter(j,T(ceil((NodesL+2)/2),ceil((NodesH+2)/2),j),'g.');
    scatter(j,T(NodesL+2,NodesH+2,j),'b.');
    axis tight; xlabel('Time Iterations');
    axis tight; ylabel('Temperature');
    legend('Red Point','Green Point ','Blue Point ','Location','northwest')
    drawnow
    hold off
end
