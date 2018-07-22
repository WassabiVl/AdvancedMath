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
T_east   =  50;                   % temperature on the upper side ( at y=0  "Dirichlet Conditions" )
T_west   =  50;                   % temperature on the lower side ( at y=H "Dirichlet Conditions" )
T_north  = 50 ;                   % temperature on the left  side ( at x=0  "Dirichlet Conditions" )
T_south  = 50 ;                   % temperature on the right side ( at x=L "Dirichlet Conditions" ) 

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
dt2 = (0.25 * min([dx/2 dy])^2) / Alpha2;                    % time step for plate 2% delta y
n_time=round(t_end/(dt1+dt2));                                                           % number of iterations for time
T_max=max([T_east T_west T_north T_south T_initial]);                      % Max T to set axes limits in plotting
T_min=min([T_east T_west T_north T_south T_initial]);                      % Min T to set axes limits in plotting
T_avg = (T_east+ T_west+ T_north+ T_south+ T_initial)/5;
mid_i = round(NodesL/2);
mid_j = round(NodesH/2);
Solution_type=questdlg('Which method you want to solve the time derivative with ?','Question','Euler','2nd order Runge-Kutte','Forward','Euler');                     % solve with 2nd order Runge Kutte in time or 2 to solve with Euler

if dt1 <= 1/(2*Alpha1*((1/dx^2)+(1/dy^2))) || dt2 <= 1/(2*Alpha2*((1/dx^2)+(1/dy^2)))                              % test the stability condition 
else 
    fprintf('Error, the stability condition is not met\nPlease return to "Inputs Section" and choose a "dt" smaller than %f \n',1/(2*Alpha1*((1/dx^2)+(1/dy^2))))
    return
end

message=msgbox('Your computer is now solving the problem, Please wait..... ');    % Busy message 
% ----------------- Initial Conditions for finite difference section ---------------
T=zeros(NodesL+2,NodesH+2,75000);                       % set max iterations 75,000 due to memory limitations (T variable takes maximum 1GB in memory)
T(:,1,:)=T_south;
T(:,NodesH+1,:)=T_north;
T(:,NodesH+2,:)=T_north;                            % Redundant, it has no effect in calculations but is required in plotting section
T(NodesL+1,:,:)=T_east;
T(NodesL+2,:,:)=T_east;                             % Redundant, it has no effect in calculations but is required in plotting section
T(1,:,:)=T_west;
T(:,:,1)=T_initial;
% ------------------- Initial Conditions for steady state section -------------------
Tss=zeros(NodesL+2,NodesH+2);        Tss2=zeros(NodesL+2,NodesH+2);
Tss(:,1)=T_south;            Tss2(:,1)=T_south;
Tss(:,NodesH+1)=T_north;         Tss2(:,NodesH+1)=T_north;
Tss(:,NodesH+2)=T_north;         Tss2(:,NodesH+2)=T_north;             % Redundant, it has no effect in calculations but is required in plotting section
Tss(NodesL+1,:)=T_east;          Tss2(NodesL+1,:)=T_east;
Tss(NodesL+2,:)=T_east;          Tss2(NodesL+2,:)=T_east;              % Redundant, it has no effect in calculations but is required in plotting section
Tss(1,:)=T_west;             Tss2(1,:)=T_west;


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
   

%% 4- Finite difference section (Using 2nd order Runge Kutte or Euler in time or forward method)

k=1;
switch Solution_type
    case '2nd order Runge-Kutte'
    err_R_k_max(k)=100;                            % initial error
    err_R_k_min(k)=100;                            % initial error
    while err_R_k_max(k)>=tolerance || err_R_k_min(k)>=tolerance
      for i=2:NodesL
        for j=2:NodesH
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
            % enforce Neumann boundary conditions insulation
            T(i,end) = T(i,end-1);
            T(end,j) = T(end-1,j);
            if T(mid_i,mid_j)>= T_avg
                dt = (dt1+dt2)/2;
                n_time=round(k/dt);
                disp('reached');
                break
            end    
        end
      end
      k=k+1;
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
     end

    case'Euler' % with central Differincal
        err_E_max(k)=100;                            % initial error
        err_E_min(k)=100;                            % initial error
        while err_E_max(k)>=tolerance || err_E_min(k)>=tolerance
          for i=2:NodesL
            for j=2:NodesH
                if i<=NodesL/2
                    T(i,j,k+1) = T(i,j,k)+dt1*Alpha1*(((T(i-1,j,k)-2*T(i,j,k)+T(i+1,j,k))/dx^2)+((T(i,j-1,k)-2*T(i,j,k)+T(i,j+1,k))/dy^2));
                else
                    T(i,j,k+1) = T(i,j,k)+dt2*Alpha2*(((T(i-1,j,k)-2*T(i,j,k)+T(i+1,j,k))/dx^2)+((T(i,j-1,k)-2*T(i,j,k)+T(i,j+1,k))/dy^2));
                end
                % enforce Neumann boundary conditions insulation
                T(i,end) = T(i,end-1);
                T(end,j) = T(end-1,j);
                if T(mid_i,mid_j)>= T_avg
                    dt = (dt1+dt2)/2;
                    n_time=round(k/dt);
                    break
                end
            end
          end
          k=k+1;
          err_E_max(k)=abs(max(max(T(:,:,k)-Tss)));        %calculate error
          err_E_min(k)=abs(min(min(T(:,:,k)-Tss)));        %calculate error
          if round(err_E_max(k),5)==round(err_E_max(k-1),5) && err_E_max(k)~= 0      % Test solution convergence
            errordlg('The solution is not converging, Please choose a larger tolerance','tolerance Error');
            close(message)
            return
          end
          if round(err_E_min(k),5)==round(err_E_min(k-1),5) && err_E_min(k)~= 0      % Test solution convergence
            errordlg('The solution is not converging, Please choose a larger tolerance','tolerance Error');
            close(message)
            return
          end
        end
     
    case'Forward'
        err_E_max(k)=100;                            % initial error
        err_E_min(k)=100;                            % initial error
        while err_E_max(k)>=tolerance || err_E_min(k)>=tolerance
          for i=2:NodesL
            for j=2:NodesH
                if i<=NodesL/2
                    T(i,j,k+1) = T(i,j,k)+dt1*Alpha1*(((T(i+2,j,k)-2*T(i+1,j,k)+T(i,j,k))/dx^2)+((T(i,j+2,k)-2*T(i,j+1,k)+T(i,j,k))/dy^2));
                else
                    T(i,j,k+1) = T(i,j,k)+dt2*Alpha2*(((T(i+2,j,k)-2*T(i+1,j,k)+T(i,j,k))/dx^2)+((T(i,j+2,k)-2*T(i,j+1,k)+T(i,j,k))/dy^2));
                end
                % enforce Neumann boundary conditions insulation
                T(i,end) = T(i,end-1);
                T(end,j) = T(end-1,j);
                if T(mid_i,mid_j)>= T_avg
                    dt = (dt1+dt2)/2;
                    n_time=round(k/dt);
                    disp('reached');
                    break
                end    
            end
          end
          k=k+1;
          err_E_max(k)=abs(max(max(T(:,:,k)-Tss)));        %calculate error
          err_E_min(k)=abs(min(min(T(:,:,k)-Tss)));        %calculate error
          if round(err_E_max(k),5)==round(err_E_max(k-1),5) && err_E_max(k)~= 0      % Test solution convergence
            errordlg('The solution is not converging, Please choose a larger tolerance','tolerance Error');
            close(message)
            return
          end
          if round(err_E_min(k),5)==round(err_E_min(k-1),5) && err_E_min(k)~= 0      % Test solution convergence
            errordlg('The solution is not converging, Please choose a larger tolerance','tolerance Error');
            close(message)
            return
          end
        end
         
    

    case []
    close(message)
    msgbox('Error, Please re-run the code and choose Euler or 2nd order Runge-Kutte to continue the solution')
    return
end
T=T(:,:,1:k);                                            % delete the unused assigned zero layers
SStime=k*dt1;                                             % steady state time
close(message)                                               % close the busy message


%% 5- Printed results section

fprintf('This is the solution of the heat equation through out a plate of dimensions %i X %i of material %s \n',L,H,name);
fprintf('The solution is  based on "Dirichlet Boundary Conditions" with initial values \n')
fprintf('T(x,0,t)=%i , T(x,%i,t)=%i , T(0,y,t)=%i , T(%i,y,t)=%i , T(x,y,0)=%i \n',T_south,H,T_north,T_west,L,T_east,T_initial)
fprintf('The plate takes %i seconds to reach steady-state temperature with tolerance %0.2f \n',round(SStime),tolerance);
fprintf('Now, Simulation is running with final time %i seconds and step %0.2f second \n',t_end,(dt1+dt2)/2)


%% 6- Plotting section

x=zeros(1,NodesL+2);y=zeros(1,NodesH+2);            %Generate the plate
for i = 1:NodesL+2                 
x(i) =(i-1)*dx; 
end
for i = 1:NodesH+2                 
y(i) =(i-1)*dy; 
end

% %%%            -------------- Constant plot ----------------

subplot(2,2,3)                           
hold on
title(sprintf('Temperature at steady state time : %i seconds ',round(SStime)))
surf(x,y,Tss)
plot3(  L/4,  H/4,T_max,'ko','markerfacecolor','r') % plot red point
plot3(  L/2,  H/2,T_max,'ko','markerfacecolor','g') % plot green point
plot3(3*L/4,3*H/4,T_max,'ko','markerfacecolor','b') % plot blue point
plot3(  L/4,  H/4,T_min,'ko','markerfacecolor','r') % plot red point
plot3(  L/2,  H/2,T_min,'ko','markerfacecolor','g') % plot green point
plot3(3*L/4,3*H/4,T_min,'ko','markerfacecolor','b') % plot blue point
cb=colorbar;
caxis([T_min T_max]);
view(90,-90);
xlim([0 L+dx]); xlabel('Length');
ylim([0 H+dy]); ylabel('Width');
zlim([T_min T_max]); zlabel('Temprature');
drawnow
hold off

 subplot(2,2,4)
hold on
title(sprintf('Temperature at steady state time : %i seconds ',round(SStime)))
scatter(k,Tss(floor(NodesL/4),floor(NodesH/4)),'ko','markerfacecolor','r'); 
val=(sprintf('  T =  %0.2f   ',Tss(floor(NodesL/4),floor(NodesH/4))));
text(k,Tss(floor(NodesL/4),floor(NodesH/4)),val,'HorizontalAlignment','Left');
scatter(k,Tss(floor(NodesL/2),floor(NodesH/2)),'ko','markerfacecolor','g'); 
val=(sprintf('  T =  %0.2f   ',Tss(floor(NodesL/2),floor(NodesH/2))));
text(k,Tss(floor(NodesL/2),floor(NodesH/2)),val,'HorizontalAlignment','right');
scatter(k,Tss(floor(3*NodesL/4),floor(3*NodesH/4)),'ko','markerfacecolor','b'); 
val=(sprintf('  T =  %0.2f   ',Tss(floor(3*NodesL/4),floor(3*NodesH/4))));
text(k,Tss(floor(3*NodesL/4),floor(3*NodesH/4)),val,'HorizontalAlignment','Left');
axis tight; xlabel('Time Iterations');
ylim([T_min T_max]); ylabel('Temperature');
legend('Red Point','Green Point ','Blue Point ','Location','northwest')
drawnow
hold off

%%%             ------------ Animated plot ----------

for j=1:n_time                               
    
subplot(2,2,1)
surf(x,y,T(:,:,j))
hold on
title(sprintf('Temperature at time : %i seconds ',round(j*((dt1+dt2)/2))))
plot3(  L/4,  H/4,T_max,'ko','markerfacecolor','r') % plot red point
plot3(  L/2,  H/2,T_max,'ko','markerfacecolor','g') % plot green point
plot3(3*L/4,3*H/4,T_max,'ko','markerfacecolor','b') % plot blue point
plot3(  L/4,  H/4,T_min,'ko','markerfacecolor','r') % plot red point
plot3(  L/2,  H/2,T_min,'ko','markerfacecolor','g') % plot green point
plot3(3*L/4,3*H/4,T_min,'ko','markerfacecolor','b') % plot blue point
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
title(sprintf('Temperature at time : %i seconds ',round(j*((dt1+dt2)/2))))
scatter(j,T(floor((NodesL+2)/4),floor((NodesH+2)/4),j),'r.'); 
scatter(j,T(ceil((NodesL+2)/2),ceil((NodesH+2)/2),j),'g.'); 
scatter(j,T(ceil(3*(NodesL+2)/4),ceil(3*(NodesH+2)/4),j),'b.'); 
axis tight; xlabel('Time Iterations');
axis tight; ylabel('Temperature');
legend('Red Point','Green Point ','Blue Point ','Location','northwest')
drawnow
hold off

end

