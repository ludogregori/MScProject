%% CONNECT PROJECT mission
% MSc Thesis
% 
% 
% Project: CONNECT
% Author: Ludovico Gregori
% URN : 6778145
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all; close all; clc; format lonG; %avoiding old data, variables, figures, ...

%Set Defaults
set(0,'DefaultTextInterpreter','latex');
set(0,'DefaultLegendInterpreter','latex');
set(0,'DefaultAxesFontSize',20);
set(0,'DefaultLineLineWidth',2);

% Plotting utility constants
global control_plot gif_plot save_plot
control_plot = 1; %to control each plot.
gif_plot = 0; %to control the gif plot.
save_plot = 0; %in order to save each plot, switch this parameter to 0 when you are trying to run it.

% Conversion factors
d2r = pi/180;   % from degrees to radians (deg^-1);
r2d = 180/pi;   % from radians to degrees (deg);

% Planet Constants Struct
planet.GM = 398600.4418;   % Earth's gravitational parameter (km3/s2)
planet.Re = 6378.137;      % Earth's equatorial radius (km)
planet.FaceAlpha = 0.50;   %Earth opacity

GM = planet.GM;
Re = planet.Re;
FaceAlpha = planet.FaceAlpha;

%% Orbit Parameters
altitude = 500; %km 
a = Re + altitude; %km
ecc     = 0.0;               %eccentricity 
I       = 15 * d2r;           %inclination angle
RAAN    = 0 * d2r;           %right ascension of ascending node
argP    = 0 * d2r;           %argument of perigee
M0      = 0 * d2r;           %Mean initial anomaly

fprintf(1, '\n');
disp('CONNECT MSc Thesis, Ludovico Gregori URN: 6778145')

n = sqrt(GM/a^3);   % mean  motion (rad/s)
P = 2*pi/n;         % Orbit period (s)

fprintf(1, '\n');
disp('Clohessy Wiltshire Orbit Integration & Data filtering')
fprintf(1, '\n');
fprintf('Mean motion of the satellite = %.4f rad/s.\n',n);
fprintf(1, '\n');
fprintf('Virtual Chief Orbital Period  = %.4f sec = %.4f hours.\n',P,P/3600);

%% ORBIT INTEGRATION

%% INPUTS 
%first set of inputs to perform the integration of the HCW equation to
%describe the chief/chasers system. 
Rform = 100; %[metres], formation radius
nc = 0.0011; %[rad/s], mean motion, nc = n of previous calcolations.

%second set of inputs to complete the design. 
prompt = "Type 1 (3 sats), 2 (4 sats) or 3 (5 sats):";
configuration = input(prompt);
if configuration ~= 1 && configuration ~= 2 && configuration ~= 3
    configuration = 1; %3 sats
elseif isempty(configuration)
    configuration = 1; %3 sats
end

%SNC = 0; %to use or not the state noise compensation, as a input in the optimised code
prompt = "Do you want to use the State Noise Compensation? [Y/N]:";
SNC = input(prompt,"s");
if SNC == 'Y'
    SNC = 1;
elseif SNC == 'N'
    SNC = 0;
elseif isempty(SNC)
    SNC = 0; %no SNC active
elseif SNC ~= 'Y' && SNC ~= 'N'
    SNC = 0;
end

%ToF = 24; %to select how long we want to integrate the system. 3 options: 1.5 hours, 24 hours, 10 days, as a input in the optimised code
prompt = "How long do you wanna propagate the mission for? [1 for 1 orbital period, 24 for 24 hours, 10 for 10 days]:";
ToF = input(prompt);

if  ToF ~= 1 && ToF ~= 24 && ToF ~= 10
    ToF = 1; %no SNC activeisempty(ToF)
elseif isempty(ToF)
    ToF = 24; %no SNC active
end

prompt = "Which type of orbit you want to select? [circular or oblate]: ";
orbit_type = input(prompt,'s'); %string
if strcmpi(orbit_type,'circular') == 0 && strcmpi(orbit_type,'oblate') == 0 %both false 
    orbit_type = 'circular';
elseif isempty(orbit_type)
    orbit_type = 'circular';
end

%% Possible Mission Configurations
if configuration == 1
    Nsat = 3;
    Ncentral = 0;
elseif configuration == 2
    Nsat = 4;
    Ncentral = 0;
elseif configuration == 3
    Nsat = 4;
    Ncentral = 1;
elseif configuration ~= 1 && configuration ~= 2 && configuration ~= 3
    Nsat = 3;
    Ncentral = 0;
end

%struct pars to move parameters faster
pars.nc = nc; %mean motion
Pc = 2*pi/nc; %[sec]
scale = 86400/Pc; %to obtain 24 hours when it is multiplied by Pc

if ToF == 1 
    orbital_period = 1; 
elseif ToF == 24
    orbital_period = scale;
elseif ToF == 10
    orbital_period = 10*scale;
end

t0 = 0; %initial time
tend = orbital_period*Pc; %final time
tvec = linspace(t0, tend + t0, 1000); %time vector of 1000 equally spaced points
tvec_design = linspace(t0, Pc + t0, 1000); %vector time to print the element along one orbit. 
theta0 = linspace(0,2*pi,Nsat+1); %initial phase difference

%equations via HCW hypotesis
%HCW equations (open loop, no control input vector u)
% xdotdot - 2nydot = 3n^2 *x
% ydotdot + 2nxdot = 0
% zdotdot          = -n^2 *z

%n = 'circular';
IC = zeros(6,Nsat); %intial conditions matrix with positions and velocities at t0
switch orbit_type 
    case 'circular'
        % Circular Orbit Integration
        [IC] = circularorbit(Nsat,Rform,nc,tvec,theta0);
    case 'oblate'
        % Oblate Orbit Integration
        [IC] = oblateorbit(Nsat,Rform,nc,tvec,theta0);
end

% %sanity check, need to finish to code this.
% X1 = Rform* [sqrt(3)/2 0 1.2 nc/2 -sqrt(3)*nc -sqrt(3)*nc/2]'; %first deputy sat
% X2 = 1.0e+02 *[0.002125447556331  -2.999995482467587  -1.348708466560883  -0.001099997515357  -0.000004675984624  -0.000669994344300]'; % second deputy sat,position 1925 of Xhcw1
% X3 = [-87.288001945011317  -2.413018973902471  12.955673201642741   0.053672839564357   0.192033604279026   0.162160112677530]'; % third deputy sat,position 3850 of Xhcw1
% Xinit = [IC(:,1); IC(:,2); IC(:,3)]; 

%% CHASERS AROUND THE VIRTUAL CHIEF IN RTN REF FRAME
%satelliti follower rispetto al virtual chief
options = odeset('RelTol',2.22045e-14,'AbsTol',1e-16); %ode options for tolerances 

[T,Xhcw] = ode45(@HCWint,tvec,IC,options,pars,Nsat); %single output vector 18x[tvec] or 24x[tvec].
Xhcw = Xhcw'; %better way to visualise it. 

% %% rotating the coordinates to ECI/ECEF
% %for sat 1, then switch to all 3 sats
% r = Xhcw(1:3,1)/vecnorm(Xhcw(1:3,1),2); %radial axis
% n = cross(Xhcw(4:6,1),Xhcw(1:3,1))/vecnorm(cross(Xhcw(4:6,1),Xhcw(1:3,1)),2); %v is positive clockwise, n should point upward, cross-track
% t = cross(n,r)/vecnorm(cross(n,r),2); %along-track
%         % 
%         % figure(2)
%         % hold on;grid on;
%         % %RTN
%         % quiver3(0, 0, 0, r(1,1), r(2,1), r(3,1), 1e4, 'b', 'LineWidth', 2);
%         % quiver3(0, 0, 0, t(1,1), t(2,1), t(3,1), 1e4, 'b', 'LineWidth', 2);
%         % quiver3(0, 0, 0, n(1,1), n(2,1), n(3,1), 1e4, 'b', 'LineWidth', 2);
%         % text(1e4*r(1,1), 1e4*r(2,1), 1e4*r(3,1), '$$\hat{r}_{rtn}$$', 'FontSize', 20, 'Interpreter', 'Latex', 'Color', 'b')
%         % text(1e4*t(1,1), 1e4*t(2,1), 1e4*t(3,1), '$$\hat{t}_{rtn}$$', 'FontSize', 20, 'Interpreter', 'Latex', 'Color', 'b')
%         % text(1e4*n(1,1), 1e4*n(2,1), 1e4*n(3,1), '$$\hat{n}_{rtn}$$', 'FontSize', 20, 'Interpreter', 'LAtex', 'Color', 'b')
%         % %ECI
%         % quiver3(0, 0, 0, 1, 0, 0, 2e4, '--r', 'LineWidth', 1);
%         % quiver3(0, 0, 0, 0, 1, 0, 2e4, '--r', 'LineWidth', 1);
%         % quiver3(0, 0, 0, 0, 0, 1, 2e4, '--r', 'LineWidth', 1);
%         % text(2e4, 0, 0, '$$\hat{I}_{ECI}$$', 'FontSize', 20, 'Interpreter', 'Latex', 'Color', 'b')
%         % text(0, 2e4, 0, '$$\hat{J}_{ECI}$$', 'FontSize', 20, 'Interpreter', 'Latex', 'Color', 'b')
%         % text(0, 0, 2e4, '$$\hat{K}_{ECI}$$', 'FontSize', 20, 'Interpreter', 'LAtex', 'Color', 'b')
%         % %ECEF
%         % %DA AGGIUNGERE
%         % view(-45,+5);
%         % axis equal;
%         % hold off;
%         % % da vedere bene, prendere esempio ecef to eci da esercizio
%         % % Rrtn2eci = [dot(r(1),(1 0 0)) dot(,(1 0 0)) dot(,(1 0 0));
%         % %             dot(r(2),(0 1 0)) dot(,(0 1 0)) dot(,(0 1 0));
%         % %             dot(r(3),(0 0 1)) dot(,(0 0 1)) dot(,(0 0 1))];

%[T,Xhcw1] = ode45(@HCWint,tvec,X1,options,pars); %first dep ode, fare tutto con un ode sola, mettere un vettore 18x1. 
%[T,Xhcw2] = ode45(@HCWint,tvec,X2,options,pars); %second dep ode
%[T,Xhcw3] = ode45(@HCWint,tvec,X3,options,pars); %third dep ode

%% PLOTTING
if control_plot == 1
    figure('name','HCW plot','units','normalized','outerposition',[0.2 0.2 0.8 0.8]); %figure 1
    view(-40,+5); %POV for the plot
    %trajectories
    hold on; grid on;
    switch configuration
        case 1 %3 sats
            %formation 1
            %[output] = plotconf1(Xhcw);
            origin = plot3(0,0,0,'kd'); %virtual chief placed in the origin.
            sat1 = plot3(Xhcw(1,:),Xhcw(2,:),Xhcw(3,:),'r-','LineWidth',1); %trajectory of sat1, if not using a gif.
            sat2 = plot3(Xhcw(7,:),Xhcw(8,:),Xhcw(9,:),'b-','LineWidth',1); %trajectory of sat2, if not using a gif.
            sat3 = plot3(Xhcw(13,:),Xhcw(14,:),Xhcw(15,:),'g-','LineWidth',1); %trajectory of sat3, if not using a gif.

            %initial and final positions
            initial_pos1 = plot3(Xhcw(1,1), Xhcw(2,1), Xhcw(3,1), 'sr', 'MarkerFaceColor', 'r','MarkerSize',10); %deputy initial position
            final_pos1   = plot3(Xhcw(1,end), Xhcw(2,end), Xhcw(3,end), 'sr', 'MarkerFaceColor', 'r'); %formation final position
            ite_sat1 = scatter3(Xhcw(3,1), Xhcw(2,1), Xhcw(3,1),[20],'or','filled','MarkerEdgeColor','r');  %first iteration for our formation

            initial_pos2 = plot3(Xhcw(7,1), Xhcw(8,1), Xhcw(9,1), 'sb', 'MarkerFaceColor', 'b','MarkerSize',10); %deputy initial position
            final_pos2   = plot3(Xhcw(7,end), Xhcw(8,end), Xhcw(9,end), 'sb', 'MarkerFaceColor', 'b'); %formation final position
            ite_sat2 = scatter3(Xhcw(7,1), Xhcw(8,1), Xhcw(9,1),[20],'ob','filled','MarkerEdgeColor','b');  %first iteration for our formationinitial_pos = plot3(Xhcw(1,1), Xhcw(7,1), Xhcw(13,1), 'sr', 'MarkerFaceColor', 'b','MarkerSize',10); %deputy initial position

            initial_pos3 = plot3(Xhcw(13,1), Xhcw(14,1), Xhcw(15,1), 'sg', 'MarkerFaceColor', 'g','MarkerSize',10); %deputy initial position
            final_pos3   = plot3(Xhcw(13,end), Xhcw(14,end), Xhcw(15,end), 'sg', 'MarkerFaceColor', 'g'); %formation final position
            ite_sat3 = scatter3(Xhcw(13,1), Xhcw(14,1), Xhcw(15,1),[20],'og','filled','MarkerEdgeColor','g');  %first iteration for our formation

            output = [origin sat1 sat2 sat3 initial_pos1 initial_pos2 initial_pos3];
            %output = [origin sat1 sat2 sat3 initial_pos1 final_pos1 ite_sat1 initial_pos2 final_pos2 ite_sat2 initial_pos3 final_pos3 ite_sat3];

        case 2 %4 sats
            %formation 2
            %[output] = plotconf2(Xhcw);
            origin = plot3(0,0,0,'kd'); %virtual chief placed in the origin.
            sat1 = plot3(Xhcw(1,:),Xhcw(2,:),Xhcw(3,:),'b-','LineWidth',1); %trajectory of sat1, if not using a gif.
            sat2 = plot3(Xhcw(7,:),Xhcw(8,:),Xhcw(9,:),'m-','LineWidth',1); %trajectory of sat2, if not using a gif.
            sat3 = plot3(Xhcw(13,:),Xhcw(14,:),Xhcw(15,:),'r-','LineWidth',1); %trajectory of sat3, if not using a gif.
            sat4 = plot3(Xhcw(19,:),Xhcw(20,:),Xhcw(21,:),'b-','LineWidth',1); %trajectory of sat1, if not using a gif.

            %initial and final positions
            initial_pos1 = plot3(Xhcw(1,1), Xhcw(2,1), Xhcw(3,1), 'sr', 'MarkerFaceColor', 'k','MarkerSize',10); %deputy initial position
            final_pos1   = plot3(Xhcw(1,end), Xhcw(2,end), Xhcw(3,end), 'sr', 'MarkerFaceColor', 'r'); %formation final position
            ite_sat1 = scatter3(Xhcw(3,1), Xhcw(2,1), Xhcw(3,1),[20],'or','filled','MarkerEdgeColor','k');  %first iteration for our formation

            initial_pos2 = plot3(Xhcw(7,1), Xhcw(8,1), Xhcw(9,1), 'sb', 'MarkerFaceColor', 'k','MarkerSize',10); %deputy initial position
            final_pos2   = plot3(Xhcw(7,end), Xhcw(8,end), Xhcw(9,end), 'sb', 'MarkerFaceColor', 'r'); %formation final position
            ite_sat2 = scatter3(Xhcw(7,1), Xhcw(8,1), Xhcw(9,1),[20],'ob','filled','MarkerEdgeColor','k');  %first iteration for our formationinitial_pos = plot3(Xhcw(1,1), Xhcw(7,1), Xhcw(13,1), 'sr', 'MarkerFaceColor', 'b','MarkerSize',10); %deputy initial position

            initial_pos3 = plot3(Xhcw(13,1), Xhcw(14,1), Xhcw(15,1), 'sg', 'MarkerFaceColor', 'k','MarkerSize',10); %deputy initial position
            final_pos3   = plot3(Xhcw(13,end), Xhcw(14,end), Xhcw(15,end), 'sg', 'MarkerFaceColor', 'r'); %formation final position
            ite_sat3 = scatter3(Xhcw(13,1), Xhcw(14,1), Xhcw(15,1),[20],'og','filled','MarkerEdgeColor','k');  %first iteration for our formation

            initial_pos4 = plot3(Xhcw(19,1), Xhcw(20,1), Xhcw(21,1), 'sm', 'MarkerFaceColor', 'k','MarkerSize',10); %deputy initial position
            final_pos4  = plot3(Xhcw(19,end), Xhcw(20,end), Xhcw(21,end), 'sm', 'MarkerFaceColor', 'r'); %formation final position
            ite_sat4 = scatter3(Xhcw(19,1), Xhcw(20,1), Xhcw(21,1),[20],'om','filled','MarkerEdgeColor','k');  %first iteration for our formation

            output = [origin sat1 sat2 sat3 sat4 initial_pos1 final_pos1 ite_sat1 initial_pos2 final_pos2 ite_sat2 initial_pos3 final_pos3 ite_sat3 initial_pos4 final_pos4 ite_sat4];

        case 3 %4sats + 1central
            %formation 3
            %[output] = plotconf3(Xhcw);
            sat0 = plot3(0,0,0,'sk', 'MarkerFaceColor', 'k','MarkerSize',10);  %virtual chief placed in the origin.
            sat1 = plot3(Xhcw(1,:),Xhcw(2,:),Xhcw(3,:),'b-','LineWidth',1);    %trajectory of sat1, if not using a gif.
            sat2 = plot3(Xhcw(7,:),Xhcw(8,:),Xhcw(9,:),'m-','LineWidth',1);    %trajectory of sat2, if not using a gif.
            sat3 = plot3(Xhcw(13,:),Xhcw(14,:),Xhcw(15,:),'r-','LineWidth',1); %trajectory of sat3, if not using a gif.
            sat4 = plot3(Xhcw(19,:),Xhcw(20,:),Xhcw(21,:),'b-','LineWidth',1); %trajectory of sat1, if not using a gif.

            %initial and final positions
            initial_pos1 = plot3(Xhcw(1,1), Xhcw(2,1), Xhcw(3,1), 'sr', 'MarkerFaceColor', 'k','MarkerSize',10); %deputy initial position
            final_pos1   = plot3(Xhcw(1,end), Xhcw(2,end), Xhcw(3,end), 'sr', 'MarkerFaceColor', 'r'); %formation final position
            ite_sat1 = scatter3(Xhcw(3,1), Xhcw(2,1), Xhcw(3,1),[20],'or','filled','MarkerEdgeColor','k');  %first iteration for our formation

            initial_pos2 = plot3(Xhcw(7,1), Xhcw(8,1), Xhcw(9,1), 'sb', 'MarkerFaceColor', 'k','MarkerSize',10); %deputy initial position
            final_pos2   = plot3(Xhcw(7,end), Xhcw(8,end), Xhcw(9,end), 'sb', 'MarkerFaceColor', 'r'); %formation final position
            ite_sat2 = scatter3(Xhcw(7,1), Xhcw(8,1), Xhcw(9,1),[20],'ob','filled','MarkerEdgeColor','k');  %first iteration for our formationinitial_pos = plot3(Xhcw(1,1), Xhcw(7,1), Xhcw(13,1), 'sr', 'MarkerFaceColor', 'b','MarkerSize',10); %deputy initial position

            initial_pos3 = plot3(Xhcw(13,1), Xhcw(14,1), Xhcw(15,1), 'sg', 'MarkerFaceColor', 'k','MarkerSize',10); %deputy initial position
            final_pos3   = plot3(Xhcw(13,end), Xhcw(14,end), Xhcw(15,end), 'sg', 'MarkerFaceColor', 'r'); %formation final position
            ite_sat3 = scatter3(Xhcw(13,1), Xhcw(14,1), Xhcw(15,1),[20],'og','filled','MarkerEdgeColor','k');  %first iteration for our formation

            initial_pos4 = plot3(Xhcw(19,1), Xhcw(20,1), Xhcw(21,1), 'sm', 'MarkerFaceColor', 'k','MarkerSize',10); %deputy initial position
            final_pos4  = plot3(Xhcw(19,end), Xhcw(20,end), Xhcw(21,end), 'sm', 'MarkerFaceColor', 'r'); %formation final position
            ite_sat4 = scatter3(Xhcw(19,1), Xhcw(20,1), Xhcw(21,1),[20],'om','filled','MarkerEdgeColor','k');  %first iteration for our formation

            output = [sat1 sat2 sat3 sat4 sat0 initial_pos1 final_pos1 ite_sat1 initial_pos2 final_pos2 ite_sat2 initial_pos3 final_pos3 ite_sat3 initial_pos4 final_pos4 ite_sat4];
    end

    % %projections on the plot
    % xcoordinate = repelem(100,length(Xhcw(1,:))); %1x1000, placing the projection at the x coordinate. 
    % yzproj = plot3(xcoordinate(1,:),Xhcw(2,:),Xhcw(3,:),'b--','LineWidth',1); %plotting along the crosstrack-alongtrack plane
    % zcoordinate = repelem(130,length(Xhcw(3,:))); %1x1000, placing the projection at the z coordinate.
    % yxptroj = plot3(Xhcw(1,:),Xhcw(2,:),zcoordinate(1,:),'b--','LineWidth',1); %plotting along the radial-alongtrack plane
    
    %% ANIMATED GIF
    if gif_plot == 1
        %%Iterating through the length of the time array
        %f = factor(round(tvec(end)));
            switch configuration
                case 1 %3 sats
                    for k = 1:6:length(tvec) %end of the time vector used for the integration of the deputies' trajectories.
                        % Updating the line
                        sat1.XData = Xhcw(1,1:k);
                        sat1.YData = Xhcw(2,1:k);
                        sat1.ZData = Xhcw(3,1:k);
                        sat2.XData = Xhcw(7,1:k);
                        sat2.YData = Xhcw(8,1:k);
                        sat2.ZData = Xhcw(9,1:k);
                        sat3.XData = Xhcw(13,1:k);
                        sat3.YData = Xhcw(14,1:k);
                        sat3.ZData = Xhcw(15,1:k);
    
                        % Updating the point
                        ite_sat1.XData = Xhcw(1,k);
                        ite_sat1.YData = Xhcw(2,k);
                        ite_sat1.ZData = Xhcw(3,k);
                        ite_sat2.XData = Xhcw(7,k);
                        ite_sat2.YData = Xhcw(8,k);
                        ite_sat2.ZData = Xhcw(9,k);
                        ite_sat3.XData = Xhcw(13,k);
                        ite_sat3.YData = Xhcw(14,k);
                        ite_sat3.ZData = Xhcw(15,k);
    
                        %Updating the title
                        title(sprintf('Trajectory of Deputy Satellites\n Time: %0.2f sec', tvec(k)),...
                            'Interpreter','Latex');
                        %view(-45,+20);
                        drawnow; %forcing to draw each instant.
    
                    end
    
                case 2 %4 sats
                    for k = 1:6:length(tvec) %end of the time vector used for the integration of the deputies' trajectories.
                        % Updating the line
                        sat1.XData = Xhcw(1,1:k);
                        sat1.YData = Xhcw(2,1:k);
                        sat1.ZData = Xhcw(3,1:k);
                        sat2.XData = Xhcw(7,1:k);
                        sat2.YData = Xhcw(8,1:k);
                        sat2.ZData = Xhcw(9,1:k);
                        sat3.XData = Xhcw(13,1:k);
                        sat3.YData = Xhcw(14,1:k);
                        sat3.ZData = Xhcw(15,1:k);
                        sat4.XData = Xhcw(19,1:k);
                        sat4.YData = Xhcw(20,1:k);
                        sat4.ZData = Xhcw(21,1:k);
                        % Updating the point
                        ite_sat1.XData = Xhcw(1,k);
                        ite_sat1.YData = Xhcw(2,k);
                        ite_sat1.ZData = Xhcw(3,k);
                        ite_sat2.XData = Xhcw(7,k);
                        ite_sat2.YData = Xhcw(8,k);
                        ite_sat2.ZData = Xhcw(9,k);
                        ite_sat3.XData = Xhcw(13,k);
                        ite_sat3.YData = Xhcw(14,k);
                        ite_sat3.ZData = Xhcw(15,k);
                        ite_sat4.XData = Xhcw(19,k);
                        ite_sat4.YData = Xhcw(20,k);
                        ite_sat4.ZData = Xhcw(21,k);
    
                        %Updating the title
                        title(sprintf('Trajectory of Deputy Satellites\n Time: %0.2f sec', tvec(k)),...
                            'Interpreter','Latex');
                        %view(-45,+20);
                        drawnow; %forcing to draw each instant.
                        %pause(100) %speed of the position update
                    end
    
                    %Delay
                    % pause(100) %speed of the position update
    
                    % % Saving the figure %provare senza pause e con questo on per
                    %         %save as a gif. 
                    %         frame = getframe(gcf);
                    %         im = frame2im(frame);
                    %         [imind,cm] = rgb2ind(im,256);
                    %         if k == 1
                    %             imwrite(imind,cm,filename,'gif', 'Loopcount',inf,...
                    %             'DelayTime',0.1);
                    %         else
                    %             imwrite(imind,cm,filename,'gif','WriteMode','append',...
                    %             'DelayTime',0.1);
                    %         end
                    %     end

            end %closing the switch
end %closing the gif_plot

%projections on the plot
xcoordinate = repelem(100,length(Xhcw(1,:))); %1x1000, placing the projection at the x coordinate.
yzproj = plot3(xcoordinate(1,:),Xhcw(2,:),Xhcw(3,:),'b--','LineWidth',1); %plotting along the crosstrack-alongtrack plane
zcoordinate = repelem(130,length(Xhcw(3,:))); %1x1000, placing the projection at the z coordinate.
yxptroj = plot3(Xhcw(1,:),Xhcw(2,:),zcoordinate(1,:),'b--','LineWidth',1); %plotting along the radial-alongtrack plane

xlabel('Radial [m]','FontSize',16);
ylabel('Along-Track [m]','FontSize',16);
zlabel('Cross-Track [m]','FontSize',16);
axis([-100 +100 -100 +100]);

savingplot(configuration,orbit_type,ToF,output); %saving each figure 1 based on the configuration and adding the title.

hold off;
  
end %closing the control_plot

%% array design antennas positions along the orbit 
%adding the third coordinate as the radial component. 
%it is changing every "elements number" iterations, cus the array is in a
%single plane.

%% AGGIUNGERE SWITCH PER PRENDERE COORDINATE DI ORBITA CIRCOLARE E ORBITA OBLATA
switch configuration
    case 1 %3 sats
        positions = [Xhcw(1,:)' Xhcw(2,:)' Xhcw(3,:)' Xhcw(7,:)' Xhcw(8,:)' Xhcw(9,:)' Xhcw(13,:)' Xhcw(14,:)' Xhcw(15,:)']; %1000x9
        %save('/Users/ludovicogregori/Desktop/SURREY/FINAL PROJECT/saving data/positions.mat','positions') %saving the sat positions 

        %% relative position
        % distance between Xhcw1-Xhcw2, Xhcw2-Xhcw3, Xhcw1-Xhcw-3
        vec12 = [positions(:,1)-positions(:,4) positions(:,2) - positions(:,5) positions(:,3)-positions(:,6) ]; %vector between sat1 and sat2
        vec13 = [positions(:,1)-positions(:,7) positions(:,2) - positions(:,8) positions(:,3)-positions(:,9) ]; %vector between sat1 and sat3
        vec23 = [positions(:,4)-positions(:,7) positions(:,5) - positions(:,8) positions(:,6)-positions(:,9) ]; %vector between sat2 and sat3
        distvec12 = vecnorm(vec12,2,2); %distance vector between sat1 and sat2
        distvec13 = vecnorm(vec13,2,2); %distance vector between sat1 and sat3
        distvec23 = vecnorm(vec23,2,2); %distance vector between sat2 and sat3

        distmatrix = [distvec12 distvec13 distvec23]; %matrix with all the distances as columns
        errormatrix = abs(diff(distmatrix,1)); %first-order difference

    case 2 %4 sats
        positions = [Xhcw(1,:)' Xhcw(2,:)' Xhcw(3,:)' Xhcw(7,:)' Xhcw(8,:)' Xhcw(9,:)' Xhcw(13,:)' Xhcw(14,:)' Xhcw(15,:)' Xhcw(19,:)' Xhcw(20,:)' Xhcw(21,:)']; %1000x9
        %save('/Users/ludovicogregori/Desktop/SURREY/FINAL PROJECT/saving data/positions.mat','positions') %saving the sat positions
        %distance vector along the orbit between sat1 and sat2

        %% relative position
        % distance between Xhcw1-Xhcw2, Xhcw2-Xhcw3,
        % Xhcw3-Xhcw-4,Xhcw4-Xhcw1,Xhcw1-Xhcw3,Xhcw2-Xhc4
        vec12 = [positions(:,1) - positions(:,4)  positions(:,2) -  positions(:,5)  positions(:,3) - positions(:,6) ];   %vector between sat1 and sat2
        vec13 = [positions(:,1) - positions(:,7)  positions(:,2) -  positions(:,8)  positions(:,3) - positions(:,9) ];   %vector between sat1 and sat3
        vec23 = [positions(:,4) - positions(:,7)  positions(:,5) -  positions(:,8)  positions(:,6) - positions(:,9) ];   %vector between sat2 and sat3
        vec24 = [positions(:,4) - positions(:,10) positions(:,5) -  positions(:,11) positions(:,6) - positions(:,12)];   %vector between sat2 and sat4
        vec34 = [positions(:,7) - positions(:,10) positions(:,8) -  positions(:,11) positions(:,9) - positions(:,12)];   %vector between sat3 and sat4
        vec41 = [positions(:,10) - positions(:,1) positions(:,11) - positions(:,2)  positions(:,12)- positions(:,3) ];   %vector between sat4 and sat1
        
        distvec12 = vecnorm(vec12,2,2); %distance vector between sat1 and sat2
        distvec23 = vecnorm(vec23,2,2); %distance vector between sat2 and sat3
        distvec34 = vecnorm(vec34,2,2); %distance vector between sat3 and sat4
        distvec41 = vecnorm(vec41,2,2); %distance vector between sat4 and sat1
        distvec24 = vecnorm(vec24,2,2); %distance vector between sat2 and sat4
        distvec13 = vecnorm(vec13,2,2); %distance vector between sat1 and sat3
        
        distmatrix = [distvec12 distvec23 distvec34 distvec41 distvec13 distvec24]; %matrix with all the distances as columns
        errormatrix = abs(diff(distmatrix,1)); %first-order difference

    case 3 %4 sats + 1sat central
        
        positions = [Xhcw(1,:)' Xhcw(2,:)' Xhcw(3,:)' Xhcw(7,:)' Xhcw(8,:)' Xhcw(9,:)' Xhcw(13,:)' Xhcw(14,:)' Xhcw(15,:)' Xhcw(19,:)' Xhcw(20,:)' Xhcw(21,:)']; %1000x9
        % distance between Xhcw1-Xhcw0, Xhcw2-Xhcw0, Xhcw3-Xhcw0, Xhcw4-Xhcw0
        vec10 = [positions(:,1) positions(:,2) positions(:,3)];
        vec20 = [positions(:,4) positions(:,5) positions(:,6)];
        vec30 = [positions(:,7) positions(:,8) positions(:,9)];
        vec40 = [positions(:,10) positions(:,11) positions(:,12)];
        distvec10 = vecnorm(vec10,2,2); %distance vector between sat1 and sat0
        distvec20 = vecnorm(vec20,2,2); %distance vector between sat2 and sat0
        distvec30 = vecnorm(vec30,2,2); %distance vector between sat3 and sat0
        distvec40 = vecnorm(vec40,2,2); %distance vector between sat1 and sat0
        
        distmatrix = [distvec10 distvec20 distvec30 distvec40]; %matrix with all the distances as columns
        errormatrix = abs(diff(distmatrix,1)); %first-order difference
end

[Centres3Dtot_reshaped] = arraydesign(positions,tvec_design); %57700x9
%save('/Users/ludovicogregori/Desktop/SURREY/FINAL PROJECT/saving data/positions.mat','positions') %saving the sat positions for Anouk & Callum
%output from the function the centres of each system
%if we need the coordinates of each smaller antennas we can output coor3Dtot. 

%% relative position 
% distance between Xhcw1-Xhcw2, Xhcw2-Xhcw3, Xhcw1-Xhcw-3
% magnitude to obtain the exact metres, try to fix that distance at 100m,
% modify the radius formation parameter to obtain a fixed 100m distance
% between deputy satellites.
            
            % %distance vector along the orbit between sat1 and sat2
            % vec12 = [positions(:,1)-positions(:,4) positions(:,2) - positions(:,5) positions(:,3)-positions(:,6) ]; %vector between sat1 and sat2
            % vec13 = [positions(:,1)-positions(:,7) positions(:,2) - positions(:,8) positions(:,3)-positions(:,9) ]; %vector between sat1 and sat3
            % vec23 = [positions(:,4)-positions(:,7) positions(:,5) - positions(:,8) positions(:,6)-positions(:,9) ]; %vector between sat2 and sat3
            % distvec12 = vecnorm(vec12,2,2); %distance vector between sat1 and sat2
            % distvec13 = vecnorm(vec13,2,2); %distance vector between sat1 and sat3
            % distvec23 = vecnorm(vec23,2,2); %distance vector between sat2 and sat3
            % distmatrix = [distvec12 distvec13 distvec23]; %matrix with all the distances as columns
            % errormatrix = abs(diff(distmatrix,1)); %first-order difference
            % %checking the max relative error in the entire constallation for the centres of the formations.
            % errors1 = errormatrix(errorposition1);
            % maxerror1 = max(errormatrix);
            
            % %distance between each elements in the different satellite.
            % elemvec12 = [Centres3Dtot_reshaped(:,1) - Centres3Dtot_reshaped(:,4) Centres3Dtot_reshaped(:,2)-Centres3Dtot_reshaped(:,5) Centres3Dtot_reshaped(:,3)-Centres3Dtot_reshaped(:,6)]; %vector between
            % elemvec13 = [Centres3Dtot_reshaped(:,1) - Centres3Dtot_reshaped(:,7) Centres3Dtot_reshaped(:,2)-Centres3Dtot_reshaped(:,8) Centres3Dtot_reshaped(:,3)-Centres3Dtot_reshaped(:,9)]; %vector between
            % elemvec23 = [Centres3Dtot_reshaped(:,4) - Centres3Dtot_reshaped(:,7) Centres3Dtot_reshaped(:,5)-Centres3Dtot_reshaped(:,8) Centres3Dtot_reshaped(:,6)-Centres3Dtot_reshaped(:,9)]; %vector between
            % elemdist12 = vecnorm(elemvec12,2,2); %distance vector between
            % elemdist13 = vecnorm(elemvec13,2,2); %distance vector between
            % elemdist23 = vecnorm(elemvec23,2,2); %distance vector between
            % elemdistmatrix = [elemdist12 elemdist13 elemdist23]; %matrix with all the distances as columns
            % elemerrormatrix = abs(diff(elemdistmatrix,1)); %first-order difference
            % %checking the max relative error in the entire constallation for each elements.
            % errors2 = elemerrormatrix(errorposition2);
            % maxerror2 = max(elemerrormatrix);

%% DATA FILTERING (Kalman) 
% get the relative distance and relative velocities in it to define the
% matrices and the new state vector with a control error. 
% % sub-centimetre and sub-centimetre/sec

%% nominal parameters
sigmaDGPS    = 0.10;    %[m], 10 cm
sigmaISL     = 10^-9;  %[m], nanometer
sigmaDoppler = 10; %[m/s]

%% DGPS, ISLS & DGPS + ISLS - combined cases 
% we are using  the coordinate of the position of each satellite plus the
% relative distance. two sensors togheter to hopefully improve the accuracy
% state vector 18x1, C matrix 12x18, Ycomp = 12x1
switch configuration
    case 1 %3 sats
        n = 18;  % no. of state variables, (6 variables per each satellite, 3 satellites in total),(x,y,z,xdot,ydot,zdot)
        I = eye(n);
        m = length(tvec);     % no. of observations
        %% nominal parameters
        sigmaDGPS = 0.10; %[m]
        sigmaISL = 10^-9; %[m]

        %% least squares solution
        prePk_1bar = diag([1e4,1e4,1e4,1,1,1]);   % a-priori covariance matrix, sigma^2
        Pk_1bar = blkdiag(prePk_1bar,prePk_1bar,prePk_1bar); %18x18

        %% Integrate Reference Trajectory & state transition matrix
        %STM = zeros(18,18); %initialising the STM matrix
        cellSTM = {};
        %Loop over measurements
        for z = 1:3
            for tt = 1:m
                %Retrieve STMs
                %for i = 1 : length(tvec)
                c = cos(nc*tvec(tt));
                s = sin(nc*tvec(tt));

                STM = [4-3*c             0    0      s/nc           (2/nc)*(1-c)      0;
                    6*(s-nc*tvec(tt)) 1    0 -2*(1-c)/nc   (4*s-3*nc*tvec(tt))/nc  0;
                    0                 0    c       0               0             s/nc;
                    3*nc*s            0    0       c              2*s              0;
                    -6*nc*(1-c)        0    0     -2*s            4*c-3             0;
                    0                 0  -nc*s     0               0               c]; %STM matrix

                newSTM(:,:,tt) = blkdiag(STM,STM,STM); %18x18
                cellSTM(tt) = {newSTM(:,:,tt)}; %storing the STM in a cell array

                if(tt == 1)
                    STMkk_1 = newSTM(:,:,tt);
                else
                    STMkk_1 = newSTM(:,:,tt) / newSTM(:,:,tt-1);  % STM(tk,tk_1) = STM(tk, 0) / STM(tk_1, 0), conventional KF
                end

                [Pkbar] = sncompensation(SNC,tt,STMkk_1,Pk_1bar,configuration); %State Noise Compensation

                %% MEASUREMENT UPDATE
                % state-observation relationship & C-matrix
                % (12x18) partial derivatvies wrt the state with the DGPS+ ISLs with the
                % relative positions rho12, rho13, rho23 (constant, so no need to evaluate over Xdot).
                Ck = zeros(12,18);
                rho12 = distvec12(tt,:); %norm of the dist to get the magn.
                rho13 = distvec13(tt,:);
                rho23 = distvec23(tt,:);

                Ck(1:3,1:3) = eye(3);   %x1,y1,z1, DGPS
                Ck(4:6,7:9) = eye(3);   %x2,y2,z2, DGPS
                Ck(7:9,13:15) = eye(3); %x3,y3,z3, DGPS
                Ck(10,1:3)   = (vec12(tt,:)/rho12); %per each iteration, ISLS
                Ck(10,7:9)   = - Ck(10,1:3); %ISLS
                Ck(11,1:3)   = (vec13(tt,:)/rho13); %per each iteration, ISLS
                Ck(11,13:15) = - Ck(11,1:3); %ISLS
                Ck(12,7:9)   = (vec23(tt,:)/rho23); %per each iteration, ISLS
                Ck(12,13:15) = -Ck(12,7:9); %ISLS

                CkDGPS = Ck(1:9,1:18);   %9x18
                CkISLS = Ck(10:12,1:18); %3x18
                ck = {CkDGPS CkISLS Ck}; %3 cases

                %measurement-noise covariance matrix
                d = [repelem(sigmaDGPS^2,9) repelem(sigmaISL^2,3)];
                Rk = diag(d); %initialising Rk, 12x12
                RkDGPS = Rk(1:9,1:9);
                RkISLS = Rk(10:12,10:12);
                rk = {RkDGPS RkISLS Rk}; %all the cases in a single cell array per each iteration
                Kk = (Pkbar*ck{z}')*((rk{z}+(ck{z}*Pkbar*ck{z}'))^-1); % (18x12) Kalman Gain

                % CKF specific to evaluate the standard deviations
                % a-posteriori covariance matrix
                Pk(:,:,tt)  = (I-(Kk*ck{z}))*Pkbar; %18x18, for the next iteration
                P12(:,:,tt) = Pk(1:6,1:6,tt)   + Pk(7:12,7:12,tt)   - Pk(1:6,7:12,tt)   - Pk(7:12,1:6,tt);   %selecting the components for the sigma^2
                P13(:,:,tt) = Pk(1:6,1:6,tt)   + Pk(13:18,13:18,tt) - Pk(1:6,13:18,tt)  - Pk(13:18,1:6,tt);  %selecting the components for the sigma^2
                P23(:,:,tt) = Pk(7:12,7:12,tt) + Pk(13:18,13:18,tt) - Pk(7:12,13:18,tt) - Pk(13:18,7:12,tt); %selecting the components for the sigma^2

                % store standard deviations
                sig12(:,tt) = diag(P12(:,:,tt)); %sigma^2
                sig13(:,tt) = diag(P13(:,:,tt)); %sigma^2
                sig23(:,tt) = diag(P23(:,:,tt)); %sigma^2
                RMS12 = 3 * sqrt(sig12(1,:) + sig12(2,:) + sig12(3,:)); %root mean square, 3-standard dev.
                RMS13 = 3 * sqrt(sig13(1,:) + sig13(2,:) + sig13(3,:)); %root mean square, 3-standard dev.
                RMS23 = 3 * sqrt(sig23(1,:) + sig23(2,:) + sig23(3,:)); %root mean square, 3-standard dev.

                RMS{z} = [RMS12' RMS13' RMS23']; %saving DGPS, ISLS , DGPS + ISLS as a cell array, each has 1000x3 elements

                Pk_1bar = Pk(:,:,tt); %18x18 , Update a-priori state deviation and covariance matrix before moving to next obs
                covcase1(tt) = {Pk(:,:,tt)}; %storing the covariance matrix in a cell array

            end % end of observation loop

            %% Plots
            %% Stand-Dev Check
            if control_plot == 1
                plotcase1(z,tvec,RMS12,RMS13,RMS23,orbit_type,ToF,SNC); %figure 2, figure 3, figure 4
            end
        end %closing z loop
        
        %RMS vector generated completely
        plotXALLconfig1(z,tvec,RMS,orbit_type,ToF,SNC) %plotting all systems per sat

        %single plot DGPS vs (DGPS & ISLs)
        if control_plot == 1
            plotcase1DGPSvsDGPSISLS(tvec,RMS,orbit_type,ToF,SNC) %figure 5
        end

    case 2 %4sats
        n = 24;  % no. of state variables, (6 variables per each satellite, 4 satellites in total),(x,y,z,xdot,ydot,zdot)
        I = eye(n);
        m = length(tvec);     % no. of observations
        %% nominal parameters
        sigmaDGPS = 0.10; %[m]
        sigmaISL = 10^-9; %[m]

        %% least squares solution
        prePk_1bar = diag([1e4,1e4,1e4,1,1,1]);   % a-priori covariance matrix, sigma^2
        Pk_1bar = blkdiag(prePk_1bar,prePk_1bar,prePk_1bar,prePk_1bar); %24x24

        %% Integrate Reference Trajectory & state transition matrix
        %STM = zeros(24,24); %initialising the STM matrix
        cellSTM = {};
        % Loop over measurements
        for z = 1:3
            for tt = 1:m

                %Retrieve STMs
                %for i = 1 : length(tvec)
                c = cos(nc*tvec(tt));
                s = sin(nc*tvec(tt));

                STM = [4-3*c             0    0      s/nc           (2/nc)*(1-c)      0;
                    6*(s-nc*tvec(tt)) 1    0 -2*(1-c)/nc   (4*s-3*nc*tvec(tt))/nc  0;
                    0                 0    c       0               0             s/nc;
                    3*nc*s            0    0       c              2*s              0;
                    -6*nc*(1-c)        0    0     -2*s            4*c-3             0;
                    0                 0  -nc*s     0               0               c]; %STM matrix

                newSTM(:,:,tt) = blkdiag(STM,STM,STM,STM); %24x24
                cellSTM(tt) = {newSTM(:,:,tt)}; %storing the STM in a cell array

                if(tt == 1)
                    STMkk_1 = newSTM(:,:,tt);
                else
                    STMkk_1 = newSTM(:,:,tt) / newSTM(:,:,tt-1);  % STM(tk,tk_1) = STM(tk, 0) / STM(tk_1, 0), conventional KF
                end

                [Pkbar] = sncompensation(SNC,tt,STMkk_1,Pk_1bar,configuration); %State Noise Compensation

                %% MEASUREMENT UPDATE
                % state-observation relationship & C-matrix
                % (12x18) partial derivatvies wrt the state with the DGPS+ ISLs with the
                % relative positions rho12, rho13, rho23, rho13, rho24 (constant, so no need to evaluate over Xdot).
                Ck = zeros(15,24); %15x24
                rho12 = distvec12(tt,:); %norm of the dist to get the magn.
                rho23 = distvec23(tt,:);
                rho34 = distvec34(tt,:);
                rho41 = distvec41(tt,:);
                rho13 = distvec13(tt,:);
                rho24 = distvec24(tt,:);

                Ck(1:3,1:3)   = eye(3);              %x1,y1,z1, DGPS
                Ck(4:6,7:9)   = eye(3);              %x2,y2,z2, DGPS
                Ck(7:9,13:15) = eye(3);              %x3,y3,z3, DGPS
                Ck(10,1:3)    = (vec12(tt,:)/rho12); %ISLS, per each iteration
                Ck(10,7:9)    = - Ck(10,1:3);        %ISLS, per each iteration
                Ck(11,7:9)    = (vec23(tt,:)/rho23); %ISLS, per each iteration
                Ck(11,13:15)  = -Ck(2,7:9);          %ISLS, per each iteration
                Ck(12,13:15)  = (vec34(tt,:)/rho34); %ISLS, per each iteration
                Ck(12,19:21)  = -Ck(3,13:15);        %ISLS, per each iteration
                Ck(13,1:3)    = -vec41(tt,:)/rho41;  %ISLS, per each iteration
                Ck(13,19:21)  = -Ck(13,1:3);         %ISLS, per each iteration
                Ck(14,1:3)    = vec13(tt,:)/rho13;   %ISLS, per each iteration
                Ck(14,13:15)  = -Ck(14,1:3);         %ISLS, per each iteration
                Ck(15,7:9)    = vec24(tt,:)/rho24;   %ISLS, per each iteration
                Ck(15,19:21)  = - Ck(15,7:9);        %ISLS, per each iteration

                CkDGPS = Ck(1:9,1:24); %DGPS
                CkISLS = Ck(10:15,1:24); %ISLS
                ck = {CkDGPS CkISLS Ck}; %total system

                %measurement-noise covariance matrix
                d = [repelem(sigmaDGPS^2,9) repelem(sigmaISL^2,6)];
                Rk = diag(d); %initialising Rk, 15x15
                RkDGPS = Rk(1:9,1:9); %9x9
                RkISLS = Rk(10:15,10:15); %6x6
                rk = {RkDGPS RkISLS Rk}; %all the cases in a single cell array per each iteration
                Kk = (Pkbar*ck{z}')*((rk{z}+(ck{z}*Pkbar*ck{z}'))^-1); % (18x12) Kalman Gain

                % CKF specific to evaluate the standard deviations
                % a-posteriori covariance matrix
                Pk(:,:,tt)  = (I-(Kk*ck{z}))*Pkbar; %24x24, for the next iteration

                P12(:,:,tt) = Pk(1:6,1:6,tt)     + Pk(7:12,7:12,tt)   - Pk(1:6,7:12,tt)    - Pk(7:12,1:6,tt);             %selecting the components for the sigma^2
                P23(:,:,tt) = Pk(7:12,7:12,tt)   + Pk(13:18,13:18,tt) - Pk(7:12,13:18,tt)  - Pk(13:18,7:12,tt);     %selecting the components for the sigma^2
                P34(:,:,tt) = Pk(13:18,13:18,tt) + Pk(19:24,19:24,tt) - Pk(13:18,19:24,tt) - Pk(19:24,13:18,tt); %selecting the components for the sigma^2
                P41(:,:,tt) = Pk(19:24,19:24,tt) + Pk(1:6,1:6,tt)     - Pk(1:6,19:24,tt)   - Pk(19:24,1:6,tt);         %selecting the components for the sigma^2
                P13(:,:,tt) = Pk(1:6,1:6,tt)     + Pk(13:18,13:18,tt) - Pk(1:6,13:18,tt)   - Pk(13:18,1:6,tt);
                P24(:,:,tt) = Pk(7:12,7:12,tt)   + Pk(19:24,19:24,tt) - Pk(1:6,19:24,tt)   - Pk(19:24,1:6,tt);

                % store standard deviations
                %sig1(:,tt)  = diag(Pk(:,:,tt));  %sigma^2
                sig12(:,tt) = diag(P12(:,:,tt)); %sigma^2
                sig23(:,tt) = diag(P23(:,:,tt)); %sigma^2
                sig34(:,tt) = diag(P34(:,:,tt)); %sigma^2
                sig41(:,tt) = diag(P41(:,:,tt)); %sigma^2
                sig13(:,tt) = diag(P13(:,:,tt)); %sigma^2
                sig24(:,tt) = diag(P24(:,:,tt)); %sigma^2

                RMS12 = 3 * sqrt(sig12(1,:) + sig12(2,:) + sig12(3,:)); %root mean square, 3-standard dev.
                RMS23 = 3 * sqrt(sig23(1,:) + sig23(2,:) + sig23(3,:)); %root mean square, 3-standard dev.
                RMS34 = 3 * sqrt(sig34(1,:) + sig34(2,:) + sig34(3,:)); %root mean square, 3-standard dev.
                RMS41 = 3 * sqrt(sig41(1,:) + sig41(2,:) + sig41(3,:)); %root mean square, 3-standard dev.
                RMS13 = 3 * sqrt(sig13(1,:) + sig13(2,:) + sig13(3,:)); %root mean square, 3-standard dev.
                RMS24 = 3 * sqrt(sig24(1,:) + sig24(2,:) + sig24(3,:)); %root mean square, 3-standard dev.

                RMS{z} = [RMS12' RMS23' RMS34' RMS41' RMS13' RMS24']; %saving DGPS, ISLS , DGPS + ISLS as a cell arrays, each has 1000x3 elements

                %% Update a-priori state deviation and covariance matrix before moving to next obs.
                Pk_1bar = Pk(:,:,tt); %24x24

            end % end of observation loop

            %% Plots
            %% Stand-Dev Check
            if control_plot == 1
                plotcase2(z,tvec,RMS12,RMS23,RMS34,RMS41,RMS13,RMS24,orbit_type,ToF,SNC) %figure 2, figure 3, figure 4
            end
        end %closing z loop
        
        plotXALLconfig2(z,tvec,RMS,orbit_type,ToF,SNC)

        %single plot ISLS vs (DGPS & ISLs)
        if control_plot == 1
            plotcase2DGPSvsDGPSISLS(tvec,RMS,orbit_type,ToF,SNC) %figure 5
        end

    case 3 %4 sat + 1 sat central that analysis the positions.

        n = 30;  % no. of state variables, (6 variables per each satellite, 5 satellites in total),(x,y,z,xdot,ydot,zdot)
        I = eye(n);
        m = length(tvec);     % no. of observations
        %% nominal parameters
        sigmaDGPS = 0.10; %[m]
        sigmaISL = 10^-9; %[m]

        %% least squares solution
        prePk_1bar = diag([1e4,1e4,1e4,1,1,1]);   % a-priori covariance matrix, sigma^2
        Pk_1bar = blkdiag(prePk_1bar,prePk_1bar,prePk_1bar,prePk_1bar,prePk_1bar); %30x30

        %% Integrate Reference Trajectory & state transition matrix
        %STM = zeros(30,30); %initialising the STM matrix
        cellSTM = {};
        % Loop over measurements
        for z = 1:3
            for tt = 1:m

                %Retrieve STMs
                %for i = 1 : length(tvec)
                c = cos(nc*tvec(tt));
                s = sin(nc*tvec(tt));

                STM = [4-3*c             0    0      s/nc           (2/nc)*(1-c)      0;
                    6*(s-nc*tvec(tt)) 1    0 -2*(1-c)/nc   (4*s-3*nc*tvec(tt))/nc  0;
                    0                 0    c       0               0             s/nc;
                    3*nc*s            0    0       c              2*s              0;
                    -6*nc*(1-c)        0    0     -2*s            4*c-3             0;
                    0                 0  -nc*s     0               0               c]; %STM matrix

                newSTM(:,:,tt) = blkdiag(STM,STM,STM,STM,STM); %30x30
                cellSTM(tt) = {newSTM(:,:,tt)}; %storing the STM in a cell array

                if(tt == 1)
                    STMkk_1 = newSTM(:,:,tt);
                else
                    STMkk_1 = newSTM(:,:,tt) / newSTM(:,:,tt-1);  % STM(tk,tk_1) = STM(tk, 0) / STM(tk_1, 0), conventional KF
                end

                [Pkbar] = sncompensation(SNC,tt,STMkk_1,Pk_1bar,configuration); %State Noise Compensation

                %% MEASUREMENT UPDATE
                % state-observation relationship & C-matrix
                % (15x30) partial derivatvies wrt the state with the DGPS + (4x30) ISLs with the
                % relative positions rho12, rho13, rho23, rho13, rho24 (constant, so no need to evaluate over Xdot).
                Ck = zeros(19,30); %19x30
                rho10 = distvec10(tt,:); %norm of the dist to get the magn.
                rho20 = distvec20(tt,:);
                rho30 = distvec30(tt,:);
                rho40 = distvec40(tt,:);

                Ck(1:3,1:3) = eye(3);      %x1,y1,z1, DGPS
                Ck(4:6,7:9) = eye(3);      %x2,y2,z2, DGPS
                Ck(7:9,13:15) = eye(3);    %x3,y3,z3, DGPS
                Ck(10:12,19:21) = eye(3);  %x4,y4,z4, DGPS
                Ck(13:15,25:27) = eye(3);  %x0,y0,z0, DGPS
                Ck(16,1:3)   = vec10(tt,:)/rho10; %ISLS
                Ck(16,25:27) = -Ck(10,1:3);       %ISLS
                Ck(17,7:9)   = vec20(tt,:)/rho20; %ISLS
                Ck(17,25:27) = -Ck(11,7:9);       %ISLS
                Ck(18,13:15) = vec30(tt,:)/rho30; %ISLS
                Ck(18,25:27) = -Ck(12,13:15);     %ISLS
                Ck(19,19:21) = vec40(tt,:)/rho40; %ISLS
                Ck(19,25:27) = -Ck(13,19:21);     %ISLS

                CkDGPS = Ck(1:15,1:30);  %15x30
                CkISLS = Ck(16:19,1:30); %4x30
                ck = {CkDGPS CkISLS Ck};

                %measurement-noise covariance matrix
                d = [repelem(sigmaDGPS^2,15) repelem(sigmaISL^2,4)];
                Rk = diag(d); %initialising Rk, 19x19
                RkDGPS = Rk(1:15,1:15); %15x15
                RkISLS = Rk(16:19,16:19); %4x4
                rk = {RkDGPS RkISLS Rk}; %all the cases in a single cell array per each iteration
                Kk = (Pkbar*ck{z}')*((rk{z}+(ck{z}*Pkbar*ck{z}'))^-1); % (18x12) Kalman Gain

                %CKF specific to evaluate the standard deviations
                %a-posteriori covariance matrix
                Pk(:,:,tt)  = (I-(Kk*ck{z}))*Pkbar; %30x30, for the next iteration

                P10(:,:,tt) = Pk(1:6,1:6,tt)     + Pk(25:30,25:30,tt) - Pk(1:6,25:30,tt)   - Pk(25:30,1:6,tt);   %selecting the components for the sigma^2
                P20(:,:,tt) = Pk(7:12,7:12,tt)   + Pk(25:30,25:30,tt) - Pk(7:12,25:30,tt)  - Pk(25:30,7:12,tt);  %selecting the components for the sigma^2
                P30(:,:,tt) = Pk(13:18,13:18,tt) + Pk(25:30,25:30,tt) - Pk(13:18,25:30,tt) - Pk(25:30,13:18,tt); %selecting the components for the sigma^2
                P40(:,:,tt) = Pk(19:24,19:24,tt) + Pk(25:30,25:30,tt) - Pk(19:24,25:30,tt) - Pk(25:30,19:24,tt); %selecting the components for the sigma^2

                %store standard deviations
                sig10(:,tt) = diag(P10(:,:,tt)); %sigma^2
                sig20(:,tt) = diag(P20(:,:,tt)); %sigma^2
                sig30(:,tt) = diag(P30(:,:,tt)); %sigma^2
                sig40(:,tt) = diag(P40(:,:,tt)); %sigma^2

                RMS10 = 3 * sqrt(sig10(1,:) + sig10(2,:) + sig10(3,:)); %root mean square, 3-standard dev.
                RMS20 = 3 * sqrt(sig20(1,:) + sig20(2,:) + sig20(3,:)); %root mean square, 3-standard dev.
                RMS30 = 3 * sqrt(sig30(1,:) + sig30(2,:) + sig30(3,:)); %root mean square, 3-standard dev.
                RMS40 = 3 * sqrt(sig40(1,:) + sig40(2,:) + sig40(3,:)); %root mean square, 3-standard dev.

                RMS{z}  = [RMS10' RMS20' RMS30' RMS40']; %saving DGPS, ISLS , DGPS + ISLS as a cell array, each has 1000x3 elements

                %% Update a-priori state deviation and covariance matrix before moving to next obs.
                Pk_1bar = Pk(:,:,tt); %24x24

            end % end of observation loop

            %% Plots
            %% Stand-Dev Check

            %single plot DGPS + (DGPS & ISLs)
            if control_plot == 1
                plotcase3(z,tvec,RMS10,RMS20,RMS30,RMS40,orbit_type,ToF,SNC) %figure3, figure 4, figure 5
            end
        end %closing z loop

        plotXALLconfig3(z,tvec,RMS,orbit_type,ToF,SNC)

        if control_plot == 1
            plotcase3DGPSvsDGPSISLS(tvec,RMS,orbit_type,ToF,SNC) %figure 6
        end
end

if control_plot == 1
    if configuration == 1   
        plotconfig1AT(tvec,RMS,orbit_type,ToF,SNC) %subplot DGPS vs ISLS vs (DGPS+ISLs)
    elseif configuration == 2
        plotconfig2AT(tvec,RMS,orbit_type,ToF,SNC) %subplot DGPS vs ISLS vs (DGPS+ISLs)
    elseif configuration == 3
        plotconfig3AT(tvec,RMS,orbit_type,ToF,SNC) %subplot DGPS vs ISLS vs (DGPS+ISLs)
    end

hold off;
end

%% Doppler - case 3 (extra) 
% future plan
%
%%% nominal parameters
%sigmaDoppler = 10; %[m/s]


%% frequency achievied with the data filtering 
% RMSDGPS, RMSISLS, RMSDGPS_ISLS
switch configuration
    case 1 %3 sats
        initialRMS = [RMS{1}(1,:)   RMS{2}(1,:)   RMS{3}(1,:)];    %1x9
        oneorbitRMS  = [RMS{1}(7,:)   RMS{2}(7,:)   RMS{3}(7,:)];    %1x9
        onedayRMS = [RMS{1}(103,:) RMS{2}(103,:) RMS{3}(103,:)];  %1x9
        finalRMS   = [RMS{1}(end,:) RMS{2}(end,:) RMS{3}(end,:) ]; %1x9
        fprintf(1,'parameters   |            DGPS           |             ISLs          |           DGPS + ISls    |\n');  
        fprintf('initial RMS    |   %0.10f %0.10f %0.10f    |   %0.10f %0.10f %0.10f    |     %.10f %.10f %.10f    |\n',initialRMS);
        disp('')
        fprintf('initial RMS    |   %0.10f %0.10f %0.10f    |   %0.10f %0.10f %0.10f    |     %.10f %.10f %.10f    |\n',oneorbitRMS);
        disp('')  
        fprintf('initial RMS    |   %0.10f %0.10f %0.10f    |   %0.10f %0.10f %0.10f    |     %.10f %.10f %.10f    |\n',onedayRMS);
        disp('')
        fprintf('final RMS      |   %0.10f %0.10f %0.10f    |   %0.10f %0.10f %0.10f    |     %.10f %.10f %.10f    |\n',finalRMS);
         
    case 2 %4 sats
        fprintf(1, '\n');
        disp(' CASE 2: 4 sats ')
        initialRMS = [RMS{1}(1,:)   RMS{2}(1,:)   RMS{3}(1,:)];    %1x12
        oneorbitRMS  = [RMS{1}(7,:)   RMS{2}(7,:)   RMS{3}(7,:)];    %1x12
        onedayRMS = [RMS{1}(103,:) RMS{2}(103,:) RMS{3}(103,:)];  %1x12
        finalRMS   = [RMS{1}(end,:) RMS{2}(end,:) RMS{3}(end,:) ]; %1x12
        fprintf(1, '\n'); 
        fprintf(1,'parameters   |              DGPS                                  |                     ISLs                     |              DGPS + ISls            |\n');        
        fprintf('initial RMS    |   %0.10f %0.10f %0.10f %0.10f %0.10f %0.10f        |   %0.10f %0.10f %0.10f %0.10f %0.10f %0.10f  |   %.10f %.10f %.10f %.10f %.10f %f  |\n',initialRMS);
        disp('')
        fprintf('1 orbit RMS    |   %0.10f %0.10f %0.10f %0.10f %0.10f %0.10f        |   %0.10f %0.10f %0.10f %0.10f %0.10f %0.10f  |   %.10f %.10f %.10f %.10f %.10f %f  |\n',oneorbitRMS);
        disp('')
        fprintf('24 hours RMS   |   %0.10f %0.10f %0.10f %0.10f %0.10f %0.10f        |   %0.10f %0.10f %0.10f %0.10f %0.10f %0.10f  |   %.10f %.10f %.10f %.10f %.10f %f  |\n',onedayRMS);
        disp('')
        fprintf('final RMS      |   %0.10f %0.10f %0.10f %0.10f %0.10f %0.10f        |    %0.10f %0.10f %0.10f %0.10f %0.10f %.10f  |   %.10f %.10f %.10f %.10f %.10f %f  |\n',finalRMS);

    case 3 %4 sats + 1 central 
        fprintf(1, '\n');
        disp(' CASE 3: 4 sats + 1 central')
        initialRMS = [RMS{1}(1,:)   RMS{2}(1,:)   RMS{3}(1,:)];    %1x12
        oneorbitRMS  = [RMS{1}(7,:)   RMS{2}(7,:)   RMS{3}(7,:)];    %1x12
        onedayRMS = [RMS{1}(103,:) RMS{2}(103,:) RMS{3}(103,:)];  %1x12
        finalRMS   = [RMS{1}(end,:) RMS{2}(end,:) RMS{3}(end,:) ]; %1x12
        fprintf(1,'parameters   |                   DGPS                   |                  ISLs                    |                 DGPS + ISls              |\n');  
        fprintf('initial RMS    |   %0.10f %0.10f %0.10f %0.10f    |   %0.10f %0.10f %0.10f %0.10f    |   %.10f %.10f %.10f %.10f    |\n',initialRMS);
        disp('')
        fprintf('1 orbit RMS    |   %0.10f %0.10f %0.10f %0.10f    |   %0.10f %0.10f %0.10f %0.10f    |   %.10f %.10f %.10f %.10f    |\n',oneorbitRMS);
        disp('')
        fprintf('24 hours RMS   |   %0.10f %0.10f %0.10f %0.10f    |   %0.10f %0.10f %0.10f %0.10f    |   %.10f %.10f %.10f %.10f    |\n',onedayRMS);
        disp('')
        fprintf('final RMS      |   %0.10f %0.10f %0.10f %0.10f    |   %0.10f %0.10f %0.10f %0.10f    |   %.10f %.10f %.10f %.10f    |\n',finalRMS);
end

% RMS for the 3 cases
initial_lambdaDGPS = 4*min(initialRMS(1:3)); %the best scenario atm with DGPS
initial_lambdaTOT = 4*min(initialRMS(7:9)); %the best scenario atm with DGPS + ISLs
final_lambdaDGPS = 4*min(finalRMS(1:3)); %the best scenario atm with DGPS
final_lambdaTOT = 4*min(real(finalRMS(7:9))); %the best scenario atm with DGPS + ISLs
initial_fDGPS = 3e8/initial_lambdaDGPS; %speed of light/wavelenght
initial_fTOT = 3e8/initial_lambdaTOT; %speed of light/wavelenght
final_fDGPS = 3e8/final_lambdaDGPS; %max frequency obtained as output of the filter
final_fTOT = 3e8/final_lambdaTOT; %max frequency obtained as output of the filter
opt_f = 4e9; %[Hz]
opt_lambda = 3.00e8/opt_f; %with a 4GHz frequency
fprintf(1, '\n');
fprintf('initial lambda DGPS          |      %f m             | final lambda DGPS               |          %f m               |\n',initial_lambdaDGPS,final_lambdaDGPS);     
fprintf('initial frequency DGPS       |      %.3f MHz         | final frequency                 |          %.3f MHz           |\n',initial_fDGPS/1e6,final_fDGPS/1e6);
fprintf('initial lambda DGPS + ISL    |      %f m             | final lambda DGPS + ISL         |          %f m               |\n',initial_lambdaTOT,final_lambdaTOT);
fprintf('initial frequency DGPS + ISL |      %.3f MHz         | final frequency DGPS + ISL      |          %.3f MHz           |\n',initial_fTOT/1e6,final_fTOT/1e6);
fprintf('optimal lambda               |      %f m             | optimal frequency               |          %.3f MHz           |\n',opt_lambda,opt_f/1e6);

% %% Gaussian Distribution
% if control_plot == 1
%     switch configuration
%         case 1
%             p = zeros(1000,length(distmatrix(1,:))); %initialising vector probability density function
%             figure('name','Gaussian Distribution (PDF)','units','normalized','outerposition',[0.2 0.2 0.8 0.8]);
%             for i = 1:length(distmatrix(1,:))
%                 meanmatrix = mean(distmatrix); %mean per columns
%                 p = 1/(sigmaDGPS*sqrt(2*pi))*exp(-(distmatrix(:,i) - meanmatrix(1,i)).^2/(2*sigmaDGPS^2)) ;
%                 hold on;grid on;
%                 plot(distmatrix(:,i),p)
%                 xline(meanmatrix,'-','$\bar{x}$','FontSize',14,'Interpreter','latex');
%                 set(gca,'ytick',[]);
%                 set(gca,'yticklabel',[]);
%                 xlim([180 187]);
%                 xlabel('relative distances (m)');
%                 ylabel('probability ');
%                 title('PDF for configuration 1');
%                 lgd = legend('p12(x)','p13(x)','p23(x)');
%                 lgd.Interpreter = 'tex';
%             end
%             hold off;
%             if save_plot == 1
%                 saveas(ancestor(gcf,'figure'), fullfile('/Users/ludovicogregori/Desktop/SURREY/FINAL PROJECT/matlab codes/orbit&datafilter/figures/distribution/config1/','gaussian1.png'));
%             end
% 
%         case 2
%             p = zeros(1000,length(distmatrix(1,:)));
%             figure('name','Gaussian Distribution (PDF)','units','normalized','outerposition',[0.2 0.2 0.8 0.8]);
%             for i = 1:length(distmatrix(1,:))
%                 meanmatrix = mean(distmatrix); %mean per columns
%                 p = 1/(sigmaDGPS*sqrt(2*pi))*exp(-(distmatrix(:,i) - meanmatrix(1,i)).^2/(2*sigmaDGPS^2)) ;
%                 hold on;grid on;
%                 plot(tvec,p)
%                 xline(meanmatrix,'-','$\bar{x}$','FontSize',14,'Interpreter','latex');
%                 set(gca,'ytick',[]);
%                 set(gca,'yticklabel',[]);
%                 xlabel('relative distances (m)');
%                 ylabel('probability ');
%                 title('PDF for configuration 1');
%                 lgd = legend('p12(x)','p23(x)','p34(x)','p41(x)','p13(x)','p24(x)');
%                 lgd.Interpreter = 'tex';
%             end
%             hold off;
%             if save_plot == 1
%                 saveas(ancestor(gcf,'figure'), fullfile('/Users/ludovicogregori/Desktop/SURREY/FINAL PROJECT/matlab codes/orbit&datafilter/figures/distribution/config2/','gaussian2.png'));
%             end
% 
%         case 3
%             p = zeros(1000,length(distmatrix(1,:)));
%             figure('name','Gaussian Distribution (PDF)','units','normalized','outerposition',[0.2 0.2 0.8 0.8]);
%             for i = 1:length(distmatrix(1,:))
%                 meanmatrix = mean(distmatrix); %mean per columns
%                 p = 1/(sigmaDGPS*sqrt(2*pi))*exp(-(distmatrix(:,i) - meanmatrix(1,i)).^2/(2*sigmaDGPS^2)) ;
%                 hold on;grid on;
%                 plot(tvec,p(i))
%                 xline(meanmatrix,'-','$\bar{x}$','FontSize',14,'Interpreter','latex');
%                 set(gca,'ytick',[]);
%                 set(gca,'yticklabel',[]);
%                 xlabel('relative distances (m)');
%                 ylabel('probability ');
%                 title('PDF for configuration 1');
%                 lgd = legend('p10(x)','p20(x)','p30(x)','p40(x)');
%                 lgd.Interpreter = 'tex';hold off;
%             end
%             hold off;
%             if save_plot == 1
%                 saveas(ancestor(gcf,'figure'), fullfile('/Users/ludovicogregori/Desktop/SURREY/FINAL PROJECT/matlab codes/orbit&datafilter/figures/distribution/config3/','gaussian3.png'));
%             end
%     end
% end

%% Analysis for the EIRP, orbit initial altitude
%hypothesis or given data

control_plot = 1; %parameter to plot or not the figures

Channel_capacity = 100e6; %Mbitps
eirp_sat = 35; %dBWatt
sat_gain = 36.3; %dBi
handset_power = 0.5; %Watt
handset_gain = 5; %dBi

BW_5G_min = 30e9 ; %GHz
BW_5G_max = 300e9; %GHz
BW_5G = [BW_5G_min BW_5G_max];

communication_BW = 100e6; %MHz

%C Band
C_band_min = 4e9; %GHz
C_band_max = 8e9; %GHz
C_band = [C_band_min C_band_max];

c = 2.998e8; %m/s, speed of light
lambda = sort(c./C_band); % wavelenght [m]
relative_position = lambda/4; % initial positions range between two close-satellites in orbit [m]

constellation_altitude = 500e3; %km
Boltz_cost = -228.6 ; % dB W/KelvinHz

atmospheric_loss = 0.49 ; % 0.49 dB, we are overstimating the loss
pointing_loss = 0.55 ; % 0.55 dB we are overstimating the loss

system_temp = 249; %Kelvin, standard value, ADD THE PAPER
system_temp_dB = 10*log10(system_temp);
G_T_rx = handset_gain - system_temp_dB; %dB/K, figure of merite of the receiver, gain over system noise temp.

% FSPL =sort(20*log10(4*pi*constellation_altitude./lambda)); %dB, just deciBel
% G_T_rx = handset_gain - system_temp_dB; %dB/K, figure of merite of the receiver, gain over system noise temp.
% SNR = eirp_sat - (FSPL + atmospheric_loss + pointing_loss) +G_T_rx - Boltz_cost - 10*log10(communication_BW); %dB, ratio of powers, carrier and noise
% Shannon_capacity = communication_BW.*log2(1+SNR);

%% SIMULATION FOR FINAL RESULTS 
h = [];

%simulation 1 
BW = 100e6; %4g bandwidth
fmin = 2.444; %Ghz
fmax = 3.814; %GHz
lambdamin = 0.122; %m for fmin
lambdamax = 0.078; %m for fmax
lambda = [lambdamin lambdamax];

FSPL =sort(20*log10(4*pi*constellation_altitude./lambda)); %dB, just deciBel
SNR = eirp_sat - (FSPL + atmospheric_loss + pointing_loss) +G_T_rx - Boltz_cost - 10*log10(BW); %dB, ratio of powers, carrier and noise
SNR1 = SNR;
capacity = BW.*log2(1+SNR); %bit per sec

s1lambdaminBW = @(EIRP)(lambdamin/(4*pi)).*10.^((-2^(Channel_capacity./BW)+1+EIRP-atmospheric_loss-pointing_loss+G_T_rx-Boltz_cost-10*log10(BW))/20);
s1lambdamaxBW = @(EIRP)(lambdamax/(4*pi)).*10.^((-2^(Channel_capacity./BW)+1+EIRP-atmospheric_loss-pointing_loss+G_T_rx-Boltz_cost-10*log10(BW))/20);

figure('name','Simulation Comparisons EIRPs and Altitudes','units','normalized','outerposition',[0.2 0.2 0.8 0.8])
subplot(1,3,1)
hold on;grid on; grid minor;
fplot(s1lambdaminBW,'.-b'); %blue curve
fplot(s1lambdamaxBW,'.-r'); %red curve
title('EIRP and Altitude: Return Link','FontSize',14)
yline(490e3,'-','initial mission altitude');
yline(512e3,'-','BlueWalker 3 altitude');
yline(398e3,'-','StarLink altitude');
xlabel('EIRP (dBW)');
ylabel('Altitude (m)')
lgd = legend('\lambda_{min}','\lambda_{max}','Location','best');
lgd.Interpreter = 'tex';
axis([0 40 0 750e3]) %limit for the axis: eirp on X, altitude in meter on Y (thats why e3=km)
title('Simulation 1');

disp(' Simulation 1')
fprintf('FSPL     |      %f dB        |\n',FSPL);     
fprintf('SNR      |      %.3f dB      |\n',SNR);
fprintf('Capacity |      %f Mbps          |\n',capacity/1e6);

int_fun3 = @(EIRP)(lambda(1)/(4*pi)).*10.^((-2.^(Channel_capacity./BW)+1+EIRP-atmospheric_loss-pointing_loss+G_T_rx-Boltz_cost-10*log10(BW))/20) - constellation_altitude;
eirp3 = fsolve(int_fun3,70);
int_fun4 = @(EIRP)(lambda(2)/(4*pi)).*10.^((-2.^(Channel_capacity./BW)+1+EIRP-atmospheric_loss-pointing_loss+G_T_rx-Boltz_cost-10*log10(BW))/20) - constellation_altitude;
eirp4 = fsolve(int_fun4,70);
eirp_sol2 = sort([eirp3 eirp4]); %intersection between eirp functions and 700km altitude for our mission.
fprintf('EIRPs |      %f dBW         |\n',eirp_sol2);

int_fun3 = @(EIRP)(lambda(1)/(4*pi)).*10.^((-2.^(Channel_capacity./BW)+1+EIRP-atmospheric_loss-pointing_loss+G_T_rx-Boltz_cost-10*log10(BW))/20) - 400e3;
eirp3 = fsolve(int_fun3,70);
int_fun4 = @(EIRP)(lambda(2)/(4*pi)).*10.^((-2.^(Channel_capacity./BW)+1+EIRP-atmospheric_loss-pointing_loss+G_T_rx-Boltz_cost-10*log10(BW))/20) - 400e3;
eirp4 = fsolve(int_fun4,70);
eirp_sol2 = sort([eirp3 eirp4]); %intersection between eirp functions and 700km altitude for our mission.
fprintf('EIRPs at 400km (simul. 1)|      %f dBW         |\n',eirp_sol2);

nt_fun3 = @(EIRP)(lambda(1)/(4*pi)).*10.^((-2.^(Channel_capacity./BW)+1+EIRP-atmospheric_loss-pointing_loss+G_T_rx-Boltz_cost-10*log10(BW))/20) - 700e3;
eirp3 = fsolve(int_fun3,70);
int_fun4 = @(EIRP)(lambda(2)/(4*pi)).*10.^((-2.^(Channel_capacity./BW)+1+EIRP-atmospheric_loss-pointing_loss+G_T_rx-Boltz_cost-10*log10(BW))/20) - 700e3;
eirp4 = fsolve(int_fun4,70);
eirp_sol2 = sort([eirp3 eirp4]); %intersection between eirp functions and 700km altitude for our mission.
fprintf('EIRPs at 700km (simul. 1)|      %f dBW         |\n',eirp_sol2);

%simulation 2 
BW = 100e6; %MHz
fmin = 4.402; %Ghz
fmax = 10.135; %GHz
lambdamin = 0.068; %m for fmin
lambdamax = 0.029; %m for fmax
lambda = [lambdamin lambdamax];

FSPL =sort(20*log10(4*pi*constellation_altitude./lambda)); %dB, just deciBel
SNR = eirp_sat - (FSPL + atmospheric_loss + pointing_loss) + G_T_rx - Boltz_cost - 10*log10(BW); %dB, ratio of powers, carrier and noise
SNR23 = SNR;
capacity = BW.*log2(1+SNR);

s1lambdaminBW = @(EIRP)(lambdamin/(4*pi)).*10.^((-2^(Channel_capacity./BW)+1+EIRP-atmospheric_loss-pointing_loss+G_T_rx-Boltz_cost-10*log10(BW))/20);
s1lambdamaxBW = @(EIRP)(lambdamax/(4*pi)).*10.^((-2^(Channel_capacity./BW)+1+EIRP-atmospheric_loss-pointing_loss+G_T_rx-Boltz_cost-10*log10(BW))/20);

subplot(1,3,2)
hold on;grid on; grid minor;
fplot(s1lambdaminBW,'.-b'); %blue curve
fplot(s1lambdamaxBW,'.-r'); %red curve
title('EIRP and Altitude: Return Link','FontSize',14)
yline(490e3,'-','initial mission altitude');
yline(512e3,'-','BlueWalker 3 altitude');
yline(398e3,'-','StarLink altitude');
xlabel('EIRP (dBW)');
ylabel('Altitude (m)')
lgd = legend('\lambda_{min}','\lambda_{max}','Location','best');
lgd.Interpreter = 'tex';
axis([0 40 0 750e3]) %limit for the axis: eirp on X, altitude in meter on Y (thats why e3=km)
title('Simulation 2')

disp(' Simulation 2')
fprintf('FSPL     |      %f dB        |\n',FSPL);     
fprintf('SNR      |      %.3f dB      |\n',SNR);
fprintf('Capacity |      %f Mbps          |\n',capacity/1e6);

int_fun3 = @(EIRP)(lambda(1)/(4*pi)).*10.^((-2.^(Channel_capacity./BW)+1+EIRP-atmospheric_loss-pointing_loss+G_T_rx-Boltz_cost-10*log10(BW))/20) - constellation_altitude;
eirp3 = fsolve(int_fun3,70);
int_fun4 = @(EIRP)(lambda(2)/(4*pi)).*10.^((-2.^(Channel_capacity./BW)+1+EIRP-atmospheric_loss-pointing_loss+G_T_rx-Boltz_cost-10*log10(BW))/20) - constellation_altitude;
eirp4 = fsolve(int_fun4,70);
eirp_sol2 = sort([eirp3 eirp4]); %intersection between eirp functions and 700km altitude for our mission.
fprintf('EIRPs |      %f dBW         |\n',eirp_sol2);

int_fun3 = @(EIRP)(lambda(1)/(4*pi)).*10.^((-2.^(Channel_capacity./BW)+1+EIRP-atmospheric_loss-pointing_loss+G_T_rx-Boltz_cost-10*log10(BW))/20) - 400e3;
eirp3 = fsolve(int_fun3,70);
int_fun4 = @(EIRP)(lambda(2)/(4*pi)).*10.^((-2.^(Channel_capacity./BW)+1+EIRP-atmospheric_loss-pointing_loss+G_T_rx-Boltz_cost-10*log10(BW))/20) - 400e3;
eirp4 = fsolve(int_fun4,70);
eirp_sol2 = sort([eirp3 eirp4]); %intersection between eirp functions and 700km altitude for our mission.
fprintf('EIRPs at 400km (simul. 2)|      %f dBW         |\n',eirp_sol2);

nt_fun3 = @(EIRP)(lambda(1)/(4*pi)).*10.^((-2.^(Channel_capacity./BW)+1+EIRP-atmospheric_loss-pointing_loss+G_T_rx-Boltz_cost-10*log10(BW))/20) - 700e3;
eirp3 = fsolve(int_fun3,70);
int_fun4 = @(EIRP)(lambda(2)/(4*pi)).*10.^((-2.^(Channel_capacity./BW)+1+EIRP-atmospheric_loss-pointing_loss+G_T_rx-Boltz_cost-10*log10(BW))/20) - 700e3;
eirp4 = fsolve(int_fun4,70);
eirp_sol2 = sort([eirp3 eirp4]); %intersection between eirp functions and 700km altitude for our mission.
fprintf('EIRPs at 700km (simul. 2)|      %f dBW         |\n',eirp_sol2);

%simulation 3
BW = 100e6; %MHz
fmin = 4.356; %Ghz
fmax = 10.206; %GHz
lambdamin = 0.068; %m for fmin
lambdamax = 0.029; %m for fmax
lambda = [lambdamin lambdamax];

FSPL =sort(20*log10(4*pi*constellation_altitude./lambda)); %dB, just deciBel
SNR = eirp_sat - (FSPL + atmospheric_loss + pointing_loss) +G_T_rx - Boltz_cost - 10*log10(BW); %dB, ratio of powers, carrier and noise
capacity = BW.*log2(1+SNR);

s1lambdaminBW = @(EIRP)(lambdamin/(4*pi)).*10.^((-2^(Channel_capacity./BW)+1+EIRP-atmospheric_loss-pointing_loss+G_T_rx-Boltz_cost-10*log10(BW))/20);
s1lambdamaxBW = @(EIRP)(lambdamax/(4*pi)).*10.^((-2^(Channel_capacity./BW)+1+EIRP-atmospheric_loss-pointing_loss+G_T_rx-Boltz_cost-10*log10(BW))/20);

subplot(1,3,3)
hold on;grid on; grid minor;
fplot(s1lambdaminBW,'.-b'); %blue curve
fplot(s1lambdamaxBW,'.-r'); %red curve
title('EIRP and Altitude: Return Link','FontSize',14)
yline(490e3,'-','initial mission altitude');
yline(512e3,'-','BlueWalker 3 altitude');
yline(398e3,'-','StarLink altitude');
xlabel('EIRP (dBW)');
ylabel('Altitude (m)')
lgd = legend('\lambda_{min}','\lambda_{max}','Location','best');
lgd.Interpreter = 'tex';
axis([0 40 0 750e3]) %limit for the axis: eirp on X, altitude in meter on Y (thats why e3=km)
title('Simulation 3')

disp(' Simulation 3')
fprintf('FSPL     |      %f dB        |\n',FSPL);     
fprintf('SNR      |      %.3f dB      |\n',SNR);
fprintf('Capacity |      %f Mbps          |\n',capacity/1e6);

int_fun3 = @(EIRP)(lambda(1)/(4*pi)).*10.^((-2.^(Channel_capacity./BW)+1+EIRP-atmospheric_loss-pointing_loss+G_T_rx-Boltz_cost-10*log10(BW))/20) - constellation_altitude;
eirp3 = fsolve(int_fun3,70);
int_fun4 = @(EIRP)(lambda(2)/(4*pi)).*10.^((-2.^(Channel_capacity./BW)+1+EIRP-atmospheric_loss-pointing_loss+G_T_rx-Boltz_cost-10*log10(BW))/20) - constellation_altitude;
eirp4 = fsolve(int_fun4,70);
eirp_sol2 = sort([eirp3 eirp4]); %intersection between eirp functions and 700km altitude for our mission.
fprintf('EIRPs |      %f dBW         |\n',eirp_sol2);

%% outrage probability P
noservice = 180; %seconds of no service per year
year2sec = 365*24*3600; %number of seconds in 1 year

% slow-fading channel
r = 100 * 10^6; % 100Mbitps, required threshold information rate.
BW1 = [2.444, 3.814] * 10^9;
BW23 = [4.356, 10.206] * 10^9;
h1  = sqrt((2.^(r./BW)-1)./(SNR1));  %fading coefficient
h23 = real(sqrt((2.^(r./BW)-1)./(SNR23))); %fading coefficient
% outrage probability via guassian distribution

% h = linspace(0,1000,10000);
% Cfading = communication_BW.*log2(1+h.^2*max(test_SNR)); %function handle with h varying.
Cfadingsim1min  = @(h)(BW1(1).*log2(1+h.^2*max(SNR1))); %function handle with h varying.
Cfadingsim1max  = @(h)(BW1(2).*log2(1+h.^2*max(SNR1))); %function handle with h varying.
Cfadingsim23min = @(h)(BW23(1).*log2(1+h.^2*max(SNR23))); %function handle with h varying.
Cfadingsim23max = @(h)(BW23(2).*log2(1+h.^2*max(SNR23))); %function handle with h varying.

figure('name','Fading Channel with different fading coefficients','units','normalized','outerposition',[0.2 0.2 0.8 0.8]);
subplot(1,2,1)
hold on;grid on; grid minor;
fplot(Cfadingsim1min,'LineStyle','--','LineWidth',2)
fplot(Cfadingsim1max,'LineStyle','--','LineWidth',2)
yline(100e6,'-','minimum capacity','FontSize',20);
xline(h1(1),'-','h at 2.44 GHz','FontSize',20);
xline(h1(2),'-','h at 3.814 GHz','FontSize',20);
% ylim([0 150e6]);
% xlim([-0.05 +0.05]);
% yticks(0:20e6:150e6);
% xticks(-0.05:0.01:+0.05);
xlabel('fading coefficient','FontSize',20);
ylabel('fading channel capacity','FontSize',20)
lgd = legend('Cfading at 2.44 GHz','Cfading at 3.814 GHz','Location','northwest','FontSize',20);
lgd.Interpreter = 'tex';

subplot(1,2,2)
hold on;grid on; grid minor;
fplot(Cfadingsim23min,'LineStyle',':','LineWidth',2)
fplot(Cfadingsim23max,'LineStyle',':','LineWidth',2)
%ylim([0 150e6]);
%xlim([-0.05 +0.05]);
%yticks(0:20e6:150e6);
%xticks(-0.05:0.01:+0.05);
yline(100e6,'-','minimum capacity','FontSize',20);
xline(h23(1),'-','h at 4.356 GHz','FontSize',20);
xline(h23(2),'-','h at 10.206 GHz','FontSize',20);
xlabel('fading coefficient','FontSize',20);
ylabel('fading channel capacity','FontSize',20)
lgd = legend('Cfading at 4.356 GHz','Cfading at 10.206 GHz','Location','northwest','FontSize',20);
lgd.Interpreter = 'tex';
hold off;

sigma = 100e6; %100Mbps
h8GHz = linspace(-0.0500,0.0500,1000);
Cfading = zeros(4,1000);
Cfading(1,:) = BW1(1).*log2(1+h8GHz.^2*max(SNR1)); %simulation C for gaussian distribution
Cfading(2,:) = BW1(2).*log2(1+h8GHz.^2*max(SNR1)); %simulation C for gaussian distribution
Cfading(3,:) = BW23(1).*log2(1+h8GHz.^2*max(SNR23)); %simulation C for gaussian distribution
Cfading(4,:) = BW23(2).*log2(1+h8GHz.^2*max(SNR23)); %simulation C for gaussian distribution

meanCfad = mean(Cfading,2); %mean

p = zeros(4,1000);
p(1,:) = (1/(sigma*sqrt(2*pi))) * exp(-(Cfading(1,:) - meanCfad(1)).^2/(2*sigma^2));
p(2,:) = (1/(sigma*sqrt(2*pi))) * exp(-(Cfading(2,:) - meanCfad(2)).^2/(2*sigma^2));
p(3,:) = (1/(sigma*sqrt(2*pi))) * exp(-(Cfading(3,:) - meanCfad(3)).^2/(2*sigma^2));
p(4,:) = (1/(sigma*sqrt(2*pi))) * exp(-(Cfading(4,:) - meanCfad(4)).^2/(2*sigma^2));

figure('name','Gaussian Distribution for the fading channel','units','normalized','outerposition',[0.2 0.2 0.8 0.8]);

subplot(1,2,1)
hold on;grid on;
plot(Cfading(1,:),p(1,:),'Color','r','LineStyle','-')
plot(Cfading(2,:),p(2,:),'Color','r','LineStyle','-')
xline([meanCfad(1) meanCfad(2)],'-',{'mean sim1 min','mean sim1 max'},'FontSize',20);
set(gca,'ytick',[])
set(gca,'yticklabel',[])
xlabel('fading channel capacity','FontSize',20);
ylabel('probability p(C)','FontSize',20);
lgd = legend('p(x) at 2.44 GHz','p(x) at 3.814 GHz','Location','northwest','FontSize',20);
lgd.Interpreter = 'tex';

subplot(1,2,2)
hold on;grid on;
plot(Cfading(3,:),p(3,:),'Color','b','LineStyle','--')
plot(Cfading(4,:),p(4,:),'Color','b','LineStyle','--')
% ylim([0 150e6]);
% xlim([-0.05 +0.05]);
% yticks(0:20e6:150e6);
% xticks(-0.05:0.01:+0.05);
% yline(100e6,'-','minimum capacity');
xline([meanCfad(3) meanCfad(4)],'-',{'mean sim2 min','mean sim2 max'},'FontSize',20);
% xline(h(2),'-','h at 8GHz');
set(gca,'ytick',[])
set(gca,'yticklabel',[])
xlabel('fading channel capacity','FontSize',20);
ylabel('probability p(C)','FontSize',20);
lgd = legend('p(x) at 4.356 GHz','p(x) at 10.206 GHz','Location','northwest','FontSize',20);
lgd.Interpreter = 'tex';
hold off;
