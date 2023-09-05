function [Centres3Dtot_reshaped] = arraydesign(positions,tvec_design)

%% PATCHED ARRAY DESIGN
%clear all; 
%close all;
%clc; 
format lonG; %avoiding old data, variables, figures, ...

%Set Defaults
set(0,'DefaultTextInterpreter','latex');
set(0,'DefaultLegendInterpreter','latex');
set(0,'DefaultAxesFontSize',20);
set(0,'DefaultLineLineWidth',2);

% Conversion factors
d2r = pi/180;   % from degrees to radians (deg^-1);
r2d = 180/pi;   % from radians to degrees (deg);

%% Elements & System
%antpar stands for antennas' parameter 
antpar.C_band_min = 4e9; %GHz
C_band_min = antpar.C_band_min;
c = 2.998e8; %m/s, speed of light
antpar.lambda = c/C_band_min; % wavelenght [m]
lambda = antpar.lambda;
antpar.relative_position = lambda/4; % initial positions range between two close-satellites in orbit [m]
relative_position = antpar.relative_position;
antpar.element_radius = lambda(1)/5; %metres
element_radius = antpar.element_radius;
antpar.elements_Callum = 553; %number of the single patched antenna
elements_Callum = antpar.elements_Callum;
antpar.single_system_radius = 0.5; %metres
single_system_radius = antpar.single_system_radius; 

%controls
global control_plot gif_plot

%circular patched antennas
extR = 0.5; %exernal radius
%theta = (0:pi/50:2*pi); %2pi
theta = linspace(0,2*pi,100); %100 elements
% xinit = 5500;
% yinit = xinit;
% zinit = xinit;

xinit = positions(1,2); %first sat, RST, S components along y axis.
yinit = positions(1,3); % T components, along z axis.  
zinit = positions(1,1); % R components, along x axis.

% for plotting a single array
xc = xinit + extR*cos(theta); %x components
yc = yinit + extR*sin(theta); %y components
zc = zinit ; %repelem(zinit,length(xc));

%for deriving all the coordinates of the centres along the orbit
%SAT 1
XC = positions(:,2) + extR * cos(theta); %1000x100
YC = positions(:,3) + extR * sin(theta); %1000x100
ZC = [];
for i = 1 : length(positions(:,1))
ZC = [ZC; repelem(positions(i,1),100)];  %1000x100
end
SAT1 = [XC YC ZC];

%SAT2 
XC2 = positions(:,2) + extR * cos(theta); %1000x100
YC2 = positions(:,3) + extR * sin(theta); %1000x100
ZC2 = [];
for i = 1 : length(positions(:,1))
ZC2 = [ZC2; repelem(positions(i,1),100)];  %1000x100
end
SAT2 = [XC2 YC2 ZC2];

%SAT3
XC3 = positions(:,2) + extR * cos(theta); %1000x100
YC3 = positions(:,3) + extR * sin(theta); %1000x100
ZC3 = [];
for i = 1 : length(positions(:,1))
ZC3 = [ZC3; repelem(positions(i,1),100)];  %1000x100
end
SAT3 = [XC3 YC3 ZC3];

POSSAT = {SAT1 SAT2 SAT3};

antpar.spacer = 2.3; %inter element distance parameter
spacer = antpar.spacer;

%cubesat initial dimension
cubesat_area =  spacer * extR; %m^2 
cubesat_side = sqrt(cubesat_area);

x1 = xinit - cubesat_side/2; %metres
x2 = xinit + cubesat_side/2; 
y1 = yinit - cubesat_side/2; %metres
y2 = yinit + cubesat_side/2; 
z1 = zinit - cubesat_side/2; %metres
z2 = zinit + cubesat_side/2;

% x1 = xinit - cubesat_side/2; %metres
% x2 = xinit + cubesat_side/2; 
% y1 = yinit - cubesat_side/2; %metres
% y2 = yinit + cubesat_side/2; 
% z1 = zinit - cubesat_side/2; %metres
% z2 = zinit + cubesat_side/2;

x = [x1, x2, x2, x1, x1, x1, x2, x2, x1, x1, x1, x1, x2, x2, x2, x2, x1, x1]; %x coordinates
y = [y1, y1, y2, y2, y1, y1, y1, y2, y2, y1, y1, y1, y1, y1, y2, y2, y2, y2]; %y coordinates
z = [z1, z1, z1, z1, z1, z2, z2, z2, z2, z2, z1, z2, z1, z2, z1, z2, z1, z2]; %z coordinates 
% figure('name','S/C top face sides 3D','units','normalized','outerposition',[0.2 0.2 0.8 0.8]); 
% cubesat = plot3(x, y, z, 'k:','LineWidth',3); %cubesat external sides
% xlabel('x [m]','FontSize',16);
% ylabel('y [m]','FontSize',16);
% zlabel('z [m]','FontSize',16);
%axis square;

%Define rows and columns based on min,max xc yc coordinates.
safety_dist = 0.05;
rows = (min(xc)+safety_dist : spacer*element_radius : max(xc)-safety_dist); %vector of the row for the x-components

%ROWS matrix is generating the x components of the antennas' positions.
ROWS = [];
for i = 1 : length(positions(:,1))
ROWS = [ROWS; (min(XC(i,:))+safety_dist : spacer*element_radius : max(XC(i,:))-safety_dist)'];  %27000x1
end

matrows = []; %zeros(3000,27);
matrowsascolumns = []; %zeros(27000,1);
for j = 1:3
    for i = 1:1000
        matrows(i,:) = (min(POSSAT{j}(i,1:100))+safety_dist : spacer*element_radius : max(POSSAT{j}(i,1:100))-safety_dist);
        matrows = [matrows; (matrows(i,:))]; 
    end
end
matrows(end,:) = []; %last vector has the same components of the first one, like in a circle 0=2pi
matrows_reshaped = reshape(matrows,81000,1);% all the x components of the 3 systems in orbit. 

antpar.maxyc = yinit + extR; %it will become max(yc)
maxyc = antpar.maxyc;
antpar.minyc = yinit - extR; %it will become min(yc)
minyc = antpar.minyc;

%minmaxYC matrix is generating the y components of the antennas' positions.
minmaxYC = [];
for i = 1 : length(YC(:,1))
minmaxYC = [minmaxYC; [min(YC(i,:)) max(YC(i,:))]];
end

minmaxYCtot = []; %3000x2
for j = 1:3
    for i = 1:1000
        minmaxYCtot(i,:) = [min(POSSAT{j}(i,101:200)) max(POSSAT{j}(i,101:200))];
        minmaxYCtot = [minmaxYCtot; (minmaxYCtot(i,:))]; 
    end
end
minmaxYCtot(end,:) = []; %last vector has the same components of the first one, like in a circle 0=2pi

[CC] = arraycolumns(antpar); %y components per column 2D
[CC3D] = arraycolumns3D(antpar,minmaxYC); %y components per column 3D
[CC3Dtot] = arraycolumns3Dtot(antpar,minmaxYCtot); %y components per column per 3 systems.

%number of elements per columns
sizeCC = []; 
for j = 1 : length(CC)
    sizeCC(j) = [length(CC{j})];
end
elements_number_check = sum(sizeCC);

%total number of elements per columns per orbit
sizeCC3D = []; %initialising the elements per column vector 
ycentres = []; %initialising the coordinates per element per column vector
for j = 1 : length(CC3D) %j from 1 to 1000 
        for k = 1 : length(CC3D{j}) % k from 1 to each single cell length
            sizeCC3D = [sizeCC3D; length(CC3D{j}{k})'];
            new_ite_y = [CC3D{j}{k}];
            ycentres = [ycentres; new_ite_y'];
        end
end
elements_number3D_check = sum(sizeCC3D);

%total number of elements per columns for the entire system of 3 sats 
sizeCC3Dtot = []; %initialising the elements per column vector 
ycentrestot = []; %initialising the coordinates per element per column vector

for h = 1:3 %h from 1 to 3 to check all the CC3Dtot elements. 
     for j = 1 : length(CC3Dtot{h}) %j from 1 to 1000 
            for k = 1 : length(CC3Dtot{h}{j}) % k from 1 to each single cell length
                sizeCC3Dtot = [sizeCC3Dtot; length(CC3Dtot{h}{j}{k})'];
                new_ite_y_tot = [CC3Dtot{h}{j}{k}];
                ycentrestot = [ycentrestot; new_ite_y_tot'];  %81000x1 cells, each cell has the y components
            end
     end
end     

elements_number3Dtot_check = sum(sizeCC3Dtot);

%patched array centres
Centres = []; %zeros(tot_elements,2), initial position 2D
Centres3D = []; %zeros(elements_number3D_check,3), all positions 3D for 1 satellite along the orbit.
Centres3Dtot = []; %zeros(elements_number3Dtot,3), all positions 3D for 3 satellites along their orbits.
%adding the third coordinate as the radial component. 
%it is changing every "elements number" iterations, cus the array is in a
%single plane.

%% 2D
for i = 1 : length(rows)
    new_ite = [repelem(rows(i) , sizeCC(i))' CC{i}'];
    Centres = [Centres; new_ite]; %all the coordinates of the centres.
end
elements_number = length(Centres);

%% 3D

%x coordinates of the centres along the orbit
for i = 1 : length(ROWS(:,1)) %i from 1 to 27000
            new_ite3D = [repelem(ROWS(i,1) , sizeCC3D(i,1)),]; %CC3D{i}']%repelem(positions(i,1)() , elements_number_check)'];
            Centres3D = [Centres3D; new_ite3D']; %all the coordinates of the centres.
end

%x coordinates of the centres for the entire system of 3 phased array
%systems
for i = 1 : length(matrows_reshaped(:,1)) %i from 1 to 81000, 27000 x each sat 
            new_ite3Dtot = [repelem(matrows_reshaped(i,1) , sizeCC3Dtot(i,1)),];
            Centres3Dtot = [Centres3Dtot; new_ite3Dtot']; %all the coordinates of the centres.
end

%z coordinates of the centres along the orbit
zcentres = [];
for i = 1:length(positions(:,1)) %i from 1 to 1000%
    zcentres = [zcentres; repelem(positions(i,1), elements_number)']; %generating the same radial component per each antennas in a system of 577 elements.
end
Centres3D = [Centres3D ycentres zcentres]; %adding y coordinates and z coordinates

%z coordinates of the centres for the entire system of 3 phased array
zcentrestot = [];
for h = 1:length(POSSAT) %h from 1 to 3, number of elements in POSSAT
     for j = 1:length(POSSAT{h}(:,282)) %i from 1 to 1000%
        zcentrestot = [zcentrestot; repelem(POSSAT{h}(j,282), elements_number)']; %generating the same radial component per each antennas in a system of 577 elements for 1000 times for 3 systems.
     end
end
Centres3Dtot = [Centres3Dtot ycentrestot zcentrestot]; %adding y coordinates and z coordinates
Centres3Dtot_reshaped = [(Centres3Dtot(1:elements_number3D_check,1:3)) (Centres3Dtot(elements_number3D_check+1:2*elements_number3D_check,1:3)) (Centres3Dtot(2*elements_number3D_check+1:3*elements_number3D_check,1:3))];

%external boundaries of the patched elements
%% 2D, dentro i plot per risparmiare ponteza computazionale
coorX = Centres(:,1) + element_radius*cos(theta);
coorY = Centres(:,2) + element_radius*sin(theta);
%coorZ = positions(:,1);

%% 3D, dentro i plot per risparmiare potenza computazionale
coorX3D = Centres3D(:,1) + element_radius*cos(theta);
coorY3D = Centres3D(:,2) + element_radius*sin(theta);
coorZ3D = Centres3D(:,3); 
coor3D = [coorX3D coorY3D coorZ3D]; 

%% 3D full constellations, placing the array in each satellite along the orbit. 
coorX3Dtot = [Centres3Dtot(:,1)] + element_radius*cos(theta);
coorY3Dtot = [Centres3Dtot(:,2)] + element_radius*sin(theta);
coorZ3Dtot = [Centres3Dtot(:,3)];
coor3Dtot  = [coorX3Dtot coorY3Dtot coorZ3Dtot]; %1731000x201, 100 (columns) coordinates to describe the circle for x and y, 1 for the z that is fixed per each time.

if control_plot == 1
    figure('name','Design of a single patched antenna in 2D & 3D','units','normalized','Position', [0.2 0.2 0.6 0.6])
    hold on; grid on;
    xlimplot = [0 xinit+0.5; xinit-0.5 xinit+0.5];
    ylimplot = [0 yinit+0.5; yinit-0.5 yinit+0.5];

    % xlimplot = [yinit/2 1.5*yinit; yinit-1.5*extR yinit+1.5*extR];
    % ylimplot = [zinit/2 1.5*zinit; zinit-1.5*extR zinit+2.0*extR];

    % for i = 1:2
    %   subplot(1,2,i);
        subplot(1,2,1); %2D plotting 
        hold on;grid on;
        %cubesat = plot(x, y, 'k:','LineWidth',3); %cubesat external sides
        patched = plot(xc,yc,'r'); %external boundary of the circular patched antennas
        centres = plot(Centres(:,1),Centres(:,2),'k.','MarkerSize',5); %plotting the centre of each patch antenna
        single_patch = plot(coorX(1:1:end),coorY(1:1:end),'b.','MarkerSize',0.5); %plot function.
        circle = plot(0,'bo','LineWidth',0.05);
        xlabel('Cross-Track [m]','FontSize',16);
        ylabel('Along-Track [m]','FontSize',16);
        xlim(xlimplot(2,:));
        ylim(ylimplot(2,:));
        % xlim(xlimplot(i,:));
        % ylim(ylimplot(i,:));
        axis square;
        legend('external boundary of the antennas','centres of each pacthed antenna','','patched antenna elements','Location','bestoutside','fontsize',15);
        hold off;

        subplot(1,2,2) %3D plotting 
        hold on; grid on
        patched = plot3(xc,yc,repelem(zinit,length(xc)),'r'); %external boundary of the circular patched antennas
        centres = plot3(Centres(:,1),Centres(:,2),repelem(zinit,length(Centres(:,1))),'k.','MarkerSize',5); %plotting the centre of each patch antenna
        single_patch = plot3(coorX(1:1:end),coorY(1:1:end),repelem(zinit,length(coorX(1:1:end))),'b.','MarkerSize',0.5); %plot function.
        xlabel('Cross-Track [m]','FontSize',16);
        ylabel('Along-Track [m]','FontSize',16);
        zlabel('Radial [m]','FontSize',16);
        axis square;
        view(-45,+40);
        hold off;

        saveas(ancestor(gcf,'figure'), fullfile('/Users/ludovicogregori/Desktop/SURREY/FINAL PROJECT/matlab codes/orbit&datafilter/figures/phasedarray', 'phasedarray.png'));
end
 fprintf(1, '\n');
 fprintf('Element radius is equal to %.5f cm.\n',element_radius);
 fprintf('Number of elements is equal to %d.\n',elements_number);
 
%distance between each elements in the same column 
%2 closest elements in the column 1 
fprintf(1, '\n');
d_same_column = sqrt((Centres(2,1)-Centres(1,1))^2+(Centres(2,2)-Centres(1,2))^2) - 2*element_radius;
fprintf('distance between elements in the same column is equal to %.5f cm\n',d_same_column);

%2 closest elements in the same row
fprintf(1, '\n');
d_same_row = sqrt((Centres(10,1)-Centres(1,1))^2+(Centres(12,2)-Centres(1,2))^2) - 2*element_radius;
fprintf('distance between elements in the same row is equal to %.5f cm\n',d_same_row);




