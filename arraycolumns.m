function [CC] = arraycolumns(antpar)

%antpar stands for antennas' parameter 
spacer = antpar.spacer;
element_radius = antpar.element_radius;
maxyc = antpar.maxyc;
minyc = antpar.minyc;

% y1  = (0.40 :spacer*element_radius: 0.70); %first column
% y2  = (0.35 :spacer*element_radius: 0.77);
% y3  = (0.30 :spacer*element_radius: 0.83);
% y4  = (0.25 :spacer*element_radius: 0.85);
% y5  = (0.20 :spacer*element_radius: 0.90);
% y6  = (0.18 :spacer*element_radius: 0.93);
% y7  = (0.15 :spacer*element_radius: 0.95);
% y8  = (0.15 :spacer*element_radius: 0.98);
% y9  = (0.15 :spacer*element_radius: 0.98);
% y10 = (0.12 :spacer*element_radius: 1.00);
% y11 = (0.10 :spacer*element_radius: 1.00);
% y12 = (0.10 :spacer*element_radius: 1.00);
% y13 = (0.10 :spacer*element_radius: 1.00);
% y14 = (0.10 :spacer*element_radius: 1.00); %central column

y1  = (minyc + 0.35 :spacer*element_radius: maxyc - 0.35); %first column
y2  = (minyc + 0.30 :spacer*element_radius: maxyc - 0.30);
y3  = (minyc + 0.25 :spacer*element_radius: maxyc - 0.25);
y4  = (minyc + 0.20 :spacer*element_radius: maxyc - 0.20);
y5  = (minyc + 0.15 :spacer*element_radius: maxyc - 0.15);
y6  = (minyc + 0.13 :spacer*element_radius: maxyc - 0.13);
y7  = (minyc + 0.10 :spacer*element_radius: maxyc - 0.10);
y8  = (minyc + 0.10 :spacer*element_radius: maxyc - 0.10);
y9  = (minyc + 0.10 :spacer*element_radius: maxyc - 0.10);
y10 = (minyc + 0.08 :spacer*element_radius: maxyc - 0.08);
y11 = (minyc + 0.05 :spacer*element_radius: maxyc - 0.05);
y12 = y11; %(0.10 :spacer*element_radius: 1.00);
y13 = y12; %(0.10 :spacer*element_radius: 1.00);
y14 = y13; %(-max(yc) + 0.10 :spacer*element_radius: +max(yc) - 1.00); %central column
y15 = y13;
y16 = y12;
y17 = y11;
y18 = y10;
y19 = y9;
y20 = y8;
y21 = y7;
y22 = y6;
y23 = y5;
y24 = y4;
y25 = y3;
y26 = y2;
y27 = y1; %last column

CC = {y1 y2 y3 y4 y5 y6 y7 y8 y9 y10 y11 y12 y13 y14 y15 y16 y17 y18 y19 y20 y21 y22 y23 y24 y25 y26 y27};