function [CC3D] = arraycolumns3D(antpar,minmaxYC)

%antpar stands for antennas' parameter 
spacer = antpar.spacer;
element_radius = antpar.element_radius;

for i = 1 : length(minmaxYC(:,1))
y1  = (minmaxYC(i,1) + 0.35 :spacer*element_radius: minmaxYC(i,2) - 0.35); %first column
y2  = (minmaxYC(i,1) + 0.30 :spacer*element_radius: minmaxYC(i,2) - 0.30);
y3  = (minmaxYC(i,1) + 0.25 :spacer*element_radius: minmaxYC(i,2) - 0.25);
y4  = (minmaxYC(i,1) + 0.20 :spacer*element_radius: minmaxYC(i,2) - 0.20);
y5  = (minmaxYC(i,1) + 0.15 :spacer*element_radius: minmaxYC(i,2) - 0.15);
y6  = (minmaxYC(i,1) + 0.13 :spacer*element_radius: minmaxYC(i,2) - 0.13);
y7  = (minmaxYC(i,1) + 0.10 :spacer*element_radius: minmaxYC(i,2) - 0.10);
y8  = (minmaxYC(i,1) + 0.10 :spacer*element_radius: minmaxYC(i,2) - 0.10);
y9  = (minmaxYC(i,1) + 0.10 :spacer*element_radius: minmaxYC(i,2) - 0.10);
y10 = (minmaxYC(i,1) + 0.08 :spacer*element_radius: minmaxYC(i,2) - 0.08);
y11 = (minmaxYC(i,1) + 0.05 :spacer*element_radius: minmaxYC(i,2) - 0.05);
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

CC3D{i} = {y1 y2 y3 y4 y5 y6 y7 y8 y9 y10 y11 y12 y13 y14 y15 y16 y17 y18 y19 y20 y21 y22 y23 y24 y25 y26 y27};

end