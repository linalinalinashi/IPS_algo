clc
clear all
close all
%---------------------ENTER PARAMETERS------------------------%
space = 1;
% Distance between tx and rx ( Meter )
heightLED = 1.48;

lx=4; ly=4; lz=2.15; % room dimension in meter

TP1=[-lx/2 lx/2 heightLED]; % transmitter position
TP2=[lx/2 -lx/2 heightLED]; % transmitter position
TP3=[-lx/2 -lx/2 heightLED]; % transmitter position
TP4=[lx/2 lx/2 heightLED]; % transmitter position
% Transmitter Semi-angle, angle of irradiance in half (Radian)
theta=70; % semi-angle at half power
phi = (theta*pi)/180;
% Speed of Light
c = 300E6; 
% Time
t = 0:0.01:4;
%---------------------END OF PARAMETERS-----------------------%

% 3D Meshgrid X-axis and Y-axis %
radius = heightLED * tan(phi);
[X,Y] = meshgrid(-4:0.2:4); 

xydist = sqrt((X+TP1(1)).^2 + (Y+TP1(2)).^2);
hdist = sqrt(xydist.^2 + heightLED.^2);
% Incidence angles of receiver according to X-Y axis % 
incidence = atand(xydist.* heightLED ^(-1));
[P, PO, Z, impulset, impulsetd, impulsef, impulsefd]=RxSNR(incidence,hdist,t,phi); % SNR in dB at each X-Y location %
 % Received Power in mW at each X-Y location %
P1=P;

xydist = sqrt((X+TP2(1)).^2 + (Y-TP2(2)).^2);
hdist = sqrt(xydist.^2 + heightLED.^2); 
incidence = atand(xydist.* heightLED ^(-1));
[P, PO, Z, impulset, impulsetd, impulsef, impulsefd]=RxSNR(incidence,hdist,t,phi); % SNR in dB at each X-Y location % 
P2=P;

xydist = sqrt((X-TP3(1)).^2 + (Y+TP3(2)).^2);
hdist = sqrt(xydist.^2 + heightLED.^2); 
incidence = atand(xydist.* heightLED ^(-1));
[P, PO, Z, impulset, impulsetd, impulsef, impulsefd]=RxSNR(incidence,hdist,t,phi); % SNR in dB at each X-Y location % 
P3=P; 

xydist = sqrt((X-TP4(1)).^2 + (Y-TP4(2)).^2);
hdist = sqrt(xydist.^2 + heightLED.^2); 
incidence = atand(xydist.* heightLED ^(-1));
[P, PO, Z, impulset, impulsetd, impulsef, impulsefd]=RxSNR(incidence,hdist,t,phi); % SNR in dB at each X-Y location % 
P4=P;

P=P1+P2+P3+P4;


%%%%%%%%%%%%%%%%%%%%%Plot Circle%%%%%%%%%%%%%%%%%%%%%%%%%%%
[X1,Y1]=circle([TP1(1),TP1(2)],radius,150);
figure;
plot(X1,Y1,'+')
xlabel('Length of surface [m]')
ylabel('Breadth of surface [m]')
title('Top view of LED light coverage')
hold on

[X1,Y1]=circle([TP2(1),TP2(2)],radius,150);
plot(X1,Y1,'+')
hold on

[X1,Y1]=circle([TP3(1),TP3(2)],radius,150);
plot(X1,Y1,'+')
hold on

[X1,Y1]=circle([TP4(1),TP4(2)],radius,150);
plot(X1,Y1,'+')
hold off
%%%%%%%%%%%%%%%%%%%%%End Circle%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%3D diagram for Received Power WITH LENS%%%%%%%%%%%%%%%%%
figure
mesh(X,Y,P,'EdgeColor','black')
xlabel('Length of room [m]')
ylabel('Width of room [m]')
zlabel('Received Power in (W)')
title('3D Plot for Room Receivered Power Distribution with Lens')
