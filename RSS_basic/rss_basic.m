%position algorithm for :
%Tx: 3 in (tx,ty,tz)
%Rx: 1 in (locx, locy, locz)
%Text point: 66 in xref.mat
%input: P_LED
%output: (locx, locy, locz) & error
% 66 reference points & 66 estimated points

clear all
close all
clc
%SNR = 50;
% LED coordinates
LED1 = [-0.50 0.60 2.195];
LED2 = [-0.50 0.25 2.195];
LED3 = [0.50  0.25 2.195];
h = 2.195;

tx = [LED1(1), LED2(1), LED3(1)];
ty = [LED1(2), LED2(2), LED3(2)];
tz = [LED1(3), LED2(3), LED3(3)];


%transmitted optical power by individual LED

global P_LED
P_LED = [21, 21.9952, 26];


figure;
plot(tx,ty,'gs',...
    'LineWidth',2,...
    'MarkerSize',10,...
    'MarkerEdgeColor','k',...
    'MarkerFaceColor','y');hold on
%axis([-1 1 0 1]);
%     LEDs={'LED1','LED2','LED3','LED4'};
%     text(x_coor(1,:)+0.15,y_coor(1,:)+0.15,z_coor(1,:),LEDs);

load('xref.mat'); % test points coordinates

for i = 1:length(xref)
    x_ref = xref(i,1);
    y_ref = xref(i,2);
    z_ref = 0.085; % height of PD
    
    plot(x_ref,y_ref,'ob','LineWidth',1.5);hold on
    %axis([-1 1 0 1]);
    %     LEDs={'LED1','LED2','LED3','LED4'};%,'LED5','LED6','LED7','LED8'};
    %     text(x_coor(1,:)+0.15,y_coor(1,:)+0.15,z_coor(1,:),LEDs);
    
    % distances between PD and LED a,b,c
    d_a = sqrt((tx(1)-x_ref)^2 + (ty(1)-y_ref)^2 + (tz(1)-z_ref)^2);
    d_b = sqrt((tx(2)-x_ref)^2 + (ty(2)-y_ref)^2 + (tz(2)-z_ref)^2);
    d_c = sqrt((tx(3)-x_ref)^2 + (ty(3)-y_ref)^2 + (tz(3)-z_ref)^2);
    
    % Calculate Hlos
    theta = 70;     % semi-angle at half power
    m = log10(2)/log10(cosd(theta));    %Lambertian order of emission
    PD_diameter = 7*1e-3; % 7mm
    Adet=pi*(PD_diameter/2)^2; % detector physical area of a PD
    %Adet = 1e-4;    %detector physical area of a PD
    h_a = tz(1) - z_ref;
    h_b = tz(2) - z_ref;
    h_c = tz(3) - z_ref;
    cosphi_a=h_a/d_a;    % angle vector
    cosphi_b=h_b/d_b;
    cosphi_c=h_c/d_c;
%     cosphi_a = 1;
%     cosphi_b = 1;
%     cosphi_c = 1;
    
    Hlos_a = (m+1)*Adet.*cosphi_a.^(m+1)./(2*pi.*d_a.^2);
    Hlos_b = (m+1)*Adet.*cosphi_b.^(m+1)./(2*pi.*d_b.^2);
    Hlos_c = (m+1)*Adet.*cosphi_c.^(m+1)./(2*pi.*d_c.^2);
    
    % received optical powers
    Pr_a = Hlos_a * P_LED(1);
    Pr_b = Hlos_b * P_LED(2);
    Pr_c = Hlos_c * P_LED(3);
    
    % estimated 3 distances
    [de_a, de_b, de_c] = D_estimated(Adet,m,Pr_a, Pr_b, Pr_c, P_LED,h);
    
    % positioning
    A = [2*tx(2)-2*tx(1),2*ty(2)-2*ty(1),2*tz(2)-2*tz(1); 2*tx(3)-2*tx(1),2*ty(3)-2*ty(1),2*tz(3)-2*tz(1)];
    B = [(de_a)^2-(de_b)^2+tx(2)^2+ty(2)^2-tx(1)^2-ty(1)^2-tz(1)^2+tz(2)^2;(de_a)^2-(de_c)^2+tx(3)^2+ty(3)^2-tx(1)^2-ty(1)^2-tz(1)^2+tz(3)^2];
    x = pinv(A' * A) * A' * B;
    
    locx = x(1);
    locy = x(2);
    locz = x(3);
    
    plot(locx,locy,'*r','LineWidth',1.5);hold on
    %axis([-1 1 0 1]);
    error(i) = sqrt((x_ref - locx)^2 + (y_ref - locy)^2 + (z_ref - locz)^2);
    
    %     title('triposition');
end




hold off
grid on
legend('LED','Test point','Estimated point');
xlabel('x (m)');
ylabel('y (m)');
title('Coordinates of estimated points and test points');

figure;
h0 = cdfplot(error(1,:)*100);%hold on
%h1 = cdfplot(error(2,:)*100);hold off
set(h0, 'color', 'g', 'LineStyle', '-.', 'LineWidth', 2);
%set(h1, 'color', 'b', 'LineStyle', '-.', 'LineWidth', 2);
%legend('SNR = 50dB','SNR = 55dB')
xlabel('Positioning Error (cm)');
ylabel('The CDF of the Positioning Error');
