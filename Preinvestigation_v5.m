%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Name: Pre-investigation for 53537–503955 pair asteroid 
% Desc: we want to verify the formation of 53537–503955 pair.
%       Back-integrate the current orbit of this pair, and test
%       if they seprated by spin-orbi resonance.
%       53537: https://ssd.jpl.nasa.gov/tools/sbdb_lookup.html#/?sstr=53537&view=OPDA
%       503955: https://ssd.jpl.nasa.gov/tools/sbdb_lookup.html#/?sstr=503955&view=OPDA
% Auth: Hai-Shuo Wang
% Time: 01/10/2023
% Version 5: A totally new method to investigate this pair.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clc;
clear;
close all;
format LONG;

% Asteroid Parameters
rA = 15D2; % meter
q = 0.01; % Mass ratio q = rB^3/rA^3
rB = rA*q^(1/3);
p = 2100D0; % There is no density information

% There is no shape information
a_bB = 1.2D0;
b_cB = 1.2D0;
a_bA = 1.2D0;
b_cA = 1.2D0;
cB = rB / (a_bB * b_cB ^ 2D0) ^ (1D0 / 3D0);
bB = cB * b_cB;
aB = bB * a_bB;
cA = rA / (a_bA * b_cA ^ 2D0) ^ (1D0 / 3D0);
bA = cA * b_cA;
aA = bA * a_bA;
% UnitL = R;
MA = 4D0 / 3D0 * pi * rA ^ 3D0 * p;
MB = 4D0 / 3D0 * pi * rB ^ 3D0 * p;
G = 6.67D-11; % N m^2 / kg^2
% Yorp effect
CY = 0.25D-1;
Bs = 2D0 / 3D0;
Fs = 1.0D17; % kg m/s^2
au = 1.496D11;
as = 2.449111518894087*au;
es = 0.07949126748318633;
ThY = Bs*Fs*rA/(as^2*sqrt(1-es^2)*MA)*CY;
% ThY = ThY * UnitT^2;


SepTime = [565-258,565+902]*1E3*365*24*60*60;
Pr1 = 72.74*60*60;
R = 2*rA:7*rA;
j = 1;
T2(j,:) = 1./ThY./j*sqrt(G.*MA./R.^3) - 2.*pi./Pr1./ThY;
T2(j,:) = T2(j,:)./3600./24./365./1000;

j = 2;
T2(j,:) = 1./ThY./j*sqrt(G.*MA./R.^3) - 2.*pi./Pr1./ThY;
T2(j,:) = T2(j,:)./3600./24./365./1000;

j = 3;
T2(j,:) = 1./ThY./j*sqrt(G.*MA./R.^3) - 2.*pi./Pr1./ThY;
T2(j,:) = T2(j,:)./3600./24./365./1000;

j = 1.5;
T2(4,:) = 1./ThY./j*sqrt(G.*MA./R.^3) - 2.*pi./Pr1./ThY;
T2(4,:) = T2(4,:)./3600./24./365./1000;

figure
hold on
% plot(R/rA,T2(1,:),LineWidth=2,DisplayName='p/q=1')
plot(R/rA,T2(2,:),LineWidth=2,DisplayName='p/q=2')
plot(R/rA,T2(3,:),LineWidth=2,DisplayName='p/q=3')
plot(R/rA,T2(4,:),LineWidth=2,DisplayName='p/q=3/2')
yline(565-258,Color='blue',LineStyle='-.',LineWidth=2,DisplayName='307 kyr')
yline(565+902,Color='red',LineStyle='-.',LineWidth=2,DisplayName='1467 kyr')
xlabel('Seperate Distance (unit: Primary Radii)')
ylabel('Seperate Time (kyr)')
grid on
set(gca,'FontSize',20,'FontWeight','bold')

