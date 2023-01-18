%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Name: Pre-investigation for 53537–503955 pair asteroid 
% Desc: we want to verify the formation of 53537–503955 pair.
%       Back-integrate the current orbit of this pair, and test
%       if they seprated by spin-orbi resonance.
%       53537: https://ssd.jpl.nasa.gov/tools/sbdb_lookup.html#/?sstr=53537&view=OPDA
%       503955: https://ssd.jpl.nasa.gov/tools/sbdb_lookup.html#/?sstr=503955&view=OPDA
% Auth: Hai-Shuo Wang
% Time: 01/10/2023
% Version: InitialPosition v1: This sub-code aims to study 
%           the initial speration position w.r.t shape and ratio of mass.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clc;
clear;
close all;
format LONG;



mu_list = 0.001:0.001:0.5;
a_bA_list = 1.001:0.001:1.5;
for i = 1:length(mu_list)
    for j = 1:length(a_bA_list)
        % Asteroid Parameters
        rA = 15D2; % meter
        q = mu_list(i); % Mass ratio q = rB^3/rA^3
        rB = rA*q^(1/3);
        p = 2100D0; % There is no density information
        
        % There is no shape information
        a_bB = 1.2D0;
        b_cB = 1.2D0;
        a_bA = a_bA_list(j);
        b_cA = 1.2D0;

        R = 2D4;
        
        % Compute Physical Parameters
        [m,A1,A2,A3,IzA,IzB,PO,ThY,re,Unit] = coeff(rA,rB,p,R,a_bB,b_cB,a_bA,b_cA);
        record_re(i,j) = re;
    end
end
absmu = abs(mu_list-0.01); 
index_mu = find(absmu == min(absmu));
Dismu = ['\mu = '  num2str(mu_list(index_mu))];
absab = abs(a_bA_list-1.2);
index_ab = find(absab == min(absab));
Disab = ['a_A/b_A = '  num2str(a_bA_list(index_ab))];
figure
yyaxis left
plot(record_re(index_mu,:)*Unit(1)/1D3,a_bA_list,LineWidth=2,DisplayName=Dismu);
xlabel('Long-equilibiurm Position (km)')
ylabel('a_A/b_A')
yyaxis right
plot(record_re(:,index_ab)*Unit(1)/1D3,mu_list,LineWidth=2,DisplayName=Disab);
ylabel('\mu')
grid on
set(gca,'FontSize',20,'FontWeight','bold')

PlotValue_ab = [1,1.1,1.2,1.3,1.4,1.5];
figure
hold on
for i=1:length(PlotValue_ab)
    absab = abs(a_bA_list-PlotValue_ab(i));
    index_ab = find(absab == min(absab));
    Disab = ['a_A/b_A = '  num2str(a_bA_list(index_ab))];
    plot(mu_list, record_re(:,index_ab)*Unit(1)/1D3,LineWidth=2,DisplayName=Disab);
end
yline(3*rA/1D3,LineWidth=2,Color='red',HandleVisibility='off');
text(0.1,3*rA/1D3,'3r_A',VerticalAlignment='bottom',FontSize=20,Color='red');
yline(8*rA/1D3,LineWidth=2,Color='red',HandleVisibility='off');
text(0.1,8*rA/1D3,'8r_A',VerticalAlignment='bottom',FontSize=20,Color='red');
xlabel('\mu')
ylabel('Long-equilibiurm Position (km)')
grid on
legend
set(gca,'FontSize',20,'FontWeight','bold')




function [m,A1,A2,A3,IzA,IzB,PO,ThY,re,Unit] = coeff(rA,rB,p,R,a_bB,b_cB,a_bA,b_cA)

    cB = rB / (a_bB * b_cB ^ 2D0) ^ (1D0 / 3D0);
    bB = cB * b_cB;
    aB = bB * a_bB;
    cA = rA / (a_bA * b_cA ^ 2D0) ^ (1D0 / 3D0);
    bA = cA * b_cA;
    aA = bA * a_bA;
    UnitL = aA+bA;
    MA = 4D0 / 3D0 * pi * rA ^ 3D0 * p;
    MB = 4D0 / 3D0 * pi * rB ^ 3D0 * p;
    UnitM = MA + MB;
    mu = MB / (MA + MB);
    m = mu * (1 - mu);
    alphaB = rB / UnitL;
    alphaA = rA / UnitL;
    J2B = (aB ^ 2D0 + bB ^ 2D0 - 2D0 * cB ^ 2D0) / (10D0 * rB ^ 2D0);
    J2A = (aA ^ 2D0 + bA ^ 2D0 - 2D0 * cA ^ 2D0) / (10D0 * rA ^ 2D0);
    J22B = (aB ^ 2D0 - bB ^ 2D0) / (20D0 * rB ^ 2D0);
    J22A = (aA ^ 2D0 - bA ^ 2D0) / (20D0 * rA ^ 2D0);
    A1 = (alphaB ^ 2D0 * J2B + alphaA ^ 2D0 * J2A)/2D0;
    A2 = 3D0 * alphaB ^ 2D0 * J22B;
    A3 = 3D0 * alphaA ^ 2D0 * J22A;
    IzB = mu * (aB ^ 2D0 + bB ^ 2D0) / (5D0 * UnitL ^ 2D0);
    IzA = (1 - mu) * (aA ^ 2D0 + bA ^ 2D0) / (5D0 * UnitL ^ 2D0);
    G = 6.67D-11; % N m^2 / kg^2
    UnitT = sqrt(UnitL ^ 3D0 / (MA + MB) / G);
%     UnitT = 1*24*60*60;
    Unit = [UnitL, UnitM, UnitT];
    PO = sqrt(4*pi^2*R^3/(G*(MA+MB)))/Unit(3);

    % Yorp effect
    CY = 0.25D-1;
    Bs = 2D0 / 3D0;
    Fs = 1.0D17;
    au = 1.496D11;
    as = 2.449111518894087*au;
    es = 0.07949126748318633;
    ThY = Bs*Fs*rA/(as^2*sqrt(1-es^2)*MA)*CY;
    ThY = ThY * UnitT^2;

    % Long-term equilibrium position
    fBY = 0.01;
    as = as/Unit(1);
    QAk = 6D5*rA/1D3;
    aAL = aA/Unit(1);
    Fs = Fs*Unit(3)^2/Unit(1)/Unit(2);
    re = (3*mu^2*aAL^5*as^2*sqrt(1-es^2)/(2*QAk*pi*Bs*alphaB^2*fBY*Fs))^(1/7);
end
