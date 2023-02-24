%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Name: Pre-investigation for 53537–503955 pair asteroid 
% Desc: we want to verify the formation of 53537–503955 pair.
%       Back-integrate the current orbit of this pair, and test
%       if they seprated by spin-orbi resonance.
%       53537: https://ssd.jpl.nasa.gov/tools/sbdb_lookup.html#/?sstr=53537&view=OPDA
%       503955: https://ssd.jpl.nasa.gov/tools/sbdb_lookup.html#/?sstr=503955&view=OPDA
% Auth: Hai-Shuo Wang
% Time: 02/23/2023
% Version 4.3: Parallel compute is employed and plot function is simplified
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

% Rlist = 45D2:5000:15D3;
Rlist = 3*rA:rA:8*rA;

% Set parallel compute parameter
% if length(Rlist)<maxNumCompThreads
%     p=parpool(length(Rlist));
% else
%     p=parpool(maxNumCompThreads);
% end
% p=parpool(7);
parfor i=1:length(Rlist)
    R = Rlist(i);
    % Compute Physical Parameters
    [m,A1,A2,A3,IzA,IzB,PO,ThY,re,Unit] = coeff(rA,rB,p,R,a_bB,b_cB,a_bA,b_cA);
    ThY = ThY*1E5;
    
    % Pr2 = PO/2;
    % Pr1 = 72.74*60*60/Unit(3);
    % SepTime = [565-258,565+902]*1E3*365*24*60*60;
    Pr2 = PO/5.123;
    Pr1 = PO/1.5;
    P2 = PO; % secondary
    SepTime = [0 (2*pi/Pr2-2*pi/Pr1)/ThY]; % 7E+10
    % SepTime = [0 1D6];
    
    % Set Initial Values
    S = R/Unit(1);  % Mutual Distance
    Theta = 0;      % Mutual Rotation Angle
    thetaA = 0;     % Primary Rotation Angle
    thetaB = 0;     % Secondary Rotation Angle
    S1 = 0;         % Radial Velocity
    Theta1 = 2*pi/PO;     % Orbit angular velocity
    thetaA1 = 2*pi/Pr2;    % Primary rotation velocity
    thetaB1 = 2*pi/P2;    % Secondary rotation velocity
    
    X0 = [S;Theta;thetaA;thetaB;S1;Theta1;thetaA1;thetaB1];
    tspan = [0 SepTime(2)*1.2]; % Integration Time Span
    coef = [m,A1,A2,A3,IzA,IzB,ThY];    % Parameters for Right Function
    
    % Because the integration time is too long, we cant recorde every data that
    % will be heavy burden for memory. 
    Num = 10000; % The integration time span will be sperated 'Num' arcs
    Tstep = (tspan(2) - tspan(1)) / Num;
    X = X0;
    [TT,yy,Energy,Momentum] = PropagationLoop(Num, Tstep, tspan, X, coef);

    Record_y{i} = yy;
    Record_T{i} = TT;
    Record_E{i} = Energy;
    Record_M{i} = Momentum;
end
% delete parallel compute
% delete(p);

save('Preinvestigation_v4_3');
% Plot Results
figure
hold on
for i=1:length(Rlist)
    DispName = ['Initial Position: ' num2str(Rlist(i)/rA) ' rA'];
    scatter(Record_T{i}*Unit(3)/(365*24*60*60), ...
        Record_y{i}(7,:)./Record_y{i}(6,:),...
        Marker=".",DisplayName=DispName);
end
yline(1.5,LineWidth=2,Color='red',HandleVisibility='off');
yline(2,LineWidth=2,Color='red',HandleVisibility='off');
xlabel('time (yr)')
ylabel('\omega_A/n')
ylim([0 5]);
grid on
set(gca,'FontSize',20,'FontWeight','bold')

%% Subfunction

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

function [Record_T,Record_y,Energy,Momentum] = PropagationLoop(Num, Tstep, tspan, X, coef)
    opts = odeset('RelTol',1e-13,'AbsTol',1e-13); % Setting Tolerance
    m = coef(1);
    A1 = coef(2);
    A2 = coef(3);
    A3 = coef(4);
    IzA = coef(5);
    IzB = coef(6);
    for i = 1:Num
        fprintf('Loop index: %i/%i \n',i,Num);
        % Integrate one arc
        tspani = [tspan(1)+Tstep*(i-1) tspan(1)+Tstep*i];
        [Ti,yi] = ode45(@(t,X) YHC4BT(t,X,coef),tspani,X,opts);
        yi = yi';
        X = yi(:,end);
        
        % Compute Energy and total Momentum
        [S,Theta,thetaA,thetaB,S1,Theta1,thetaA1,thetaB1] = X2S(yi(:,end));
        
        Energy(i) = m / 2D0 * (S1 ^ 2D0 + S ^ 2D0 * Theta1 ^ 2D0)...
            + 1D0 / 2D0 * IzB * thetaB1 ^ 2D0...
            + 1D0 / 2D0 * IzA * thetaA1 ^ 2D0 - m / S...
            - m / (2D0 * S ^ 3D0) * (A1 + A2 * cos(2D0 * Theta - 2D0 * thetaB)...
            + A3 * cos(2D0 * Theta - 2D0 * thetaA));
        Momentum(i) = m * S ^ 2D0 * Theta1 + IzA * thetaA1 + IzB * thetaB1;
    
        % Record data
        Record_y(:,i) = yi(:,end);
        Record_T(i) = Ti(end);
    end

end

function dydt = YHC4BT(t,X,coef)
    % Right Function for 53537–503955 pair asteroid

    m = coef(1);
    A1 = coef(2);
    A2 = coef(3);
    A3 = coef(4);
    IzA = coef(5);
    IzB = coef(6);
    ThY = coef(7);
    
    S = X(1);
    Theta = X(2);
    thetaA = X(3);
    thetaB = X(4);
    S1 = X(5);
    Theta1 = X(6);
    thetaA1 = X(7);
    thetaB1 = X(8);

    Y(1) = S1;
    Y(2) = Theta1;
    Y(3) = thetaA1;
    Y(4) = thetaB1;
    Y(5) = S * Theta1^2 - 1D0 / S^2D0...
        - (3D0 / 2D0) * (A1 + A2 * cos(2D0 * Theta - 2D0 * thetaB)...
        + A3 * cos(2D0 * Theta - 2D0 * thetaA)) / S ^ 4D0;
    Y(6) = -A2 * sin(2D0 * Theta - 2D0 * thetaB) / S ^ 5D0...
        - A3 * sin(2D0 * Theta - 2D0 * thetaA) / S ^ 5D0...
        - 2D0 * S1 * Theta1 / S;
    Y(7) = m * A3 * sin(2D0 * Theta - 2D0 * thetaA) / (IzA * S ^ 3D0)-ThY;
    Y(8) = m * A2 * sin(2D0 * Theta - 2D0 * thetaB) / (IzB * S ^ 3D0);
    
    dydt = Y';
end

function [S,Theta,thetaA,thetaB,S1,Theta1,thetaA1,thetaB1] = X2S(X)
    S = X(1);
    Theta = X(2);
    thetaA = X(3);
    thetaB = X(4);
    S1 = X(5);
    Theta1 = X(6);
    thetaA1 = X(7);
    thetaB1 = X(8);
end

function Plot_func(Record_T,Record_y,Energy,Momentum)
    figure
    yyaxis left
    plot(Record_T,Energy,LineWidth=2);
    xlabel('time')
    ylabel('Energy')
    yyaxis right
    plot(Record_T,Momentum,LineWidth=2);
    ylabel('Momentum')
    grid on
    set(gca,'FontSize',20,'FontWeight','bold')
    
    figure
    yyaxis left
    plot(Record_T,Record_y(7,:),LineWidth=2);
    xlabel('time')
    ylabel('\omega_A')
    yyaxis right
    plot(Record_T,Record_y(6,:),LineWidth=2);
    ylabel('n')
    grid on
    set(gca,'FontSize',20,'FontWeight','bold')
    
    figure
    scatter(Record_T, Record_y(7,:)./Record_y(6,:),'.','black');
    xlabel('time')
    ylabel('\omega_A/n')
    grid on
    set(gca,'FontSize',20,'FontWeight','bold')
    
    figure
    scatter(Record_T, Record_y(1,:),'.','black');
    xlabel('time')
    ylabel('S')
    grid on
    set(gca,'FontSize',20,'FontWeight','bold')
    
    figure
    scatter(Record_y(1,:).*cos(Record_y(2,:)), Record_y(1,:).*sin(Record_y(2,:)),'.','black');
    xlabel('x')
    ylabel('y')
    grid on
    set(gca,'FontSize',20,'FontWeight','bold')
end

