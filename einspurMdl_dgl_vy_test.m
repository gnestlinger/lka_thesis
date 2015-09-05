% LKA Reglerentwurf: Testfunktion für Einspurmodell-DGL
% 
% Subject: lka
% Author: georgnoname
% Date: 24.09.2012 - 25.02.2013

clc
% clear all

% run init file
initFile_;

% Zustandsvektor: x = [sy_V; vy_V; psi; psiDot; sx_0; sy_0]
% initial condition
xInit = [0;0;0;0;0;0;0];

% interval of integration
tend = 35.2;
tspan = [0,tend];

% input u = delta
k = deg2rad(1); u = @(t,x) k*(t>=5); % [k] = °
% u = @(t,x) interp1(sol.LQR_2int.x,sol.LQR_2int.addData.u.value,t);

% simulation: vehicle parameter
stringSim = 'paramFile_SingleTrackMdl_BMW5';
pin.VehicleModel.singleTrack.parameterFile = stringSim;
pin.VehicleModel.singleTrack.parameter = loadParameter(stringSim,'about');

% Längsgeschwindigkeit [m/s]
pin.vx = 25;  

% solver
options = odeset('MaxStep',0.3);
sol1 = ode45(@(t,x) einspurMdl_dgl_vy(t,x,u,pin),tspan,xInit,options);
sol1.simparameter = pin;


%% plot

figure
subplot(2,1,1)
plot( sol1.y(5,:),sol1.y(6,:) ); grid on;
xlabel('x_{0} [m]'); ylabel('y_{0} [m]');
title(['v_x = ',num2str(sol1.simparameter.vx),' m/s']);
axis equal

subplot(2,1,2)
p = plot( sol1.x,[sol1.y([2,4],:);u(sol1.x,0)]); grid on;
set(p(end),'LineStyle','--')
legend('v_y','Giergeschw.','u');
xlabel('t [s]');

