% test_yL_epsL
% 
% Associated Simulink-Models: 'test_yL_epsL.mdl'
% 
% Subject: 
% Author: georgnoname
% Date: 18.03.2013 - 02.05.2013

clc
clear all
% close all

% run init file
initFile_;

% load Simulink model
mdlName = 'test_yL_epsL';
load_system(mdlName);


%% LKA-controller design

% LKA-controller design: vehicle parameter
stringC = 'paramFile_SingleTrackMdl_BMW5';

% LKA-controller design: longitudinal velocity vx [m/s]
vxC = 30;

% LKA-controller design: look-ahead distance lad [m]
ladC = 10;

% LKA-controler design: calculate controller with upper parameters
contr.t = lkaController_t(stringC,vxC,ladC);


%% Simulation parameter

% simulation: vehicle parameter
stringSim = 'paramFile_SingleTrackMdl_BMW5';
pin.VehicleModel.parameterFile = stringSim;
pin.VehicleModel.parameter = loadParameter(stringSim,'about');

% simulation: longitudinal velocity vx [m/s]
pin.vx = vxC;

% simulation: look-ahead distance lad [m]
pin.lad = ladC;

% load intended trajectory
[pin.traj,trajErr] = lka_trajectory_01(0,0.05);

% simulation: interval of integration
% tend = pin.traj.s(end)/pin.vx - 2*pin.lad/pin.vx;
tend = (max(pin.traj.x)-2*pin.lad)/pin.vx;
tend = floor(tend);

% reduce 'MaxStep' or enlarge 'delta' of lka_trajectory_xx(sw,delta) if
% simulation fails
simset('Solver','ode45','MaxStep',0.3);

% vehicle model: initial condition
pin.VehicleModel.initialValue.sy = 0;
pin.VehicleModel.initialValue.vy = 0;
pin.VehicleModel.initialValue.psi = deg2rad(0);
pin.VehicleModel.initialValue.psiDot = 0;
pin.VehicleModel.initialValue.x_g = 0;
pin.VehicleModel.initialValue.y_g = 0;
pin.VehicleModel.initialValue.s_g = 0;
pin.VehicleModel.initialValue.yL = 0;
pin.VehicleModel.initialValue.epsL = 0;

% deltat = simget(mdlName,'FixedStep');

% Lenkwinkelvorgabe für Trajektorie lka_trajectory_08 (Klothoide)
% t0 = 5;
% t1 = 8;
% t2 = 10;
% tend = tend+6;
% delta = @(u) 0.005*(u-t0)*(u>t0) - 0.001*(u-t1)*(u>t1) - 0.005*(u-t2)*(u>t2);

% Lenkwinkelvorgabe für Trajektorie lka_trajectory_01 (Viertelkreis), vx=45m/s
t0 = 2.5;
t1 = 11;
t2 = 20;
tend = 20;
delta = @(u) 0.0025*(u-t0)*(u>t0) - 0.005*(u-t1)*(u>t1) + 0.0025*(u-t2)*(u>t2);

% Lenkwinkelvorgabe für Trajektorie lka_trajectory_01 (Viertelkreis), vx=30m/s
t0 = 4;
t1 = 16;
tend = 30;
delta = @(u) 0.0015*(u-t0)*(u>t0) - 0.0035*(u-t1)*(u>t1);


%% set solver parameter

% für Sollbahn mit sprungf. Krümmungsverlauf (kapL wird richtig ermittelt)
set_param(mdlName,'SolverType','Fixed-step')
set_param(mdlName,'Solver','ode3')
set_param(mdlName,'FixedStep','0.01')

% sonst
set_param(mdlName,'SolverType','Variable-step')
set_param(mdlName,'Solver','ode23')
set_param(mdlName,'MaxStep','auto')


%% run simulation

% sim
tic
sim('test_yL_epsL');
toc

% get lka-controller-string used in simulation
% lbl = 'none';

% Simulationsprogramm
sol.simProg = 'Simulink';

% Zeitstempel
sol.simDate = datestr(now);

% input parameter
sol.simIn = pin;


%% post-processing

% time
sol.simOut.t = simoutState.time;

% Vehicle State
sol.simOut.vehicleState.sy = simoutState.signals.values(:,1);
sol.simOut.vehicleState.syDot = simoutState.signals.values(:,2);
sol.simOut.vehicleState.psi = simoutState.signals.values(:,3);
sol.simOut.vehicleState.psiDot = simoutState.signals.values(:,4);
sol.simOut.vehicleState.x_g = simoutState.signals.values(:,5);
sol.simOut.vehicleState.y_g = simoutState.signals.values(:,6);
sol.simOut.vehicleState.s_g = simoutState.signals.values(:,7);
sol.simOut.vehicleState.yL_linMdl = simoutState.signals.values(:,8);
sol.simOut.vehicleState.epsL_linMdl = simoutState.signals.values(:,9);
sol.simOut.vehicleState.yL_exactMdl = simoutState.signals.values(:,10);
sol.simOut.vehicleState.epsL_exactMdl = simoutState.signals.values(:,11);

% steering wheel function handle
sol.simIn.controlInp.steeringAngle.fh = delta;
sol.simIn.controlInp.steeringAngle.tend = sol.simOut.t(end);

% Vehicle Video Sensor
sol.simOut.lkaSensor.yL = simoutSensor.signals.values(:,1);
sol.simOut.lkaSensor.epsL = simoutSensor.signals.values(:,2);
sol.simOut.lkaSensor.kapL = simoutSensor.signals.values(:,3);

% steering wheel input
sol.simOut.controlInp.steeringAngle = simoutSteerAngle.signals.values;

sol.simPost = calc_relPos(sol,0);


%% plot

figure;

% yL
subplot(2,1,1);
plot(sol.simOut.t,[...
    sol.simOut.lkaSensor.yL,...
    sol.simOut.vehicleState.yL_linMdl,...
    sol.simOut.vehicleState.yL_exactMdl]);
grid on
title('y_L')

% epsL
subplot(2,1,2);
plot(sol.simOut.t,[...
    sol.simOut.lkaSensor.epsL,...
    sol.simOut.vehicleState.epsL_linMdl,...
    sol.simOut.vehicleState.epsL_exactMdl]);
grid on
title('eps_L')

legend('Sensor','lineares Modell','exaktes Modell')


%% export

% % filename
% filename = ['vglRelpos_exakt_modell','_traj',num2str(pin.traj.ID.nbr)];
% filename = [filename,'_vx',num2str(vxC)];
% 
% % save as mat-file
% save(filename,'sol');
% 
% % save as ascii-file
% saveFigure_txt(filename,...
%     97,...
%     'xtraj',sol.simIn.traj.x,...
%     'ytraj',sol.simIn.traj.y,...
%     't',sol.simOut.t,...
%     'xg',sol.simOut.vehicleState.x_g,...
%     'yg',sol.simOut.vehicleState.y_g,...
%     'yL_Sensor',sol.simOut.lkaSensor.yL,...
%     'yL_linMdl',sol.simOut.vehicleState.yL_linMdl,...
%     'yL_exactMdl',sol.simOut.vehicleState.yL_exactMdl,...
%     'epsL_Sensor',sol.simOut.lkaSensor.epsL,...
%     'epsL_linMdl',sol.simOut.vehicleState.epsL_linMdl,...
%     'epsL_exactMdl',sol.simOut.vehicleState.epsL_exactMdl,...
%     'u',rad2deg(sol.simOut.controlInp.steeringAngle));

