% LKA-Simulation: Simulink
% 
% Associated Simulink-Models: 'LKA.mdl'
% 
% Subject: lka
% Author: georgnoname
% Date: 08.11.2012 - 29.04.2013

clc
clear
% close all

% load Simulink model
mdlName = 'LKA';
load_system(mdlName);


%% LKA-controller design

paramFile_SingleTrackMdl_BMW5

% nicht steuerbar bei
vx0 = sqrt(-cs_front*(Izz-l_front*l_rear*m)*(l_front+l_rear)/(l_front^2*m^2));

% LKA-controller design: vehicle parameter
stringCvehicle = 'paramFile_SingleTrackMdl_BMW5';

% LKA-controller design: longitudinal velocity vx [m/s]
% vxC = vx0;
vx_Ctrl = 20;

% LKA-controller design: look-ahead distance lad [m]
LAD_Ctrl = 10;

% LKA-controler design: calculate controller with upper parameters
contr.t = lkaController_t(stringCvehicle,vx_Ctrl,LAD_Ctrl);
addpath('algSynth');
% contr.s = lkaController_s(stringCvehicle,vx_Ctrl,LAD_Ctrl);


%% Simulation parameter

% simulation: vehicle parameter
stringSim = 'paramFile_SingleTrackMdl_BMW5';
pin.VehicleModel.parameterFile = stringSim;
pin.VehicleModel.parameter = paramFile2Struct(stringSim);

% simulation: longitudinal velocity vx [m/s]
pin.vx = vx_Ctrl;

% simulation: look-ahead distance lad [m]
pin.LAD = LAD_Ctrl;

% simulation: vehicle model state-space model
% pin.VehicleModel.ssMdl = ssMdl_SingleTrack('stvis',stringSim,pin.vx,pin.lad);

% load desired path
% [pin.traj,trajErr] = lka_trajectory_07(0,0.05);
pin.traj = ...
	LkPathStraight(0.05,50,0) + ...
	LkPathCircle(0.05,-pi/2,0,200);
pin.traj_sd = pin.traj.pathData;

% simulation: interval of integration
tend = pin.traj.length/pin.vx - 2*pin.LAD/pin.vx;

% vehicle model: initial condition
pin.VehicleModel.initialValue.sy = 0;
pin.VehicleModel.initialValue.vy = 0;
pin.VehicleModel.initialValue.psi = 0;
pin.VehicleModel.initialValue.psiDot = 0;
pin.VehicleModel.initialValue.x_g = 0;
pin.VehicleModel.initialValue.y_g = 0;
pin.VehicleModel.initialValue.s_g = 0;
pin.VehicleModel.initialValue.yL = 0;
pin.VehicleModel.initialValue.epsL = 0;


%% set solver parameter

% für Sollbahn mit sprungf. Krümmungsverlauf (kapL wird richtig ermittelt)
set_param(mdlName,'SolverType','fixed-step')
set_param(mdlName,'solver','ode3')
set_param(mdlName,'fixedstep','0.01')

% sonst
set_param(mdlName,'SolverType','Variable-step')
set_param(mdlName,'Solver','ode23')
set_param(mdlName,'MaxStep','auto')


%% run simulation

% sim
tic
sim(mdlName);
toc

% Simulationsprogramm
soli.simProg = 'Simulink';

% Zeitstempel
soli.simDate = datestr(now);

% input parameter
sol.simIn = pin;

% assign simulink to workspace output
sol.simOut = simout;

%% post-processing

% sol.(lbl) = lkaPostProcessing(soli.(lbl),pin,[],simoutState,simoutSensor,simoutSteerAngle);


%% plot

figure;

subplot(3,1,1);
plot(sol.simOut.LaneSensor.curvature_LAD_Sensor);
grid on

subplot(3,1,2);
plot(sol.simOut.Control);
grid on

subplot(3,1,3);
plot(sol.simOut.LaneSensor.lateralOff_LAD_Sensor);
hold all
plot(sol.simOut.LaneSensor.poly_lateralOff_LAD_Sensor.Time,...
	squeeze(sol.simOut.LaneSensor.poly_lateralOff_LAD_Sensor.Data(:,end,:)));
grid on
legend(['LAD = ',num2str(pin.LAD)],'LAD = 0');


%% export

% filename = ['lkaSim_ohneLenkung','_traj',num2str(pin.traj.ID.nbr),'_',lbl];
% 
% % save as mat-file
% save(filename,'sol');
% 
% % save as ascii-file
% saveFigure_txt(filename,...
%     300,...
%     'xtraj',sol.(lbl).simIn.traj.x,...
%     'ytraj',sol.(lbl).simIn.traj.y,...
%     't',sol.(lbl).simOut.t,...
%     'xg',sol.(lbl).simOut.vehicleState.x_g,...
%     'yg',sol.(lbl).simOut.vehicleState.y_g,...
%     'yL_Sensor',sol.(lbl).simOut.lkaSensor.yL,...
%     'epsL_Sensor',sol.(lbl).simOut.lkaSensor.epsL,...
%     'epsL_exactMdl',sol.(lbl).simOut.vehicleState.epsL_exactMdl,...
%     'u',rad2deg(sol.(lbl).simOut.controlInp.steeringAngle));

