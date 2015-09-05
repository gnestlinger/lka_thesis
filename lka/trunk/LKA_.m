% LKA-Simulation: Simulink
% 
% Associated Simulink-Models: 'LKA.mdl'
% 
% Subject: lka
% Author: georgnoname
% Date: 08.11.2012 - 29.04.2013

clc
clear all
% close all

% run init file
initFile_;

% load Simulink model
mdlName = 'LKA';
load_system(mdlName);


%% LKA-controller design

% LKA-controller design: vehicle parameter
stringCvehicle = 'paramFile_SingleTrackMdl_BMW5';

% LKA-controller design: longitudinal velocity vx [m/s]
vxC = 20;

% LKA-controller design: look-ahead distance lad [m]
ladC = 10;

% LKA-controler design: calculate controller with upper parameters
contr.t = lkaController_t(stringCvehicle,vxC,ladC);
addpath('algSynth');
contr.s = lkaController_s(stringCvehicle,vxC,ladC);


%% Simulation parameter

% Configureable Subsystem: select controller
% set_param([mdlName,'/LKA-Controller'],'BlockChoice','LQR');
% set_param([mdlName,'/LKA-Controller'],'BlockChoice','LQR_2int');
lbl = get_param([mdlName,'/LKA-Controller'],'BlockChoice');

% store selected control strategy in parameter-in structure pin
try
    pin.Controller.lka = contr.t.(lbl);
catch exception
    pin.Controller.lka = contr.s.(lbl);
end

% simulation: vehicle parameter
stringSim = 'paramFile_SingleTrackMdl_BMW5';
pin.VehicleModel.parameterFile = stringSim;
pin.VehicleModel.parameter = loadParameter(stringSim,'about');

% simulation: longitudinal velocity vx [m/s]
pin.vx = vxC;

% simulation: look-ahead distance lad [m]
pin.lad = ladC;

% simulation: vehicle model state-space model
% pin.VehicleModel.ssMdl = ssMdl_SingleTrack('stvis',stringSim,pin.vx,pin.lad);

% add folder 'trajectories' to search path and load intended trajectory
addpath('trajectories');
[pin.traj,trajErr] = lka_trajectory_07(0,0.05);

% simulation: interval of integration
tend = pin.traj.s(end)/pin.vx - 2*pin.lad/pin.vx;

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
% set_param(mdlName,'SolverType','Fixed-step')
% set_param(mdlName,'Solver','ode3')
% set_param(mdlName,'FixedStep','0.01')

% sonst
set_param(mdlName,'SolverType','Variable-step')
set_param(mdlName,'Solver','ode23')
set_param(mdlName,'MaxStep','auto')


%% run simulation

% sim
tic
sim(mdlName);
toc

% get lka-controller-string used in simulation
lbl = get_param([mdlName,'/LKA-Controller'],'BlockChoice');

% Simulationsprogramm
soli.(lbl).simProg = 'Simulink';

% Zeitstempel
soli.(lbl).simDate = datestr(now);


%% post-processing

sol.(lbl) = lkaPostProcessing(soli.(lbl),pin,[],simoutState,simoutSensor,simoutSteerAngle);


%% plot

% figure
% lkaPlot(sol,0,'traj');
lkaPlot(sol,0,'add');
% lkaPlot(sol,pin.lad,'add');
% lkaPlot(sol,0);


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

