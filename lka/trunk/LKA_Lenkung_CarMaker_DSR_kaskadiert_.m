% LKA-Simulation mit CarMaker-Lenkmodell: Simulink
% 
% Associated Simulink-Models: 
%   .) Lenkung_CarMaker_DSR.mdl
%   .) LKA_Lenkung_CarMaker_DSR_kaskadiert.mdl
% 
% Subject: lka
% Author: georgnoname
% Date: 11.12.2012 - 26.04.2013

clc
clear all
% close all

% run init file
initFile_;

% load Simulink model
mdlName = 'LKA_Lenkung_CarMaker_DSR_kaskadiert';
load_system(mdlName);


%% steering model controller design

% parameter file
stringC_steer = 'paramFile_SteeringMdl_CarMaker_DSR';

% state-space model
sysSteeringModel = ssMdl_EPS('CarMakerDSR',stringC_steer);

% transfer function
P = tf(sysSteeringModel);
P = P(1);

% Wunsch: Überschwingen (Mp) und Anstiegszeit (tr)
Mp = 1.0;
tr = 0.01;

% Tabelle (vgl. RT Horn/Dourdoumas S. 347)
xi = interp1([1.00;1.02;1.05;1.10;1.16;1.25;1.37;1.53;1.73],...
	[0.9,0.8,0.7,0.6,0.5,0.4,0.3,0.2,0.1],Mp,'spline','extrap');
wntr = interp1([1.00;1.02;1.05;1.10;1.16;1.25;1.37;1.53;1.73],...
    [2.85,2.55,2.28,2.04,1.84,1.66,1.50,1.37,1.26],Mp,'spline','extrap');

% Soll-Führungsübertragungsfunktion (System mit dominantem Polpaar)
s = tf('s');
Tsoll = (wntr/tr)^2/(s^2 + 2*xi*wntr/tr*s + (wntr/tr)^2);

% Ordnungserhöhung von Tsoll
ws = s+10; % Hurwitz-Polynom!
Tsoll = Tsoll*ws/ws;

% steering model: controller design (algebraische Synthese)
pin.Controller.steerSys.desc = 'Algebraische Synthese';
[R,V] = algsynth(P,Tsoll);
pin.Controller.steerSys.R.tf = R;
pin.Controller.steerSys.R.ss = ssNormForm(R,'obs');
pin.Controller.steerSys.V.tf = V;
pin.Controller.steerSys.V.ss = ssNormForm(V,'obs');



%% LKA-controller design

% LKA-controller design: vehicle parameter
stringC = 'paramFile_SingleTrackMdl_BMW5';

% LKA-controller design: longitudinal velocity vx [m/s]
vxC = 25;

% LKA-controller design: look-ahead distance lad [m]
ladC = 10;

% LKA-controler design: calculate controler with upper parameters
contr.t = lkaController_t(stringC,vxC,ladC);
contr.s = lkaController_s(stringC,vxC,ladC);



%% Simulation parameter

% Configureable Subsystem: choose controller
set_param([mdlName,'/LKA-Controller'],'BlockChoice','LQR_2int');
lbl = get_param([mdlName,'/LKA-Controller'],'BlockChoice');

% store selected control strategy in parameter-in structure pin
try
    pin.Controller.lka = contr.t.(lbl);
catch exception
    pin.Controller.lka = contr.s.(lbl);
end

% simulation: steering parameter
stringSim = 'paramFile_SteeringMdl_CarMaker_DSR';
pin.SteeringModel.parameterFile = stringSim;
pin.SteeringModel.parameter = loadParameter(stringSim,'about');

% simulation: vehicle parameter
stringSim = 'paramFile_SingleTrackMdl_BMW5';
pin.VehicleModel.parameterFile = stringSim;
pin.VehicleModel.parameter = loadParameter(stringSim,'about');

% simulation: longitudinal velocity vx [m/s]
pin.vx = vxC;

% simulation: look-ahead-distance lad [m]
pin.lad = ladC;

% load intended trajectory
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

% load Lenkuebersetzung_Demo_BMW_5_spline;


%% set solver parameter

% set_param(mdlName,'SolverType','Fixed-step')
% set_param(mdlName,'Solver','ode3')
% set_param(mdlName,'FixedStep','0.01')

set_param(mdlName,'SolverType','Variable-step')
set_param(mdlName,'Solver','ode23')
set_param(mdlName,'MaxStep','auto')


%% run simulation

% run
tic
sim('LKA_Lenkung_CarMaker_DSR_kaskadiert');
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
lkaPlot(sol,0,'traj');
% lkaPlot(sol,0,'add');
% lkaPlot(sol,pin.lad,'add');
% lkaPlot(sol,0);


