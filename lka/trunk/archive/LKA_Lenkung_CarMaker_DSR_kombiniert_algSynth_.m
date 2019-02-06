% LKA-Simulation mit CarMaker-Lenkmodell: Simulink
% 
% Associated Simulink-Models: 
%   .) LKA_Lenkung_CarMaker_DSR_kombiniert.mdl
% 
% Subject: LKA
% $Author$
% $LastChangedDate$
% $Revision$
clc
clear all
% close all

% load Simulink model
mdlName = 'LKA_Lenkung_CarMaker_DSR_kombiniert_algSynth';
load_system(mdlName);


%% LKA-controller design

% LKA-controller design: control plant parameter
stringCvehicle = 'paramFile_SingleTrackMdl_BMW5';
stringCsteer = 'paramFile_SteeringMdl_CarMaker_DSR';

% LKA-controller design: longitudinal velocity vx [m/s]
vxC = 20;

% LKA-controller design: look-ahead distance lad [m]
ladC = 10;

% load state-space model
sysSingleTrackVisDSR = ...
    ssMdl_SingleTrack('stvisdsr',stringCvehicle,vxC,ladC,stringCsteer);

% Streckenübertragungsfunktion
P = tf(sysSingleTrackVisDSR);
P = P(3,1);

% Wunsch: Überschwingen (Mp) und Anstiegszeit (tr)
Mp = 1;
tr = 2;

% Tabelle (vgl. RT Horn/Dourdoumas S. 347)
xi = interp1([1.00;1.02;1.05;1.10;1.16;1.25;1.37;1.53;1.73],...
	[0.9,0.8,0.7,0.6,0.5,0.4,0.3,0.2,0.1],Mp,'spline','extrap');
wntr = interp1([1.00;1.02;1.05;1.10;1.16;1.25;1.37;1.53;1.73],...
    [2.85,2.55,2.28,2.04,1.84,1.66,1.50,1.37,1.26],Mp,'spline','extrap');

% Soll-Führungsübertragungsfunktion (System mit dominantem Polpaar)
s = tf('s');
Tsoll = (wntr/tr)^2/(s^2 + 2*xi*wntr/tr*s + (wntr/tr)^2);

% Polüberschuss von T anpassen
linfact = 50/(s+50);
Tsoll = Tsoll*linfact^2;

% Ordnungserhöhung von Tsoll
ws = s+30; % Hurwitz-Polynom!
Tsoll = Tsoll*(ws/ws)^(5+1);

% Reglerberechnung
[R,V] = algsynth(P,Tsoll,0);

% controller design (algebraische Synthese)
% design parameter
contr.s.AlgSynth.label = 'AlgSynth';
contr.s.AlgSynth.desc = 'Algebraische Synthese';
contr.s.AlgSynth.designParam.file.vehicle = stringCvehicle;
contr.s.AlgSynth.designParam.file.steer = stringCsteer;
contr.s.AlgSynth.designParam.P = P;
contr.s.AlgSynth.designParam.vx = vxC;
contr.s.AlgSynth.designParam.lad = ladC;
contr.s.AlgSynth.designParam.Mp = Mp;
contr.s.AlgSynth.designParam.tr = tr;
contr.s.AlgSynth.designParam.Tsoll = Tsoll;

% controller
contr.s.AlgSynth.R.tf = R;
contr.s.AlgSynth.R.ss = ssNormForm(R,'obs');
contr.s.AlgSynth.V.tf = V;
contr.s.AlgSynth.V.ss = ssNormForm(V,'obs');


%% Simulation parameter

% store selected control strategy in parameter-in structure pin
pin.Controller.lka = contr.s.AlgSynth;

% simulation: steering parameter
stringSimSteer = 'paramFile_SteeringMdl_CarMaker_DSR';
pin.SteeringModel.parameterFile = stringSimSteer;
pin.SteeringModel.parameter = paramFile2Struct(stringSimSteer);

% simulation: vehicle parameter
stringSimVehicle = 'paramFile_SingleTrackMdl_BMW5';
pin.VehicleModel.parameterFile = stringSimVehicle;
pin.VehicleModel.parameter = paramFile2Struct(stringSimVehicle);

% simulation: longitudinal velocity vx [m/s]
pin.vx = vxC;

% simulation: look-ahead-distance lad [m]
pin.LAD = ladC;

% load intended trajectory
% [pin.traj,trajErr] = lka_trajectory_08(0,0.05);
pin.traj = ...
	lkaSegmentStraight(0.05,50,0) + ...
	lkaSegmentClothoid(0.05,0,0.005,0,400);
pin.traj_sd = pin.traj.segmentData;

% simulation: interval of integration
tend = pin.traj.segmentData.s(end)/pin.vx - 2*pin.LAD/pin.vx;

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

% run
tic
sim(mdlName);
toc

% lka-controller-string used in simulation
lbl = 'AlgSynth';

% Simulationsprogramm
soli.(lbl).simProg = 'Simulink';

% Zeitstempel
soli.(lbl).simDate = datestr(now);


%% post-processing

% sol.(lbl) = lkaPostProcessing(soli.(lbl),pin,[],simoutState,simoutSensor,simoutSteerAngle);


%% plot

% figure
% lkaPlot(sol,0,'traj');
% lkaPlot(sol,0,'add');
% lkaPlot(sol,pin.lad,'add');
% lkaPlot(sol,0);


%% export

% filename = ['lkaSim_mitLenkung','_traj',num2str(pin.traj.ID.nbr),'_',lbl];
% 
% % save as mat-file
% save(filename,'sol');
% 
% % save as ascii-file
% saveFigure_txt(filename,...
%     1000,...
%     'xtraj',sol.(lbl).simIn.traj.x,...
%     'ytraj',sol.(lbl).simIn.traj.y,...
%     't',sol.(lbl).simOut.t,...
%     'xg',sol.(lbl).simOut.vehicleState.x_g,...
%     'yg',sol.(lbl).simOut.vehicleState.y_g,...
%     'yL_Sensor',sol.(lbl).simOut.lkaSensor.yL,...
%     'epsL_Sensor',sol.(lbl).simOut.lkaSensor.epsL,...
%     'epsL_exactMdl',sol.(lbl).simOut.vehicleState.epsL_exactMdl,...
%     'u',sol.(lbl).simOut.controlInp);
