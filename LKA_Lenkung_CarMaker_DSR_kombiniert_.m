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
clear
% close all

% load Simulink model
mdlName = 'LKA_Lenkung_CarMaker_DSR_kombiniert';
load_system(mdlName);


%% LKA-controller design

% LKA-controller design: control plant parameter
paramsCvehicle	= paramFile2Struct('paramFile_SingleTrackMdl_CarMaker_BMW5');
paramsCsteer	= paramFile2Struct('paramFile_SteeringMdl_CarMaker_DSR');

% LKA-controller design: longitudinal velocity vx [m/s]
vxC = 20;

% LKA-controller design: look-ahead distance lad [m]
ladC = 10;

% load state-space model
sysSingleTrackVisDSR = ...
    getMergedSSMdl('stvisdsr',paramsCvehicle,vxC,ladC,paramsCsteer);
sysSingleTrackVisDSR_1int = ss(...
    [[sysSingleTrackVisDSR.a,zeros(6,1)];[0 0 1 0 0 0 0]],...
    [sysSingleTrackVisDSR.b(:,1);0],...
    [sysSingleTrackVisDSR.c(3,:),0],...
    0);
sysSingleTrackVisDSR_2int = ...
    getMergedSSMdl('stvisdsr_2int',paramsCvehicle,vxC,ladC,paramsCsteer);

R = 10;

% 1fach integrierend
Q0int = diag([0 0 1 0 0 0]);
[k0int,~,~] = lqr(sysSingleTrackVisDSR(:,1),Q0int,R);

% 1fach integrierend
Q1int = diag([0 0 1 0 0 0 1]);
[k1int,~,~] = lqr(sysSingleTrackVisDSR_1int,Q1int,R);

% 2fach integrierend
Q2int = diag([0 0 1 0 0 0 1 1]);
% Q = [0 0 0 0 0 0 0 0;...
%      0 0 0 0 0 0 0 0;...
%      0 0 1 0 0 0 0 0;...
%      0 0 0 0 0 0 0 0;...
%      0 0 0 0 0 0 0 0;...
%      0 0 0 0 0 0 0 0;...
%      0 0 0 0 0 0 1 0;...
%      0 0 0 0 0 0 0 1];
[k2int,~,~] = lqr(sysSingleTrackVisDSR_2int(:,1),Q2int,R);

% overall design parameter
contr.t.LQR_DSR.designParam.file.vehicle = paramsCvehicle;
contr.t.LQR_DSR.designParam.file.steer = paramsCsteer;
contr.t.LQR_DSR.designParam.vx = vxC;
contr.t.LQR_DSR.designParam.lad = ladC;

% design parameter and controller values: 0fach int.
contr.t.LQR_DSR.yL0int.designParam.Q = Q2int;
contr.t.LQR_DSR.yL0int.designParam.R = R;
contr.t.LQR_DSR.yL0int.k = [k0int,0,0];

% design parameter and controller values: 1fach int.
contr.t.LQR_DSR.yL1int.designParam.yL1int.Q = Q2int;
contr.t.LQR_DSR.yL1int.designParam.yL1int.R = R;
contr.t.LQR_DSR.yL1int.k = [k1int(1:end-1),0,k1int(end)];

% design parameter and controller values: 2fach int.
contr.t.LQR_DSR.yL2int.designParam.Q = Q2int;
contr.t.LQR_DSR.yL2int.designParam.R = R;
contr.t.LQR_DSR.yL2int.k = k2int;

% wähle Rückführvektor
contr.t.LQR_DSR.k = contr.t.LQR_DSR.yL2int.k;
disp('Verwendeter Rückführvektor:')
contr.t.LQR_DSR.k


%% Simulation parameter

% store selected control strategy in parameter-in structure pin
pin.Controller.lka = contr.t.LQR_DSR;

% simulation: steering parameter
stringSim = 'paramFile_SteeringMdl_CarMaker_DSR';
pin.SteeringModel.parameterFile = stringSim;
pin.SteeringModel.parameter = paramFile2Struct(stringSim);

% simulation: vehicle parameter
stringSim = 'paramFile_SingleTrackMdl_CarMaker_BMW5';
pin.VehicleModel.parameterFile = stringSim;
pin.VehicleModel.parameter = paramFile2Struct(stringSim);

% simulation: longitudinal velocity vx [m/s]
pin.vx = vxC;

% simulation: look-ahead-distance lad [m]
pin.LAD = ladC;

% load intended trajectory
% [pin.traj,trajErr] = lka_trajectory_08(0,0.05);
pin.traj = ...
	LkPathStraight(0.05,50,0) + ...
	LkPathClothoid(0.05,0,0.005,0,400);
pin.traj_sd = pin.traj.pathData;

% simulation: interval of integration
tend = pin.traj.pathData.s(end)/pin.vx - 2*pin.LAD/pin.vx;

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
lbl = 'LQR_2int';

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
