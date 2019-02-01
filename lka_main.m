% LKA-Simulation: Matlab
% 
% Subject: lka
% Author: georgnoname
% Date: 24.09.2012 - 24.02.2013

clc
clear all
% close all


%% LKA-controller design

% controller design: vehicle parameter
vehicleParameterC = 'paramFile_SingleTrackMdl_BMW5';

% controller design: longitudinal velocity vx [m/s]
vxC = 20;

% controller design: look-ahead distance lad [m]
ladC = 5;

% controler design: calculate controler with upper parameters
contr.t = lkaController_t(vehicleParameterC,vxC,ladC);

% get list of control strategies
fNamesContr = fieldnames(contr.t);
disp('Available controllers:');
disp(fNamesContr);


%% Simulation parameter

% simulation: vehicle parameter
string = 'paramFile_SingleTrackMdl_BMW5';
pin.VehicleModel.singleTrack.parameterFile = string;
pin.VehicleModel.singleTrack.parameter = loadParameter(string,'about');

% load intended trajectory
[pin.traj,err] = lka_trajectory_07(0,0.05);

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


%% run simulation

%       1    2    3       4          5
lbl = {'P','ZR','LQR','ZR_2int','LQR_2int'};

for i = 1:5

    disp(['Controller: ',lbl{i}]);

    % store selected control strategy in parameter-in structure pin
    pin.Controller.lka = contr.t.(lbl{i});

    % simulation: longitudinal velocity vx [m/s]
    pin.vx = contr.t.(lbl{i}).designParam.vx;
    % pin.vx = 30;

    % simulation: look-ahead-distance lad [m]
    pin.lad = contr.t.(lbl{i}).designParam.lad;
    % pin.lad = 5;

    % interval of integration (make sure: ladC des Reglers darf bei t->tend
    % nicht ¸ber das Ende der Trajektorie hinausschauen => -x*pin.lad/pin.vx
    % entsprechend groﬂ w‰hlen) 
    tend = pin.traj.s(end)/pin.vx - 3*pin.lad/pin.vx;
    tspan = [0,tend];

    % choose controller (function handle)
    % u_fh = contr.t.(lbl{i}).u_fh;
    uin = @(t,x) einspurMdl_u(t,x,pin,contr.t.(lbl{i}).u_fh);


    % solver
    % opts = odeset('RelTol',1e-6);
    opts = odeset('MaxStep',0.1);


    tic
    try
        xInit = [pin.VehicleModel.initialValue.sy;...
            pin.VehicleModel.initialValue.vy;...
            pin.VehicleModel.initialValue.psi;...
            pin.VehicleModel.initialValue.psiDot;...
            pin.VehicleModel.initialValue.x_g;...
            pin.VehicleModel.initialValue.y_g;...
            pin.VehicleModel.initialValue.s_g];
        soli.(lbl{i}) = ode45(@(t,x) einspurMdl_dgl_vy(t,x,uin,pin),...
            tspan,xInit,opts);

    catch expression
    %     disp(expression.message);
        xInit = [pin.VehicleModel.initialValue.sy;...
            pin.VehicleModel.initialValue.vy;...
            pin.VehicleModel.initialValue.psi;...
            pin.VehicleModel.initialValue.psiDot;...
            pin.VehicleModel.initialValue.x_g;...
            pin.VehicleModel.initialValue.y_g;...
            pin.VehicleModel.initialValue.s_g;0;0];
        soli.(lbl{i}) = ode45(@(t,x) einspurMdl_dgl_vy_int(t,x,uin,pin),...
            tspan,xInit,opts);

    end%try
    toc

    %%% post-processing

    % Simulationsprogramm
    soli.(lbl{i}).simProg = 'Matlab';

    % Zeitstempel
    soli.(lbl{i}).simDate = datestr(now);

    % post-processing
    sol.(lbl{i}) = lkaPostProcessing(soli.(lbl{i}),pin,uin);
    
    
end%for



%% plot

% figure
% lkaPlot(sol,0,'traj');
% lkaPlot(sol,0,'add');
lkaPlot(sol,ladC,'add');
% lkaPlot(sol,0);

%% export

% contrStrat = fNamesContr{i};
% 
% % filename =
% % ['sol_',datestr(datenum(sol.(contrStrat).simDate),'yymmdd_HHMM')];
% filename = ['sol_','traj1'];
% filename = [filename,'_',contrStrat];
% filename = [filename,'_dummy'];
% 
% % save as mat-file
% save(filename,'sol');
% 
% % save as ascii-file
% saveFigure_txt(filename,...
%     length(sol.(contrStrat).x),...
%     'xtraj',sol.(contrStrat).simParameter.traj.x,...
%     'ytraj',sol.(contrStrat).simParameter.traj.y,...
%     't',sol.(contrStrat).x,...
%     'x0',sol.(fNamesContr{i}).y(5,:),...
%     'y0',sol.(fNamesContr{i}).y(6,:),...
%     'yL_0',sol.(fNamesContr{i}).addData.lad_0.yL.value,...
%     'yL_ladC',sol.(fNamesContr{i}).addData.(['lad_',num2str(ladC)]).yL.value,...
%     'u',rad2deg(sol.(fNamesContr{i}).addData.u.value),...
%     's',sol.(fNamesContr{i}).addData.s.value);


% saveFigure_txt(filename,...
%     'xtraj',sol.(contrStrat).simParameter.traj.x,...
%     'ytraj',sol.(contrStrat).simParameter.traj.y,...
%     't',sol.(contrStrat).x,...
%     'x0',sol.(fNamesContr{i}).y(5,:),...
%     'y0',sol.(fNamesContr{i}).y(6,:),...
%     'yL_0',sol.(fNamesContr{i}).addData.lad_0.yL.value,...
%     'yL_ladC',sol.(fNamesContr{i}).addData.(['lad_',num2str(ladC)]).yL.value,...
%     'u',rad2deg(sol.(fNamesContr{i}).addData.u.value),...
%     's',sol.(fNamesContr{i}).addData.s.value);

