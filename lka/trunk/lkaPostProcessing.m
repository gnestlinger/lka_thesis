function solOut = lkaPostProcessing(soli,pIn,uIn,state,sensor,contrin)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here


switch soli.simProg
    case {'Matlab'}
        % keep simulation program
        solOut.simProg = soli.simProg;

        % keep timestamp
        solOut.simDate = soli.simDate;

        % create simIn field (simulation parameter)
        solOut.simIn = simIn_Matlab(pIn,uIn);

        % create simOut field (simulation output)
        solOut.simOut = simOut_Matlab( soli );

        % create simPost field (calculations with simulation output)
        solOut.simPost = simPost_Matlab(solOut,pIn,soli.y);

    case {'Simulink'}
        % keep simulation program
        solOut.simProg = soli.simProg;

        % keep timestamp
        solOut.simDate = soli.simDate;

        % create simIn field (simulation parameter)
        solOut.simIn = simIn_Simulink(pIn);

        % create simOut field (simulation output)
        solOut.simOut = simOut_Simulink(state,sensor,contrin);

        % create simPost field (calculations with simulation output)
        solOut.simPost = simPost_Simulink(solOut,pIn);

    otherwise
        error('')

end%switch


end%fcn



function simIn = simIn_Matlab(pin,uin)

%%% sim in %%%
% sim in: Simulationsparameter
simIn = pin;

% sim in: Regler
simIn.Controller.lka.uIn = uin;

end%fcn


function simIn = simIn_Simulink(pin)

%%% sim in %%%
% sim in: Simulationsparameter
simIn = pin;

% sim in: Regler
% simIn.Controller.lka.u = pin.Controller.lka;

end%fcn



function simOut = simOut_Matlab(solin)

%%% sim out %%%
% sim out: vehicle states
simOut.t = solin.x';

state.signals.values = solin.y';
simOut.vehicleState.sy = state.signals.values(:,1)';
simOut.vehicleState.vy = state.signals.values(:,2)';
simOut.vehicleState.psi = state.signals.values(:,3)';
simOut.vehicleState.psiDot = state.signals.values(:,4)';
simOut.vehicleState.x_g = state.signals.values(:,5)';
simOut.vehicleState.y_g = state.signals.values(:,6)';
simOut.vehicleState.s_g = state.signals.values(:,7)';
try
    simOut.vehicleState.yL = state.signals.values(:,8)';
    simOut.vehicleState.epsL = state.signals.values(:,9)';
catch exception
    simOut.vehicleState.yL = [];
    simOut.vehicleState.epsL = [];    
end%try

% LKA-Sensor data % nur bei Simulink-Simulation vorhanden
% solOut.(lbli).simOut.lkaSensor.yL = sensor.signals.values(:,1)';
% solOut.(lbli).simOut.lkaSensor.epsL = sensor.signals.values(:,2)';
% solOut.(lbli).simOut.lkaSensor.kapL = sensor.signals.values(:,3)';

% control input % nur bei Simulink-Simulation vorhanden
% solOut.(lbli).simOut.controlInp.steeringAngle = contrin.signals.values;

end%fcn


function simOut = simOut_Simulink(state,sensor,contrin)

%%% sim out %%%
    % sim out: vehicle states
    simOut.t = state.time';
    simOut.vehicleState.sy = state.signals.values(:,1);
    simOut.vehicleState.vy = state.signals.values(:,2);
    simOut.vehicleState.psi = state.signals.values(:,3);
    simOut.vehicleState.psiDot = state.signals.values(:,4);
    simOut.vehicleState.x_g = state.signals.values(:,5);
    simOut.vehicleState.y_g = state.signals.values(:,6);
    simOut.vehicleState.s_g = state.signals.values(:,7);
    try
        simOut.vehicleState.yL_linMdl = state.signals.values(:,8)';
        simOut.vehicleState.epsL_linMdl = state.signals.values(:,9)';
        simOut.vehicleState.yL_exactMdl = state.signals.values(:,10)';
        simOut.vehicleState.epsL_exactMdl = state.signals.values(:,11)';
    catch exception
        simOut.vehicleState.yL = [];
        simOut.vehicleState.epsL = [];    
    end%try

    % LKA-Sensor data
    simOut.lkaSensor.yL = sensor.signals.values(:,1)';
    simOut.lkaSensor.epsL = sensor.signals.values(:,2)';
    simOut.lkaSensor.kapL = sensor.signals.values(:,3)';

    % control input
    simOut.controlInp = contrin.signals.values;

end%fcn



function simPost = simPost_Matlab(solin,pin,state)

ladC = solin.simIn.Controller.lka.designParam.lad;
simPost.relPos = calc_relPos(solin,[0,ladC,pin.lad]);
% sol.(lbl{i}).simPost.relPos = ret;

fh = solin.simIn.Controller.lka.uIn;
ret = calc_ladIndependentData(solin.simOut.t,state,fh);
simPost.controlInp.steeringAngle = ret.u;

end%fcn


function simPost = simPost_Simulink(solIn,pin)

ladC = solIn.simIn.Controller.lka.designParam.lad;
% sim post: Querversatz yL, Relativwinkel epsL
simPost.relPos = calc_relPos(solIn,[0,ladC,pin.lad]);

end%fcn
