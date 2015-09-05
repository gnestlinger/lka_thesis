% LKA-Simulation: Simulink
% 
% Bewaehrte Solver-Einstellungen
% 
% Subject: lka
% Author: georgnoname
% Date: 30.04.2013


%% set solver parameter

% für Sollbahn mit sprungf. Krümmungsverlauf (kapL wird richtig ermittelt)
set_param(mdlName,'SolverType','Fixed-step')
set_param(mdlName,'Solver','ode3')
set_param(mdlName,'FixedStep','0.01')

% sonst
set_param(mdlName,'SolverType','Variable-step')
set_param(mdlName,'Solver','ode23')
set_param(mdlName,'MaxStep','auto')