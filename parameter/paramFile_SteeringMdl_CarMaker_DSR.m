% PARAMETER FILE
% 
% Steering model 
% CarMaker "Dynamic Steer Ratio"
% 
% Source: IPG Carmaker 4.0
% 

% Subject: 
% $Author: georgnoname@gmail.com $
% $LastChangedDate: 2015-08-12 11:19:04 +0200 (Mi, 12 Aug 2015) $
% $Revision: 196 $


about.steering = 'CarMaker "Dynamic Steer Ratio"';
about.source = 'IPG Carmaker 4.0';

% steering system parameter
drot = 0.1; % rotational damping [Nm*s/rad]
drack = 50; % damping of steering rack [Nm*s/rad]
iHR = 100;  % steering angle/rack displacement [rad/m]
J = 0.001;  % Inertia of steering wheel [kg*m^2]
mr = 10;    % mass of steering rack [kg]
mL = 10;    % steering mass left [kg]
mR = 10;    % steering mass right [kg]
V = 3;      % Amplification [1]

alph = 13;	% Übersetzung Lenkrad-/Radlenkwinkel [1]

xi = mL + mR + J*iHR^2 + mr;    % shortcut [kg]

