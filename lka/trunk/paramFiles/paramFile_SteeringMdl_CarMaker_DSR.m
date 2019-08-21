% PARAMETER FILE
% 
% Steering model 
% CarMaker "Dynamic Steer Ratio"
% 
% Source: IPG Carmaker 4.0
% 

% Subject: Parameter file
% $Author$
% $LastChangedDate$
% $Revision$


about.steering = 'CarMaker "Dynamic Steer Ratio"';
about.source = 'IPG Carmaker 4.0';

% steering system parameter
d_rot	= 0.1;	% rotational damping [Nm*s/rad]
d_rack	= 50;	% damping of steering rack [N*s/m]
i_HR	= 100;  % steering angle/rack displacement [rad/m]
I_rot	= 0.001;  % Inertia of steering wheel [kg*m^2]
m_rack	= 10;	% mass of steering rack [kg]
m_left	= 10;   % steering mass left [kg]
m_right = 10;   % steering mass right [kg]
V		= 3;	% Amplification [-]

steeringRatio = 13;	% Übersetzung Lenkrad-/Radlenkwinkel [-]
