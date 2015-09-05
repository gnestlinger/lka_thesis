function sys = ssMdl_EPS(sw,paramFile)
% ssMdl_EPS     state-space model of Electric Power Steering System
%   _______
%   Syntax:
%   sys = ssMdl_EPS(sw,paramFile)
%   ________________
%   Input arguments:
%   sw .......... string to choose state space model
%   paramFile ... m-file filename containing parameters
% 
% Source: see corresponding subfunction.
% 
% Subject: lka
% Author: georgnoname
% Date: 16.11.2012 - 26.06.2013


% check input arguments
if ~ischar(sw); error('1st input argument not of type char'); end
% if ~ischar(paramFile); error('2nd input argument not of type char'); end
    

% call subfunction
switch sw
    case {'ParmarHung'}
        sys = parmarHung();
        
    case {'CarMakerDSR'}
        sys = carMakerDSR(paramFile);

    case {'user'}
        sys = user();
        
    otherwise
        error(['Unknown string ''',sw,''' '])
end%switch

end%fcn



function ret = parmarHung()
% Source: "A Sensorless Optimal Control System for an Automotive Electric 
% Power Assist Steering System".
%   ________________
%   State variables:
%   x1 = Winkel Lenksäule
%   x2 = Winkelgeschw. Lenksäule
%   x3 = Winkel Motorsäule
%   x4 = Winkelgeschw. Motorsäule
%   x5 = Verschiebung Zahnstange
%   x6 = Geschw. Verschiebung Zahnstange
%   x7 = Motorstrom
% 
% Author: georgnoname
% Date: 16.11.2012 - 12.12.2012


% EPS-System Parameter
% steering wheel
Jc = 0.04;     % Steering wheel moment of inertia [Kg*m^2]

% steering column
Kc = 172;       % Steering column torsional stiffness [N*m/rad]
Bc = 0.0225;    % Steering column damping [N*m/(rad/s)]
%Jls = 0.9;%Annahme! % Steering column/pinion/motor moment of inertia [Kg*m^2]

% rack
Br = 3920;    % Rack damping [N/(m/s)]
Kt = 23900;   % Tire or rack centering spring rate [N/m]

% motor
Gm = 0.4686;    % Motor gear ratio [1]
Jm = 4.52e-4;   % Motor moment of inertia [Kg*m^2]
Km = 625;       % Motor and gearbox torsional stiffness [N*m(rad/s)]
Bm = 3.339e-3;  % Motor and gearbox damping [N/(m/s)]
km = 0.0345;    % Motor torque and voltage constant [N*m/A]
Lm = 9.06e-5;   % Motor inductance [H]
Rm = 0.035;     % Motor resistance [Ohm]

% other
Mr = 32;       % Rack and weel assembly mass [kg]
rp = 0.0071;  % Pinion radius [m]


% init state space model
sys = ss([],[],[],[]);

% set state space elements
sys.A = [0 1 0 0 0 0 0;...
         -Kc/Jc, -Bc/Jc, 0, 0, Kc/(Jc*rp), 0, 0;...
         0 0 0 1 0 0 0;...
         0, 0, -Km/Jm, -Bm/Jm, Km*Gm/(Jm*rp), 0, km/Jm;...
         0 0 0 0 0 1 0;...
         Kc/(Mr*rp), 0, Km*Gm/(Mr*rp), 0,...
            -(Kt/Mr+Kc/(Mr*rp^2)+Km*Gm^2/(Mr*rp^2)), -Br/Mr, 0;...
         0 0 0 -km/Lm 0 0 -Rm/Lm];
sys.B = [0 0; ...
         1/Jc 0; ...
         0, 0; ...
         0 0
         0 0
         0 0
         0, 1/Lm];
sys.C = [0 0 1 0 0 0 0];
sys.D = 0;

% info
sys.userdata.about = 'State-Space model of EPS-System [Parmar/Hung]';

% sys.userdata.vx.about = 'Längsgeschwindigkeit';
% sys.userdata.vx.value = vx;
% sys.userdata.vx.unit = 'm/s';
% 
% sys.userdata.lad.about = 'look-ahead-distance';
% sys.userdata.lad.value = lad;
% sys.userdata.lad.unit = 'm';

% output argument
ret = sys;

end%fcn


function ret = carMakerDSR(paramFile)
% Source: "CarMaker Reference Manual Version 4.0", S. 140.
%   ________________
%   State variables:
%   (x1 = rack displacement)
%   (x2 = rack displacement velocity)
%   x1 = steering wheel angle
%   x2 = steering wheel velocity
%   ______________
%   Input signals:
%   u = Moment am Lenkrad
% 
% Author: georgnoname
% Date: 11.12.2012 - 15.03.2013


% load parameter
eval(paramFile);

% init state space model
ret = ss([],[],[],[]);

% set state space elements
ret.A = [0 1;...
         0, -1/xi*(iHR^2*drot+drack)];
ret.B = [0; 1/xi*iHR^2*V];
ret.C = [1,0];
ret.D = 0;

% set state/input/output names
ret.StateName = {'SW angle','SW angle Dot'};
ret.Inputname = {'SW torque'};
% ret.OutputName = {};


% info
ret.userdata.about = 'State-Space model of steering-System [CarMaker,13.4]';
ret.userdata.i.about = 'steering angle/rack displacement';
ret.userdata.i.value = iHR;
ret.userdata.i.unit = 'rad/m';
ret.userdata.xi.about = 'shortcut';
ret.userdata.xi.value = xi;
ret.userdata.xi.unit = 'kg';

end%fcn


function ret = user()
% Source: 
%   ________________
%   State variables:
%   (x1 = rack displacement)
%   (x2 = rack displacement velocity)
%   x1 = steering wheel angle
%   x2 = steering wheel velocity
% 
% Author: georgnoname
% Date: 08.01.2013

% steering-System Parameter
J = 0.001;  % Inertia of steering wheel [kg*m^2]
dLS = 0.1;  % rotational damping [Nm*s/rad]
dr = 50;    % translational damping [Nm*s/rad]
mr = 10;    % mass of steering rack [kg]
V = 3;      % Amplification [1]

i = 100;    % steering angle/rack displacement [rad/m]
mL = 10;    % steering mass left [kg]
mR = 10;    % steering mass right [kg]

% motor
J_MC = 0;
iMC = 0.4686;    % Motor gear ratio [1]
J_Mot = 4.52e-4;   % Motor moment of inertia [Kg*m^2]
dMC = 625;       % Motor and gearbox torsional stiffness [N*m(rad/s)]
cMC = 3.339e-3;  % Motor and gearbox damping [N/(m/s)]
k = 0.0345;    % Motor torque and voltage constant [N*m/A]
L = 9.06e-5;   % Motor inductance [H]
R = 0.035;     % Motor resistance [Ohm]


% init state space model
sys = ss([],[],[],[]);

% denominator N
N = 1/(mL+mR+J*i^2+mr+J_MC*iMC^2);


% steering wheel angle
sys.A = [0 1 0 0 0;...
         0, -(dr/i+dLS*i)*N, 0, 0, 0;...
         0 0 0 1 0;...
         cMC*iMC/J_Mot, dMC*iMC/J_Mot, -cMC/J_Mot, -dMC/J_Mot, k/J_Mot;...
         0 0 0 -k/L -R/L];
sys.B = [0 0; V*i*N 0; 0 0; 0 0; 0 1/L];

% sys.C = eye(size(sys.A));
sys.D = 0;

% y = k*i = T_Mot
sys.C = [0 0 0 0 k];
sys.D = 0;

% info
% sys.userdata.about = 'State-Space model of steering-System [CarMaker,13.4]';
% sys.userdata.i.about = 'steering angle/rack displacement';
% sys.userdata.i.value = i;
% sys.userdata.i.unit = '1';

% output argument
ret = sys;

end%fcn
