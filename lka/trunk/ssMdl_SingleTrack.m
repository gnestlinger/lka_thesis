function [sys] = ssMdl_SingleTrack(sw,paramFile,vx,LAD,varargin)
% SSMDL_SINGLETRACK		returns single-track-based state-space models
%   sys = ssMdl_SingleTrack(SW,PARAMFILE,VX,LAD,varargin)
%   ________________
%   Input arguments:
%   sw .......... string to choose state space model (st/stvis/. see below)
%   paramFile ... m-file filename containing vehicle parameters
%   vx .......... longitudinal velocity [m/s]
%   lad ......... look-ahead distance [m]
%   varargin .... opt. input arguments (e.g. additional parameter files)
%   ________________________________________________
%   State variables: (Index V...Vehicle, g...global, *...additional state) 
%   'st'
%   x1 = Weg in Fahrzeugquerrichtung (sy_V) (*)
%   x2 = Fahrzeugquergeschwindigkeit (vy_V)
%   x3 = Gierwinkel (psi) (*)
%   x4 = Giergeschwindigkeit (psiDot)
% 
%   'stvis'
%   x1 = Fahrzeugquergeschwindigkeit (vy_V)
%   x2 = Giergeschwindigkeit (psiDot)
%   x3 = Querversatz (yL)
%   x4 = Relativwinkel (epsL)
% 
%   'stvis_2int'
%   x1 = Fahrzeugquergeschwindigkeit (vy_V)
%   x2 = Giergeschwindigkeit (psiDot)
%   x3 = Querversatz (yL)
%   x4 = Relativwinkel (epsL)
%   x5 = Int{Int{yL}} (*)
%   x6 = Int{yL} (*)
% 
%   'stdsr'
%   x1 = Weg in Fahrzeugquerrichtung (sy_V) (*)
%   x2 = Fahrzeugquergeschwindigkeit (vy_V)
%   x3 = Gierwinkel (psi) (*)
%   x4 = Giergeschwindigkeit (psiDot)
%   x5 = Lenkradwinkel (deltaH)
%   x6 = Lenkradwinkelgeschwindigkeit (deltaHDot)
% 
%   'stvisdsr'
%   x1 = Fahrzeugquergeschwindigkeit (vy_V)
%   x2 = Giergeschwindigkeit (psiDot)
%   x3 = Querversatz (yL)
%   x4 = Relativwinkel (epsL)
%   x5 = Lenkradwinkel (deltaH)
%   x6 = Lenkradwinkelgeschwindigkeit (deltaHDot)
% 
%   'stvisdsr_2int'
%   x1 = Fahrzeugquergeschwindigkeit (vy_V)
%   x2 = Giergeschwindigkeit (psiDot)
%   x3 = Querversatz (yL)
%   x4 = Relativwinkel (epsL)
%   x5 = Lenkradwinkel (deltaH)
%   x6 = Lenkradwinkelgeschwindigkeit (deltaHDot)
%   x7 = Int{Int{yL}} (*)
%   x8 = Int{yL} (*)
% 
% Source: see subfunctions
% 
% Subject: Diplomarbeit - LKA
% $Author: georgnestlinger $
% $LastChangedDate: 2017-05-24 08:16:45 +0200 (Mi, 24 Mai 2017) $
% $Revision: 132 $


% check input arguments
if ~ischar(sw); error('1st input argument not of type char'); end
if ~ischar(paramFile); error('2nd input argument not of type char'); end
if numel(vx) > 1; error('Dimension of 3rd input argument > 1'); end


% subfunction shortcuts
% st ....... single track
% lt ....... lane tracking
% yL2int ... 2fach integrierend bezüglich yL
% DSR ...... Dynamic Steer Ratio (steering model)

% call subfunction
switch lower(sw)
    
    %%% ohne Lenkmodell %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    case {'st'} % Einspurmodell
        sys = singleTrack(paramFile,vx);
        
    case {'stvis'} % Einspurmdl. + lane tracking
        sys = singleTrack_lt(paramFile,vx,LAD);
        
    case {'stvis_2int'} % Einspurmdl. + lane tracking + 2fach int. bzgl. yL
        sys = singleTrack_lt_yL2int(paramFile,vx,LAD);
        
    %%% mit Lenkmodell %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    case {'stdsr'} % Einspurmdl. + Lenkmodell CarMaker DSR
        sys = singleTrack_DSR(paramFile,vx,varargin{1});
    
    case {'stvisdsr'} % Einspurmdl. + lane tracking + Lenkmdl. CarMaker DSR
        sys = singleTrack_lt_DSR(paramFile,vx,LAD,varargin{1});
        
    case {'stvisdsr_2int'}
        sys = singleTrack_lt_DSR_yL2int(paramFile,vx,LAD,varargin{1});
    
    %%% otherwise %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
    otherwise
        error('Unknown string sw');
        
end%switch


end%fcn



function [sys] = singleTrack(paramFile,vx)
% singleTrack     state-space model of single track model
% 
%   Syntax:
%   ret = singleTrack(paramFile,vx)
% 
%   Input arguments:
%   paramFile ... m-file filename containing vehicle parameters
%   vx .......... longitudinal velocity
% 
% 	State variables: (Index V...Vehicle, g...global)
%   x1 = Weg in Fahrzeugquerrichtung (sy_V) (*)
%   x2 = Fahrzeugquergeschwindigkeit (vy_V)
%   x3 = Gierwinkel (psi) (*)
%   x4 = Giergeschwindigkeit (psiDot)
% 
% Source: A Comparative Study of Vision-Based Lateral Control Strategies 
% for Autonomous Highway Driving. (* ... additional state)
% 


% load parameter
eval(paramFile);

% set state space elements
A = [0 1 0 0;...
	0, -(csh+csv)/(m*vx), 0, (csh*lh-csv*lv)/(m*vx) - vx;...
	0 0 0 1;...
	0, (csh*lh-csv*lv)/(Iz*vx), 0, -(csh*lh^2+csv*lv^2)/(Iz*vx)];
B = [0; csv/m; 0; csv*lv/Iz];
C = eye(size(A));
D = 0;

% create state space model
sys = ss(A,B,C,D);

% set state/input/output names
sys.StateName = {'sy','vy','yaw angle','yaw rate'};
sys.StateUnit = {'m','m/s','rad','rad/s'};
sys.InputName = {'front wheel angle'};
sys.InputUnit = {'rad'};
sys.OutputName = sys.StateName;
sys.OutputUnit = sys.StateUnit;

% info
sys.Name = 'Single Track Model';
sys.UserData.vx.about = 'longitudinal velocity';
sys.UserData.vx.value = vx;
sys.UserData.vx.unit = 'm/s';

end%fcn



function sys = singleTrack_lt(paramFile,vx,LAD)
% singleTrack_lt    state-space model of single track model + 
% lane tracking
%   
%   Syntax:
%   ret = singleTrack_lt(paramFile,vx,lad)
%   
%   Input arguments:
%   paramFile ... m-file filename containing vehicle parameters
%   vx .......... longitudinal velocity
%   lad ......... look-ahead distance
%
%   State variables: (Index V...Vehicle, g...global)
%   x1 = Fahrzeugquergeschwindigkeit (vy_V)
%   x2 = Giergeschwindigkeit (psiDot)
%   x3 = Querversatz (yL)
%   x4 = Relativwinkel (epsL)
% 
% Source: A Comparative Study of Vision-Based Lateral Control Strategies 
% for Autonomous Highway Driving.
% 



% load parameter
eval(paramFile);

% set state space elements
A = [-(csh+csv)/(m*vx), (csh*lh-csv*lv)/(m*vx) - vx, 0, 0;...
    (csh*lh-csv*lv)/(Iz*vx), -(csh*lh^2+csv*lv^2)/(Iz*vx), 0, 0;    
    -1, -LAD, 0, vx;...
    0, -1, 0, 0];
B = [[csv/m;csv*lv/Iz;0;0], [0;0;0;vx]];
C = eye(size(A));
D = 0;

% create state space model
sys = ss(A,B,C,D);

% set state/input/output names
sys.StateName = {'vy','yawRate','lateralOff','angularDev'};
sys.StateUnit = {'m/s','rad/s','m','rad'};
sys.InputName = {'steer angle','road curvature at LAD'};
sys.InputUnit = {'rad','1/m'};
sys.OutputName = sys.StateName;
sys.OutputUnit = sys.StateUnit;


% info
sys.Name = 'Single Track + Lane Tracking Model';
sys.UserData.vx.about = 'longitudinal velocity';
sys.UserData.vx.value = vx;
sys.UserData.vx.unit = 'm/s';
sys.UserData.lad.about = 'look-ahead distance';
sys.UserData.lad.value = LAD;
sys.UserData.lad.unit = 'm';

end%fcn



function sys = singleTrack_lt_yL2int(paramFile,vx,LAD)
% singleTrack_lt_yL2int     state-space model of single track modell + 
% lane tracking + internal model of yL
% 
%   Syntax:
%   ret = singleTrack_lt_yL2int(paramFile,vx,lad)
% 
%   Input arguments:
%   paramFile ... m-file filename containing vehicle parameters
%   vx .......... longitudinal velocity
%   lad ......... look-ahead distance
% 
%   State variables: (Index V...Vehicle, g...global) 
%   x1 = Fahrzeugquergeschwindigkeit (vy_V)
%   x2 = Giergeschwindigkeit (psiDot)
%   x3 = Querversatz (yL)
%   x4 = Relativwinkel (epsL)
%   x5 = Int{Int{yL}} (*)
%   x6 = Int{yL} (*)
%   (* ... additional state)
% 
% Source: A Comparative Study of Vision-Based Lateral Control Strategies 
% for Autonomous Highway Driving.
% Source int. Model: Dynamic Controller for Lane Keeping and Obstacle 
% Avoidance Assistance System.
% 


% load parameter
eval(paramFile);

% set state space elements
A = [-(csh+csv)/(m*vx),(csh*lh-csv*lv)/(m*vx) - vx,0,0,0,0;...
    (csh*lh-csv*lv)/(Iz*vx),-(csh*lh^2+csv*lv^2)/(Iz*vx),0,0,0,0;    
    -1,-LAD,0,vx,0,0;...
    0,-1,0,0,0,0;...
    0,0,0,0,0,1;...
    0,0,1,0,0,0];
B = [[csv/m;csv*lv/Iz;0;0;0;0], [0;0;0;vx;0;0]];
C = eye(size(A));
D = 0;

% create state space model
sys = ss(A,B,C,D);

% set state/input/output names
sys.StateName = {'vy','yawRate','lateralOff','angularDev','IntInt{lateralOff}','Int{lateralOff}'};
sys.StateUnit = {'m/s','rad/s','m','rad','m*s^2','m*s'};
sys.InputName = {'steer angle','road curvature at LAD'};
sys.InputUnit = {'rad','1/m'};
sys.OutputName = sys.StateName;
sys.OutputUnit = sys.StateUnit;


% info
sys.Name = 'Single Track + Lane Tracking + double integral action on lateral offset';
sys.UserData.vx.about = 'longitudinal velocity';
sys.UserData.vx.value = vx;
sys.UserData.vx.unit = 'm/s';
sys.UserData.lad.about = 'look-ahead distance';
sys.UserData.lad.value = LAD;
sys.UserData.lad.unit = 'm';

end%fcn



function sys = singleTrack_DSR(paramFile,vx,paramFileSteer)
% singleTrack_DSR   state-space model of single track modell + steering
% model DSR
% 
%   Syntax:
%   ret = singleTrack_DSR(paramFile,vx,paramFileSteer)
% 
%   Input arguments:
%   paramFile ........ m-file filename containing vehicle parameters
%   vx ............... longitudinal velocity
%   paramFileSteer ... m-file filename containing steering parameters
% 
% 	State variables: (Index V...Vehicle, g...global)
%   x1 = Weg in Fahrzeugquerrichtung (sy_V) (*)
%   x2 = Fahrzeugquergeschwindigkeit (vy_V)
%   x3 = Gierwinkel (psi) (*)
%   x4 = Giergeschwindigkeit (psiDot)
%   x5 = Lenkradwinkel (deltaH)
%   x6 = Lenkradwinkelgeschwindigkeit (deltaHDot)
% 
% Based on: 'A Comparative Study of Vision-Based Lateral Control Strategies 
% for Autonomous Highway Driving' and CarMaker Manual (see 'Dynamic Steer
% Ratio').
% 


% check input arguments
if ~ischar(paramFileSteer); 
    error('Input argument steering-parameter-filename not of type char'); 
end


% load parameter
eval(paramFile);
eval(paramFileSteer);


% init state space model
sys = ss([],[],[],[]);

% set state space elements
A = [0 1 0 0;...
    0, -(csh+csv)/(m*vx), 0, (csh*lh-csv*lv)/(m*vx) - vx;...
    0 0 0 1;...
    0, (csh*lh-csv*lv)/(Iz*vx), 0, -(csh*lh^2+csv*lv^2)/(Iz*vx)];
B = [0; csv/m; 0; csv*lv/Iz];

sys.A = [A,B/alph,[0;0;0;0];...
    0,0,0,0,0,1;...
    0,0,0,0,0,-1/xi*(drot*iHR^2+drack)];
sys.B = [0; 0; 0; 0; 0; iHR^2*V/xi];
sys.C = eye(size(sys.A));
sys.D = 0;

% set state/input/output names
sys.StateName = {'sy','vy','yawAngle','yawRate','SWAngle','SWAngleDot'};
sys.StateUnit = {'m','m/s','rad','rad/s','rad','rad/s'};
sys.InputName = {'SWTorque'};
sys.InputUnit = {'Nm'};
sys.OutputName = sys.StateName;
sys.OutputUnit = sys.StateUnit;


% info
sys.Name = 'Einspurmodell + Lenkmodell CarMaker DSR';
sys.UserData.vx.about = 'longitudinal velocity';
sys.UserData.vx.value = vx;
sys.UserData.vx.unit = 'm/s';

end%fcn



function sys = singleTrack_lt_DSR(paramFile,vx,lad,paramFileSteer)
% singleTrack_lt_DSR    state-space model of single track modell + lane
% tracking + steering model DSR
% 
%   Syntax:
%   ret = singleTrack_lt_DSR(paramFile,vx,lad,paramFileSteer)
% 
%   Input arguments:
%   paramFile ........ m-file filename containing vehicle parameters
%   vx ............... longitudinal velocity
%   lad .............. look-ahead distance
%   paramFileSteer ... m-file filename containing steering parameters
% 
% 	State variables: (Index V...Vehicle, g...global)
%   x1 = Fahrzeugquergeschwindigkeit (vy_V)
%   x2 = Giergeschwindigkeit (psiDot)
%   x3 = Querversatz (yL)
%   x4 = Relativwinkel (epsL)
%   x5 = Lenkradwinkel (deltaH)
%   x6 = Lenkradwinkelgeschwindigkeit (deltaHDot)
% 
% Based on: 'A Comparative Study of Vision-Based Lateral Control Strategies 
% for Autonomous Highway Driving' and CarMaker Manual (see 'Dynamic Steer
% Ratio').
% 


% check input arguments
if ~ischar(paramFileSteer); 
    error('Input argument steering-parameter-filename not of type char'); 
end


% load parameter
eval(paramFile);
eval(paramFileSteer);

% init state space model
sys = ss([],[],[],[]);

% set state space elements
A = [-(csh+csv)/(m*vx),(csh*lh-csv*lv)/(m*vx) - vx,0,0;...
    (csh*lh-csv*lv)/(Iz*vx),-(csh*lh^2+csv*lv^2)/(Iz*vx),0,0;    
    -1,-lad,0,vx;...
    0,-1,0,0];
B = [csv/m; csv*lv/Iz; 0; 0];

sys.A = [A,B/alph,[0;0;0;0];...
    0,0,0,0,0,1;...
    0,0,0,0,0,-1/xi*(drot*iHR^2+drack)];
sys.B = [...
	[0;0;0;0;0;iHR^2*V/xi],...
	[0;0;0;vx;0;0],...
	[0;0;0;0;0;iHR/xi],...
	[0;0;0;0;0;-iHR/xi],...
	];
sys.C = eye(size(sys.A));
sys.D = 0;

% set state/input/output names
sys.StateName = {'vy','yawRate','lateralOff','angularDev','SWAngle','SWAngleDot'};
sys.StateUnit = {'m/s','rad/s','m','rad','rad','rad/s'};
sys.InputName = {...
	'SWTorque',...
	'road curvature at LAD',...
	'left tie rod force',...
	'right tie rod force',...
	};
sys.InputUnit = {'Nm','1/m','N','N'};
sys.OutputName = sys.StateName;
sys.OutputUnit = sys.StateUnit;


% info
sys.Name = 'Single Track + Lane Tracking + Steering Model CarMaker DSR';
sys.UserData.vx.about = 'longitudinal velocity';
sys.UserData.vx.value = vx;
sys.UserData.vx.unit = 'm/s';
sys.UserData.lad.about = 'look-ahead distance';
sys.UserData.lad.value = lad;
sys.UserData.lad.unit = 'm';

end%fcn



function sys = singleTrack_lt_DSR_yL2int(paramFile,vx,lad,paramFileSteer)
% singleTrack_lt_DSR_yL2int     state-space model of single track modell + 
% lane tracking + steering model + internal model of yL
% 
%   Syntax:
%   ret = singleTrack_lt_DSR_yL2int(paramFile,vx,lad)
% 
%   Input arguments:
%   paramFile ........ m-file filename containing vehicle parameters
%   vx ............... longitudinal velocity
%   lad .............. look-ahead distance
%   paramFileSteer ... m-file filename containing steering parameters
% 
% 	State variables: (Index V...Vehicle, g...global)
%   x1 = Fahrzeugquergeschwindigkeit (vy_V)
%   x2 = Giergeschwindigkeit (psiDot)
%   x3 = Querversatz (yL)
%   x4 = Relativwinkel (epsL)
%   x5 = Lenkradwinkel (deltaH)
%   x6 = Lenkradwinkelgeschwindigkeit (deltaHDot)
%   x7 = Int{x8} = Int{Int{yL}}
%   x8 = Int{yL}
% 
% Based on: 'A Comparative Study of Vision-Based Lateral Control Strategies 
% for Autonomous Highway Driving' and CarMaker Manual (see 'Dynamic Steer
% Ratio').
% Source int. Model: 'Dynamic Controller for Lane Keeping and Obstacle 
% Avoidance Assistance System'.
% 


% check input arguments
if ~ischar(paramFileSteer); 
    error('Input argument steering-parameter-filename not of type char'); 
end


% load parameter
eval(paramFile);
eval(paramFileSteer);

% init state space model
sys = ss([],[],[],[]);

% set state space elements
A = [-(csh+csv)/(m*vx), (csh*lh-csv*lv)/(m*vx) - vx, 0, 0;...
    (csh*lh-csv*lv)/(Iz*vx), -(csh*lh^2+csv*lv^2)/(Iz*vx), 0, 0;    
    -1, -lad, 0, vx;...
    0 -1 0 0];
B = [csv/m; csv*lv/Iz; 0; 0];

sys.A = [A,B/alph,[0;0;0;0],zeros(4,2);...
    0 0 0 0 0 1 0 0;...
    0,0,0,0,0,-1/xi*(drot*iHR^2+drack),0,0;...
    0 0 0 0 0 0 0 1;...
    0 0 1 0 0 0 0 0];
sys.B = [...
	[0;0;0;0;0;iHR^2*V/xi;0;0],...
	[0;0;0;vx;0;0;0;0],...
	[0;0;0;0;0;iHR/xi;0;0],...
	[0;0;0;0;0;-iHR/xi;0;0],...
	];
sys.C = eye(size(sys.A));
sys.D = 0;

% set state/input/output names
sys.StateName = {'vy','yawRate','lateralOff','angularDev','SWAngle','SWAngleDot',...
    'IntInt{x3}','Int{x3}'};
sys.StateUnit = {'m/s','rad/s','m','rad','rad','rad/s','m*s^2','m*s'};
sys.InputName = {...
	'SWTorque',...
	'road curvature at LAD',...
	'left tie rod force',...
	'right tie rod force',...
	};
sys.InputUnit = {'Nm','1/m','N','N'};
sys.OutputName = sys.StateName;
sys.OutputUnit = sys.StateUnit;


% info
sys.Name = 'Single Track + Lane Tracking + Steering Model CarMaker DSR + 2fach int. bzgl. yL';
sys.UserData.vx.about = 'longitudinal velocity';
sys.UserData.vx.value = vx;
sys.UserData.vx.unit = 'm/s';
sys.UserData.lad.about = 'look-ahead distance';
sys.UserData.lad.value = lad;
sys.UserData.lad.unit = 'm';

end%fcn
