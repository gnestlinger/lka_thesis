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
% $Author$
% $LastChangedDate$
% $Revision$


% check input arguments
if ~ischar(sw); error('1st input argument not of type char'); end
if ~ischar(paramFile); error('2nd input argument not of type char'); end
if numel(vx) > 1; error('Dimension of 3rd input argument > 1'); end


% subfunction shortcuts
% st ....... single track
% lt ....... lane tracking
% yL2int ... 2fach integrierend bezüglich yL
% DSR ...... Dynamic Steer Ratio (steering model)

% tunable parameter
if isempty(vx);		vx	= realp('vx',10); end%if
if isempty(LAD);	LAD = realp('LAD',0); end%if

% get state space data
switch lower(sw)
    
    %%% ohne Lenkmodell %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    case {'st'} % Einspurmodell
        [A,B,C,D,IN,IU,SN,SU,ON,OU,UD] = STM(paramFile,vx);
		Name = 'Single Track Model';
        
    case {'stvis'} % Einspurmdl. + lane tracking
        [A,B,C,D,IN,IU,SN,SU,ON,OU,UD] = STM_LT(paramFile,vx,LAD);
		Name = 'Single Track + Lane Tracking Model';
        
    case {'stvis_2int'} % Einspurmdl. + lane tracking + 2fach int. bzgl. yL
        [A,B,C,D,IN,IU,SN,SU,ON,OU,UD] = STM_LT_yL2int(paramFile,vx,LAD);
		Name = 'Single Track + Lane Tracking + double integral action on lateral offset';
		
	
    %%% mit Lenkmodell %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    case {'stdsr'} % Einspurmdl. + Lenkmodell CarMaker DSR
        [A,B,C,D,IN,IU,SN,SU,ON,OU,UD] = STM_DSR(paramFile,vx,varargin{1});
		Name = 'Single Track + Steering Model CarMaker DSR';
    
    case {'stvisdsr'} % Einspurmdl. + lane tracking + Lenkmdl. CarMaker DSR
        [A,B,C,D,IN,IU,SN,SU,ON,OU,UD] = STM_LT_DSR(paramFile,vx,LAD,varargin{1});
		Name = 'Single Track + Lane Tracking + Steering Model CarMaker DSR';
        
    case {'stvisdsr_2int'}
        [A,B,C,D,IN,IU,SN,SU,ON,OU,UD] = STM_LT_DSR_yL2int(paramFile,vx,LAD,varargin{1});
		Name = 'Single Track + Lane Tracking + Steering Model CarMaker DSR + double integral action on lateral offset';
    
		
    %%% otherwise %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
    otherwise
        error('Unknown string sw');
        
end%switch


% create state space model
sys = ss(A,B,C,D);

% set input/state/output names
sys.InputName = IN;
sys.InputUnit = IU;
try
	% In Matlab R2012a this does not work!
	sys.StateName = SN;
	sys.StateUnit = SU;
catch exc
	warning('This MATLAB release does not support properties "StateName" or "StateUnit" for class GENSS.');
end
sys.OutputName = ON;
sys.OutputUnit = OU;

% set info
sys.Name = Name;
sys.UserData = UD;

end%fcn



function [A,B,C,D,InputName,InputUnit,StateName,StateUnit,OutputName,OutputUnit,UD] = ...
	STM(paramFile,vx)
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

% set state/input/output names
StateName = {'sy','vy','yaw angle','yaw rate'};
StateUnit = {'m','m/s','rad','rad/s'};
InputName = {'front wheel angle'};
InputUnit = {'rad'};
OutputName = StateName;
OutputUnit = StateUnit;

% info
UD.vx.about = 'longitudinal velocity';
UD.vx.value = vx;
UD.vx.unit = 'm/s';

end%fcn



function [A,B,C,D,InputName,InputUnit,StateName,StateUnit,OutputName,OutputUnit,UD] = ...
	STM_LT(paramFile,vx,LAD)
% singleTrack_lt    state-space model of single track model + lane tracking
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

% set state/input/output names
StateName = {'vy','yawRate','lateralOff','angularDev'};
StateUnit = {'m/s','rad/s','m','rad'};
InputName = {'steer angle','road curvature at LAD'};
InputUnit = {'rad','1/m'};
OutputName = StateName;
OutputUnit = StateUnit;

% info
UD.vx.about = 'longitudinal velocity';
UD.vx.value = vx;
UD.vx.unit = 'm/s';
UD.LAD.about = 'look-ahead distance';
UD.LAD.value = LAD;
UD.LAD.unit = 'm';

end%fcn



function [A,B,C,D,InputName,InputUnit,StateName,StateUnit,OutputName,OutputUnit,UD] = ...
	STM_LT_yL2int(paramFile,vx,LAD)
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

% set state/input/output names
StateName = {'vy','yawRate','lateralOff','angularDev','IntInt{lateralOff}','Int{lateralOff}'};
StateUnit = {'m/s','rad/s','m','rad','m*s^2','m*s'};
InputName = {'steer angle','road curvature at LAD'};
InputUnit = {'rad','1/m'};
OutputName = StateName;
OutputUnit = StateUnit;

% info
UD.vx.about = 'longitudinal velocity';
UD.vx.value = vx;
UD.vx.unit = 'm/s';
UD.LAD.about = 'look-ahead distance';
UD.LAD.value = LAD;
UD.LAD.unit = 'm';

end%fcn



function [A,B,C,D,InputName,InputUnit,StateName,StateUnit,OutputName,OutputUnit,UD] = ...
	STM_DSR(paramFile,vx,paramFileSteer)
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

% set state space elements
A = [0 1 0 0;...
    0, -(csh+csv)/(m*vx), 0, (csh*lh-csv*lv)/(m*vx) - vx;...
    0 0 0 1;...
    0, (csh*lh-csv*lv)/(Iz*vx), 0, -(csh*lh^2+csv*lv^2)/(Iz*vx)];
B = [0; csv/m; 0; csv*lv/Iz];

A = [A,B/alph,[0;0;0;0];...
    0,0,0,0,0,1;...
    0,0,0,0,0,-1/xi*(drot*iHR^2+drack)];
B = [0; 0; 0; 0; 0; iHR^2*V/xi];
C = eye(size(A));
D = 0;

% set state/input/output names
StateName = {'sy','vy','yawAngle','yawRate','SWAngle','SWAngleDot'};
StateUnit = {'m','m/s','rad','rad/s','rad','rad/s'};
InputName = {'SWTorque'};
InputUnit = {'Nm'};
OutputName = StateName;
OutputUnit = StateUnit;

% info
UD.vx.about = 'longitudinal velocity';
UD.vx.value = vx;
UD.vx.unit = 'm/s';

end%fcn



function [A,B,C,D,InputName,InputUnit,StateName,StateUnit,OutputName,OutputUnit,UD] = ...
	STM_LT_DSR(paramFile,vx,LAD,paramFileSteer)
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

% set state space elements
A = [-(csh+csv)/(m*vx),(csh*lh-csv*lv)/(m*vx) - vx,0,0;...
    (csh*lh-csv*lv)/(Iz*vx),-(csh*lh^2+csv*lv^2)/(Iz*vx),0,0;    
    -1,-LAD,0,vx;...
    0,-1,0,0];
B = [csv/m; csv*lv/Iz; 0; 0];

A = [A,B/alph,[0;0;0;0];...
    0,0,0,0,0,1;...
    0,0,0,0,0,-1/xi*(drot*iHR^2+drack)];
B = [...
	[0;0;0;0;0;iHR^2*V/xi],...
	[0;0;0;vx;0;0],...
	[0;0;0;0;0;iHR/xi],...
	[0;0;0;0;0;-iHR/xi],...
	];
C = eye(size(A));
D = 0;

% set state/input/output names
StateName = {'vy','yawRate','lateralOff','angularDev','SWAngle','SWAngleDot'};
StateUnit = {'m/s','rad/s','m','rad','rad','rad/s'};
InputName = {...
	'SWTorque',...
	'road curvature at LAD',...
	'left tie rod force',...
	'right tie rod force',...
	};
InputUnit = {'Nm','1/m','N','N'};
OutputName = StateName;
OutputUnit = StateUnit;

% info
UD.vx.about = 'longitudinal velocity';
UD.vx.value = vx;
UD.vx.unit = 'm/s';
UD.LAD.about = 'look-ahead distance';
UD.LAD.value = LAD;
UD.LAD.unit = 'm';

end%fcn



function [A,B,C,D,InputName,InputUnit,StateName,StateUnit,OutputName,OutputUnit,UD] = ...
	STM_LT_DSR_yL2int(paramFile,vx,LAD,paramFileSteer)
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

% set state space elements
A = [-(csh+csv)/(m*vx), (csh*lh-csv*lv)/(m*vx) - vx, 0, 0;...
    (csh*lh-csv*lv)/(Iz*vx), -(csh*lh^2+csv*lv^2)/(Iz*vx), 0, 0;    
    -1, -LAD, 0, vx;...
    0 -1 0 0];
B = [csv/m; csv*lv/Iz; 0; 0];

A = [A,B/alph,[0;0;0;0],zeros(4,2);...
    0 0 0 0 0 1 0 0;...
    0,0,0,0,0,-1/xi*(drot*iHR^2+drack),0,0;...
    0 0 0 0 0 0 0 1;...
    0 0 1 0 0 0 0 0];
B = [...
	[0;0;0;0;0;iHR^2*V/xi;0;0],...
	[0;0;0;vx;0;0;0;0],...
	[0;0;0;0;0;iHR/xi;0;0],...
	[0;0;0;0;0;-iHR/xi;0;0],...
	];
C = eye(size(A));
D = 0;

% set state/input/output names
StateName = {'vy','yawRate','lateralOff','angularDev','SWAngle','SWAngleDot',...
    'IntInt{x3}','Int{x3}'};
StateUnit = {'m/s','rad/s','m','rad','rad','rad/s','m*s^2','m*s'};
InputName = {...
	'SWTorque',...
	'road curvature at LAD',...
	'left tie rod force',...
	'right tie rod force',...
	};
InputUnit = {'Nm','1/m','N','N'};
OutputName = StateName;
OutputUnit = StateUnit;

% info
UD.vx.about = 'longitudinal velocity';
UD.vx.value = vx;
UD.vx.unit = 'm/s';
UD.LAD.about = 'look-ahead distance';
UD.LAD.value = LAD;
UD.LAD.unit = 'm';

end%fcn
