function [sys] = ssMdl_SingleTrack(sw,pFile_Vhcl,vx,LAD,pFile_Steer)
% SSMDL_SINGLETRACK		returns single-track-based state-space models
%   sys = ssMdl_SingleTrack(SW,PARAMFILE,VX,LAD,varargin)
%   ________________
%   Input arguments:
%   SW ............ string to choose state space model (st/stvis/. see below)
%   PFILE_VCHL .... m-file filename containing vehicle parameters
%   VX ............ longitudinal velocity [m/s]
%   LAD ........... look-ahead distance [m]
%   PFILE_STEER ... opt. input arguments (e.g. additional parameter files)
%   ________________________________________________
% 	SW can take the following values:
%     - 'st'
%     - 'stvis'
%     - 'stvis_2int'
%     - 'stdsr'
%     - 'stvisdsr'
%     - 'stvisdsr_2int'
% 
% Source: see subfunctions
% 
% Subject: Diplomarbeit - LKA
% $Author$
% $LastChangedDate$
% $Revision$


%%% handle input arguments
if nargin < 5
	pFile_Steer = '';
end%if


%%% check input arguments
if ~ischar(sw)
	error('Input argument SW must be of class char!');
end%if
if ~ischar(pFile_Vhcl)
	error('Input argument PFILE_VHCL must be of class char!');
end%if
if numel(vx) > 1
	error('Input argument VX must be scalar!');
end%if
if numel(LAD) > 1
	error('Input argument LAD must be scalar!');
end%if
if ~ischar(pFile_Steer)
	error('Input argument PFILE_STEER must be of class char!');
end%if
	

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
        [A,B,C,D,IN,IU,SN,SU,ON,OU,UD] = STM(pFile_Vhcl,vx);
		Name = 'Single Track Model';
        
    case {'stvis'} % Einspurmdl. + lane tracking
        [A,B,C,D,IN,IU,SN,SU,ON,OU,UD] = STM_LT(pFile_Vhcl,vx,LAD);
		Name = 'Single Track + Lane Tracking Model';
        
    case {'stvis_2int'} % Einspurmdl. + lane tracking + 2fach int. bzgl. yL
        [A,B,C,D,IN,IU,SN,SU,ON,OU,UD] = STM_LT_yL2int(pFile_Vhcl,vx,LAD);
		Name = 'Single Track + Lane Tracking + double integral action on lateral offset';
		
	
    %%% mit Lenkmodell %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    case {'stdsr'} % Einspurmdl. + Lenkmodell CarMaker DSR
        [A,B,C,D,IN,IU,SN,SU,ON,OU,UD] = STM_DSR(pFile_Vhcl,vx,pFile_Steer);
		Name = 'Single Track + Steering Model CarMaker DSR';
    
    case {'stvisdsr'} % Einspurmdl. + lane tracking + Lenkmdl. CarMaker DSR
        [A,B,C,D,IN,IU,SN,SU,ON,OU,UD] = STM_LT_DSR(pFile_Vhcl,vx,LAD,pFile_Steer);
		Name = 'Single Track + Lane Tracking + Steering Model CarMaker DSR';
        
    case {'stvisdsr_2int'}
        [A,B,C,D,IN,IU,SN,SU,ON,OU,UD] = STM_LT_DSR_yL2int(pFile_Vhcl,vx,LAD,pFile_Steer);
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
sys.UserData.params_Vehicle = paramFile2Struct(pFile_Vhcl);
sys.UserData.params_Steer	= paramFile2Struct(pFile_Steer);

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
