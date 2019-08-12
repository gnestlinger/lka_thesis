function [sys] = ssMdl_SingleTrack(sw,paramsVhcl,vx,LAD,pFile_Steer)
% SSMDL_SINGLETRACK		returns single-track-based state-space models
%   sys = ssMdl_SingleTrack(SW,PARAMSVHCL,VX,LAD,varargin)
%   ________________
%   Input arguments:
%   SW ............ string to choose state space model (st/stvis/. see below)
%   PARAMSVHCL .... Struct of vehicle parameters
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
% tunable parameter
if nargin < 4 || isempty(LAD)
	LAD = realp('LAD',0);
end%if
if nargin < 3 || isempty(vx)
	vx	= realp('vx',10);
end%if



%%% check input arguments
if ~ischar(sw)
	error('Input argument SW must be of class char!');
end%if
if ~isstruct(paramsVhcl)
	error('Input argument PARAMVHCL must be of class struct!');
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
% STM ...... Single Track Model
% LT ....... Lane Tracking
% yL1int ... 1fach integrierend bezüglich yL
% yL2int ... 2fach integrierend bezüglich yL
% DSR ...... Dynamic Steer Ratio (steering model)


% get state space data
switch lower(sw)
    
    %%% ohne Lenkmodell %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    case {'st'} % Einspurmodell
        [A,B,C,D,InDesc,StateDesc,OutDesc,UserData] = STM(paramsVhcl,vx);
		Name = 'Single Track Model';
        
    case {'stvis'} % Einspurmdl. + lane tracking
        [A,B,C,D,InDesc,StateDesc,OutDesc,UserData] = STM_LT(paramsVhcl,vx,LAD);
		Name = 'Single Track + Lane Tracking Model';
	
	case {'stvis_1int'}
        [A,B,C,D,InDesc,StateDesc,OutDesc,UserData] = STM_LT_yL1int(paramsVhcl,vx,LAD);
		Name = 'Single Track + Lane Tracking + single integral action on lateral offset';
        
    case {'stvis_2int'} % Einspurmdl. + lane tracking + 2fach int. bzgl. yL
        [A,B,C,D,InDesc,StateDesc,OutDesc,UserData] = STM_LT_yL2int(paramsVhcl,vx,LAD);
		Name = 'Single Track + Lane Tracking + double integral action on lateral offset';
		
	
    %%% mit Lenkmodell %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    case {'stdsr'} % Einspurmdl. + Lenkmodell CarMaker DSR
        [A,B,C,D,InDesc,StateDesc,OutDesc,UserData] = STM_DSR(paramsVhcl,vx,pFile_Steer);
		Name = 'Single Track + Steering Model CarMaker DSR';
    
    case {'stvisdsr'} % Einspurmdl. + lane tracking + Lenkmdl. CarMaker DSR
        [A,B,C,D,InDesc,StateDesc,OutDesc,UserData] = STM_LT_DSR(paramsVhcl,vx,LAD,pFile_Steer);
		Name = 'Single Track + Lane Tracking + Steering Model CarMaker DSR';
    
	case {'stvisdsr_1int'}
        [A,B,C,D,InDesc,StateDesc,OutDesc,UserData] = STM_LT_DSR_yL1int(paramsVhcl,vx,LAD,pFile_Steer);
		Name = 'Single Track + Lane Tracking + Steering Model CarMaker DSR + single integral action on lateral offset';
		
    case {'stvisdsr_2int'}
        [A,B,C,D,InDesc,StateDesc,OutDesc,UserData] = STM_LT_DSR_yL2int(paramsVhcl,vx,LAD,pFile_Steer);
		Name = 'Single Track + Lane Tracking + Steering Model CarMaker DSR + double integral action on lateral offset';
    
	
    %%% otherwise %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
    otherwise
        error('Unknown string sw');
        
end%switch


% create state space model
sys = ss(A,B,C,D);

% set input/state/output names
sys.InputName = InDesc(1,:);
sys.InputUnit = InDesc(2,:);
try
	% In Matlab R2012a this does not work!
	sys.StateName = StateDesc(1,:);
	sys.StateUnit = StateDesc(2,:);
catch exc
	warning('This MATLAB release does not support properties "StateName" or "StateUnit" for class GENSS.');
end
sys.OutputName = OutDesc(1,:);
sys.OutputUnit = OutDesc(2,:);

% set info
sys.Name = Name;
sys.UserData = UserData;
sys.UserData.params_Vehicle = paramsVhcl;
sys.UserData.params_Steer	= paramFile2Struct(pFile_Steer);

end%fcn


function [A,B,C,D,InputDesc,StateDesc,OutputDesc,UD] = ...
	STM(params,vx)
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
[csv,csh,lv,lh,m,Iz] = getParams_STM(params);

% set state space elements
A = [-(csh+csv)/(m*vx), (csh*lh-csv*lv)/(m*vx) - vx;...
	(csh*lh-csv*lv)/(Iz*vx), -(csh*lh^2+csv*lv^2)/(Iz*vx)];
B = [csv/m; csv*lv/Iz];
C = eye(size(A));
D = 0;

% set state/input/output names
StateDesc = {...
	'vy','yaw rate';
	'm/s','rad/s'};
InputDesc = {...
	'front wheel angle';
	'rad'};
OutputDesc = StateDesc;

% info
UD.vx.about = 'longitudinal velocity';
UD.vx.value = vx;
UD.vx.unit = 'm/s';

end%fcn

function [A,B,C,D,InputDesc,StateDesc,OutputDesc,UD] = ...
	STM_LT(params,vx,LAD)
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
[csv,csh,lv,lh,m,Iz] = getParams_STM(params);

% set state space elements
A = [-(csh+csv)/(m*vx), (csh*lh-csv*lv)/(m*vx) - vx, 0, 0;...
    (csh*lh-csv*lv)/(Iz*vx), -(csh*lh^2+csv*lv^2)/(Iz*vx), 0, 0;    
    -1, -LAD, 0, vx;...
    0, -1, 0, 0];
B = [[csv/m;csv*lv/Iz;0;0], [0;0;0;vx]];
C = eye(size(A));
D = 0;

% set state/input/output names
StateDesc = {...
	'vy','yawRate','lateralOff','angularDev';
	'm/s','rad/s','m','rad'};
InputDesc = {...
	'steer angle','road curvature at LAD';
	'rad','1/m'};
OutputDesc = StateDesc;

% info
UD.vx.about = 'longitudinal velocity';
UD.vx.value = vx;
UD.vx.unit = 'm/s';
UD.LAD.about = 'look-ahead distance';
UD.LAD.value = LAD;
UD.LAD.unit = 'm';

end%fcn

function [A,B,C,D,InputDesc,StateDesc,OutputDesc,UD] = ...
	STM_LT_yL1int(params,vx,LAD)
% singleTrack_lt_yL1int     state-space model of single track modell + 
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
[csv,csh,lv,lh,m,Iz] = getParams_STM(params);

% set state space elements
A = [-(csh+csv)/(m*vx),(csh*lh-csv*lv)/(m*vx) - vx,0,0,0;...
    (csh*lh-csv*lv)/(Iz*vx),-(csh*lh^2+csv*lv^2)/(Iz*vx),0,0,0;    
    -1,-LAD,0,vx,0;...
    0,-1,0,0,0;...
    0,0,1,0,0];
B = [[csv/m;csv*lv/Iz;0;0;0], [0;0;0;vx;0]];
C = eye(size(A));
D = 0;

% set state/input/output names
StateDesc = {...
	'vy','yawRate','lateralOff','angularDev','Int{lateralOff}';
	'm/s','rad/s','m','rad','m*s'};
InputDesc = {...
	'steer angle','road curvature at LAD';
	'rad','1/m'};
OutputDesc = StateDesc;

% info
UD.vx.about = 'longitudinal velocity';
UD.vx.value = vx;
UD.vx.unit = 'm/s';
UD.LAD.about = 'look-ahead distance';
UD.LAD.value = LAD;
UD.LAD.unit = 'm';

end%fcn

function [A,B,C,D,InputDesc,StateDesc,OutputDesc,UD] = ...
	STM_LT_yL2int(params,vx,LAD)
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
[csv,csh,lv,lh,m,Iz] = getParams_STM(params);

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
StateDesc = {...
	'vy','yawRate','lateralOff','angularDev','IntInt{lateralOff}','Int{lateralOff}';
	'm/s','rad/s','m','rad','m*s^2','m*s'};
InputDesc = {...
	'steer angle','road curvature at LAD';
	'rad','1/m'};
OutputDesc = StateDesc;

% info
UD.vx.about = 'longitudinal velocity';
UD.vx.value = vx;
UD.vx.unit = 'm/s';
UD.LAD.about = 'look-ahead distance';
UD.LAD.value = LAD;
UD.LAD.unit = 'm';

end%fcn



function [A,B,C,D,InputDesc,StateDesc,OutputDesc,UD] = ...
	STM_DSR(paramsVhcl,vx,paramFileSteer)
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

% load parameter
[csv,csh,lv,lh,m,Iz] = getParams_STM(paramsVhcl);
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
StateDesc = {...
	'sy','vy','yawAngle','yawRate','SWAngle','SWAngleDot';
	'm','m/s','rad','rad/s','rad','rad/s'};
InputDesc = {...
	'SWTorque';
	'Nm'};
OutputDesc = StateDesc;

% info
UD.vx.about = 'longitudinal velocity';
UD.vx.value = vx;
UD.vx.unit = 'm/s';

end%fcn

function [A,B,C,D,InputDesc,StateDesc,OutputDesc,UD] = ...
	STM_LT_DSR(paramsVhcl,vx,LAD,paramFileSteer)
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

% load parameter
[csv,csh,lv,lh,m,Iz] = getParams_STM(paramsVhcl);
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
StateDesc = {...
	'vy','yawRate','lateralOff','angularDev','SWAngle','SWAngleDot';
	'm/s','rad/s','m','rad','rad','rad/s'};
InputDesc = {...
	'SWTorque','road curvature at LAD','left tie rod force','right tie rod force';
	'Nm','1/m','N','N'};
OutputDesc = StateDesc;

% info
UD.vx.about = 'longitudinal velocity';
UD.vx.value = vx;
UD.vx.unit = 'm/s';
UD.LAD.about = 'look-ahead distance';
UD.LAD.value = LAD;
UD.LAD.unit = 'm';

end%fcn

function [A,B,C,D,InputDesc,StateDesc,OutputDesc,UD] = ...
	STM_LT_DSR_yL1int(paramsVhcl,vx,LAD,paramFileSteer)
% singleTrack_lt_DSR_yL1int     state-space model of single track modell + 
% lane tracking + steering model + internal model of yL
% 
%   Syntax:
%   ret = singleTrack_lt_DSR_yL1int(paramFile,vx,lad)
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
%   x7 = Int{yL}
% 
% Based on: 'A Comparative Study of Vision-Based Lateral Control Strategies 
% for Autonomous Highway Driving' and CarMaker Manual (see 'Dynamic Steer
% Ratio').
% Source int. Model: 'Dynamic Controller for Lane Keeping and Obstacle 
% Avoidance Assistance System'.
% 

% load parameter
[csv,csh,lv,lh,m,Iz] = getParams_STM(paramsVhcl);
eval(paramFileSteer);

% set state space elements
A = [-(csh+csv)/(m*vx), (csh*lh-csv*lv)/(m*vx) - vx, 0, 0;...
    (csh*lh-csv*lv)/(Iz*vx), -(csh*lh^2+csv*lv^2)/(Iz*vx), 0, 0;    
    -1, -LAD, 0, vx;...
    0 -1 0 0];
B = [csv/m; csv*lv/Iz; 0; 0];

A = [A,B/alph,[0;0;0;0],zeros(4,1);...
    0 0 0 0 0 1 0;...
    0,0,0,0,0,-1/xi*(drot*iHR^2+drack),0;...
    0 0 1 0 0 0 0];
B = [...
	[0;0;0;0;0;iHR^2*V/xi;0],...
	[0;0;0;vx;0;0;0],...
	[0;0;0;0;0;iHR/xi;0],...
	[0;0;0;0;0;-iHR/xi;0],...
	];
C = eye(size(A));
D = 0;

% set state/input/output names
StateDesc = {...
	'vy','yawRate','lateralOff','angularDev','SWAngle','SWAngleDot','Int{x3}';
	'm/s','rad/s','m','rad','rad','rad/s','m*s'};
InputDesc = {...
	'SWTorque','road curvature at LAD','left tie rod force','right tie rod force';
	'Nm','1/m','N','N'};
OutputDesc = StateDesc;

% info
UD.vx.about = 'longitudinal velocity';
UD.vx.value = vx;
UD.vx.unit = 'm/s';
UD.LAD.about = 'look-ahead distance';
UD.LAD.value = LAD;
UD.LAD.unit = 'm';

end%fcn

function [A,B,C,D,InputDesc,StateDesc,OutputDesc,UD] = ...
	STM_LT_DSR_yL2int(paramsVhcl,vx,LAD,paramFileSteer)
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

% load parameter
[csv,csh,lv,lh,m,Iz] = getParams_STM(paramsVhcl);
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
StateDesc = {...
	'vy','yawRate','lateralOff','angularDev','SWAngle','SWAngleDot','IntInt{x3}','Int{x3}';
	'm/s','rad/s','m','rad','rad','rad/s','m*s^2','m*s'};
InputDesc = {...
	'SWTorque','road curvature at LAD','left tie rod force','right tie rod force';
	'Nm','1/m','N','N'};
OutputDesc = StateDesc;

% info
UD.vx.about = 'longitudinal velocity';
UD.vx.value = vx;
UD.vx.unit = 'm/s';
UD.LAD.about = 'look-ahead distance';
UD.LAD.value = LAD;
UD.LAD.unit = 'm';

end%fcn


% --- Get parameters from structure ------------------------------------- %
function [csv,csr,lv,lr,m,Izz] = getParams_STM(S)
csv = S.cs_front;
csr = S.cs_rear;
lv	= S.l_front;
lr	= S.l_rear;
m	= S.m;
Izz	= S.Izz;
end%fcn

function [J,V,alph,drack,drot,iHR,mL,mR,mr,xi] = getParams_CarMakerDSR(S)
J		= S.J;
V		= S.V;
alph	= S.alph;
drack	= S.drack;
drot	= S.drot;
iHR		= S.iHR;
mL		= S.mL;
mR		= S.mR;
mr		= S.mr;
xi		= S.xi;
end%fcn
