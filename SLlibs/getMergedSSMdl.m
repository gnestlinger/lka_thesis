function [sys] = getMergedSSMdl(sw,paramsVhcl,vx,LAD,paramsSteer)
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
	paramsSteer = struct();
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
if ~isstruct(paramsSteer)
	error('Input argument PFILE_STEER must be of class struct!');
end%if


% subfunction shortcuts
% STM ...... Single Track Model
% LT ....... Lane Tracking
% yL1int ... 1fach integrierend bezüglich yL
% yL2int ... 2fach integrierend bezüglich yL
% DSR ...... Dynamic Steer Ratio (steering model)


% get state space data
switch lower(sw)
        
    case {'stvis'} % Einspurmdl. + lane tracking
        [sys,UserData] = STM_LT(paramsVhcl,vx,LAD);
		Name = 'Single Track + Lane Tracking Model';
	
	case {'stvis_1int'} % Einspurmdl. + lane tracking + 1fach int. bzgl. yL
        [sys,UserData] = STM_LT_yL1int(paramsVhcl,vx,LAD);
		Name = 'Single Track + Lane Tracking + single integral action on lateral offset';
        
    case {'stvis_2int'} % Einspurmdl. + lane tracking + 2fach int. bzgl. yL
        [sys,UserData] = STM_LT_yL2int(paramsVhcl,vx,LAD);
		Name = 'Single Track + Lane Tracking + double integral action on lateral offset';
		
	
    %%% mit Lenkmodell %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    case {'stdsr'} % Einspurmdl. + Lenkmodell CarMaker DSR
        [sys,UserData] = STM_DSR(paramsVhcl,vx,paramsSteer);
		Name = 'Single Track + Steering Model CarMaker DSR';
    
    case {'stvisdsr'} % Einspurmdl. + lane tracking + Lenkmdl. CarMaker DSR
        [sys,UserData] = STM_LT_DSR(paramsVhcl,vx,LAD,paramsSteer);
		Name = 'Single Track + Lane Tracking + Steering Model CarMaker DSR';
    
	case {'stvisdsr_1int'}
        [sys,UserData] = STM_LT_DSR_yL1int(paramsVhcl,vx,LAD,paramsSteer);
		Name = 'Single Track + Lane Tracking + Steering Model CarMaker DSR + single integral action on lateral offset';
	
    case {'stvisdsr_2int'}
        [sys,UserData] = STM_LT_DSR_yL2int(paramsVhcl,vx,LAD,paramsSteer);
		Name = 'Single Track + Lane Tracking + Steering Model CarMaker DSR + double integral action on lateral offset';
    
	
    %%% otherwise %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
    otherwise
        error('Unknown string sw');
        
end%switch


% % create state space model
% sys = ss(A,B,C,D);
% 
% % set input/state/output names
% sys.InputName = InDesc(1,:);
% sys.InputUnit = InDesc(2,:);
% try
% 	% In Matlab R2012a this does not work!
% 	sys.StateName = StateDesc(1,:);
% 	sys.StateUnit = StateDesc(2,:);
% catch exc
% 	warning('This MATLAB release does not support properties "StateName" or "StateUnit" for class GENSS.');
% end
% sys.OutputName = OutDesc(1,:);
% sys.OutputUnit = OutDesc(2,:);

% set info
sys.Name		= Name;
sys.UserData	= UserData;
% sys.UserData.params_Vehicle = paramsVhcl;
% sys.UserData.params_Steer	= paramFile2Struct(pFile_Steer);

end%fcn


function [sys,UserData] = STM_LT(params,vx,LAD)
%  State-space model combined single track & lane tracking model.
%   SYS = singleTrack_lt(paramFile,vx,lad)
%   
%   Inputs:
%   paramFile ... m-file filename containing vehicle parameters
%   vx .......... longitudinal velocity
%   lad ......... look-ahead distance
%
%   State variables: (Index V...Vehicle, g...global)
%   x1 = lateral velocity (vy_V)
%   x2 = yaw rate (psiDot)
%   x3 = lateral offset @ LAD (latOff)
%   x4 = angular deviation @ LAD (angDev)
% 
% Source: A Comparative Study of Vision-Based Lateral Control Strategies 
% for Autonomous Highway Driving.
% 

% % load parameter
% [csv,csh,lv,lh,m,Iz] = getParams_STM(params);
% % set state space elements
% A = [-(csh+csv)/(m*vx), (csh*lh-csv*lv)/(m*vx) - vx, 0, 0;...
%     (csh*lh-csv*lv)/(Iz*vx), -(csh*lh^2+csv*lv^2)/(Iz*vx), 0, 0;    
%     -1, -LAD, 0, vx;...
%     0, -1, 0, 0];
% B = [[csv/m;csv*lv/Iz;0;0], [0;0;0;vx]];

% load single track model
stm = ssMdl_SingleTrack('stm', params, vx);

% load lane tracking model
ltm = ssMdl_laneTracking('ltm', vx, LAD);

% series connection of STM -> LTM
sys = xperm(connect(ltm, stm, ...
	{'front wheel angle' 'curv(LAD)'}, ...
	{'v_y' 'yaw rate' 'latOff(LAD)' 'angDev(LAD)'}), ...
	[3 4 1 2]);


% info
UserData.vx		= stm.UserData.vx;
UserData.params = stm.UserData.params;
UserData.LAD	= ltm.UserData.LAD;

end%fcn

function [sys,UserData] = STM_LT_yL1int(params,vx,LAD)
% State-space model of combined single track, lane tracking model &
% integral action on lateral offset.
% 
%   SYS = singleTrack_lt_yL1int(paramFile,vx,lad)
% 
%   Inputs:
%   paramFile ... m-file filename containing vehicle parameters
%   vx .......... longitudinal velocity
%   lad ......... look-ahead distance
% 
%   State variables: (Index V...Vehicle, g...global) 
%   x1 = lateral velocity (vy_V)
%   x2 = yaw rate (psiDot)
%   x3 = lateral offset @ LAD (latOff)
%   x4 = angular deviation @ LAD (angDev)
%   x5 = Int{yL} (*)
% 
% Source: A Comparative Study of Vision-Based Lateral Control Strategies 
% for Autonomous Highway Driving.
% Source int. Model: Dynamic Controller for Lane Keeping and Obstacle 
% Avoidance Assistance System.
% 

% % load parameter
% [csv,csh,lv,lh,m,Iz] = getParams_STM(params);
% % set state space elements
% A = [-(csh+csv)/(m*vx),(csh*lh-csv*lv)/(m*vx) - vx,0,0,0;...
%     (csh*lh-csv*lv)/(Iz*vx),-(csh*lh^2+csv*lv^2)/(Iz*vx),0,0,0;    
%     -1,-LAD,0,vx,0;...
%     0,-1,0,0,0;...
%     0,0,1,0,0];
% B = [[csv/m;csv*lv/Iz;0;0;0], [0;0;0;vx;0]];

sys = getMergedSSMdl('stvis', params, vx, LAD);
UserData = sys.UserData;

intsys = ss(0,1,1,0,...
	'InputName','latOff(LAD)',...
	'StateName','Int{latOff}',...
	'OutputName','Int{latOff}');
sys = connect(sys, intsys, sys.u, [sys.y; intsys.y]);

end%fcn

function [sys,UserData] = STM_LT_yL2int(params,vx,LAD)
% State-space model of combined single track, lane tracking model & double
% integral action on lateral offset.
%   ret = singleTrack_lt_yL2int(paramFile,vx,lad)
% 
%   Inputs:
%   paramFile ... m-file filename containing vehicle parameters
%   vx .......... longitudinal velocity
%   lad ......... look-ahead distance
% 
%   State variables: (Index V...Vehicle, g...global) 
%   x1 = lateral velocity (vy_V)
%   x2 = yaw rate (psiDot)
%   x3 = lateral offset @ LAD (latOff)
%   x4 = angular deviation @ LAD (angDev)
%   x5 = Int{Int{yL}}
%   x6 = Int{yL}
% 
% Source: A Comparative Study of Vision-Based Lateral Control Strategies 
% for Autonomous Highway Driving.
% Source int. Model: Dynamic Controller for Lane Keeping and Obstacle 
% Avoidance Assistance System.
% 

% % load parameter
% [csv,csh,lv,lh,m,Iz] = getParams_STM(params);
% % set state space elements
% A = [-(csh+csv)/(m*vx),(csh*lh-csv*lv)/(m*vx) - vx,0,0,0,0;...
%     (csh*lh-csv*lv)/(Iz*vx),-(csh*lh^2+csv*lv^2)/(Iz*vx),0,0,0,0;    
%     -1,-LAD,0,vx,0,0;...
%     0,-1,0,0,0,0;...
%     0,0,0,0,0,1;...
%     0,0,1,0,0,0];
% B = [[csv/m;csv*lv/Iz;0;0;0;0], [0;0;0;vx;0;0]];

sys = getMergedSSMdl('stvis_1int', params, vx, LAD);
UserData = sys.UserData;

intsys = ss(0,1,1,0,...
	'InputName','Int{latOff}',...
	'StateName','Int2{latOff}',...
	'OutputName','Int2{latOff}');

% make sure outputs are in the desired order
sys = connect(sys, intsys, sys.u, [sys.y(1:4); intsys.y; sys.y(5)]);

% make sure states are in the desired order
sys = xperm(sys, [1:4,6,5]);

end%fcn



function [sys,UserData] = STM_DSR(paramsVhcl,vx,paramsSteer)
% State-space model of combined single track model & steering model DSR
%   SYS = singleTrack_DSR(paramFile,vx,paramFileSteer)
% 
%   Input arguments:
%   paramFile ........ m-file filename containing vehicle parameters
%   vx ............... longitudinal velocity
%   paramFileSteer ... m-file filename containing steering parameters
% 
% 	State variables:
%   x1 = lateral velocity (v_y)
%   x2 = yaw rate
%	x3 = steering wheel angle
%	x4 = steering wheel angular velocity
% 
% Based on: 'A Comparative Study of Vision-Based Lateral Control Strategies 
% for Autonomous Highway Driving' and CarMaker Manual (see 'Dynamic Steer
% Ratio').
% 

% % load parameter
% [csv,csh,lv,lh,m,Iz] = getParams_STM(paramsVhcl);
[J,V,alph,drack,drot,iHR,mL,mR,mr,xi] = getParams_CarMakerDSR(paramsSteer);
% 
% % set state space elements
% A = [-(csh+csv)/(m*vx), (csh*lh-csv*lv)/(m*vx) - vx;...
%      (csh*lh-csv*lv)/(Iz*vx), -(csh*lh^2+csv*lv^2)/(Iz*vx)];
% B = [csv/m; csv*lv/Iz];

% load single track model
stm = ssMdl_SingleTrack('stm', paramsVhcl, vx);

% Create series connection with steering model manually
A = [stm.A,stm.B/alph,[0;0];...
    0,0,0,1;...
    0,0,0,-1/xi*(drot*iHR^2+drack)];
B = [...
	0 0 0;
	0 0 0;
	0 0 0;
	iHR^2*V/xi, iHR/xi, -iHR/xi,...
	];
C = eye(size(A));
D = 0;


% set state/input/output names
StateDesc = [...
	stm.StateName','SWAngle','SWVelocity';...
	{'m/s','rad/s','rad','rad/s'}];
InputDesc = {...
	'SWTorque' 'F_tieRod_left' 'F_tieRod_right';
	'Nm' 'N' 'N'};
OutputDesc = StateDesc;


%%% create state space model
sys = ss(A,B,C,D);

% set input/state/output names
sys.InputName = InputDesc(1,:);
sys.InputUnit = InputDesc(2,:);
try
	% In Matlab R2012a this does not work!
	sys.StateName = StateDesc(1,:);
	sys.StateUnit = StateDesc(2,:);
catch exc
	warning('This MATLAB release does not support properties "StateName" or "StateUnit" for class GENSS.');
end
sys.OutputName = OutputDesc(1,:);
sys.OutputUnit = OutputDesc(2,:);

% info
UserData = stm.UserData;
UserData.paramsSteer = paramsSteer;

end%fcn

function [sys,UserData] = STM_LT_DSR(paramsVhcl,vx,LAD,paramsSteer)
% State space model of combined single track model, lane tracking model &
% steering model DSR.
%   SYS = singleTrack_lt_DSR(paramFile,vx,lad,paramFileSteer)
% 

% % load parameter
% [csv,csh,lv,lh,m,Iz] = getParams_STM(paramsVhcl);
% [J,V,alph,drack,drot,iHR,mL,mR,mr,xi] = getParams_CarMakerDSR(paramsSteer);
% % set state space elements
% A = [-(csh+csv)/(m*vx),(csh*lh-csv*lv)/(m*vx) - vx,0,0;...
%     (csh*lh-csv*lv)/(Iz*vx),-(csh*lh^2+csv*lv^2)/(Iz*vx),0,0;    
%     -1,-LAD,0,vx;...
%     0,-1,0,0];
% B = [csv/m; csv*lv/Iz; 0; 0];
% 
% A = [A,B/alph,[0;0;0;0];...
%     0,0,0,0,0,1;...
%     0,0,0,0,0,-1/xi*(drot*iHR^2+drack)];
% B = [...
% 	[0;0;0;0;0;iHR^2*V/xi],...
% 	[0;0;0;vx;0;0],...
% 	[0;0;0;0;0;iHR/xi],...
% 	[0;0;0;0;0;-iHR/xi],...
% 	];

% load single track model with steering model
stm = getMergedSSMdl('stdsr', paramsVhcl, vx, LAD, paramsSteer);

% load lane tracking model
ltm = ssMdl_laneTracking('ltm', vx, LAD);

% series connection of STM -> LTM
sys = xperm(connect(ltm, stm, ...
	{'SWTorque' 'curv(LAD)' 'F_tieRod_left' 'F_tieRod_right'}, ...
	{'v_y' 'yaw rate' 'latOff(LAD)' 'angDev(LAD)' 'SWAngle' 'SWVelocity'}), ...
	[3 4 1 2 5 6]);

% set state/input/output names
InputDesc = {...
	'SWTorque','road curvature at LAD','left tie rod force','right tie rod force';
	'Nm','1/m','N','N'};

% info
UserData.vx			= stm.UserData.vx;
UserData.paramsVhcl = stm.UserData.params;
UserData.paramsSteer= stm.UserData.paramsSteer;
UserData.LAD		= ltm.UserData.LAD;

end%fcn

function [sys,UserData] = STM_LT_DSR_yL1int(paramsVhcl,vx,LAD,paramsSteer)
% State-space model of combined single track modell, lane tracking model,
% steering model & integral action on lateral offset.
%   SYS = singleTrack_lt_DSR_yL1int(paramFile,vx,lad)
% 
% Source int. Model: 'Dynamic Controller for Lane Keeping and Obstacle 
% Avoidance Assistance System'.
% 

% % load parameter
% [csv,csh,lv,lh,m,Iz] = getParams_STM(paramsVhcl);
% [J,V,alph,drack,drot,iHR,mL,mR,mr,xi] = getParams_CarMakerDSR(paramsSteer);
% % set state space elements
% A = [-(csh+csv)/(m*vx), (csh*lh-csv*lv)/(m*vx) - vx, 0, 0;...
%     (csh*lh-csv*lv)/(Iz*vx), -(csh*lh^2+csv*lv^2)/(Iz*vx), 0, 0;    
%     -1, -LAD, 0, vx;...
%     0 -1 0 0];
% B = [csv/m; csv*lv/Iz; 0; 0];
% 
% A = [A,B/alph,[0;0;0;0],zeros(4,1);...
%     0 0 0 0 0 1 0;...
%     0,0,0,0,0,-1/xi*(drot*iHR^2+drack),0;...
%     0 0 1 0 0 0 0];
% B = [...
% 	[0;0;0;0;0;iHR^2*V/xi;0],...
% 	[0;0;0;vx;0;0;0],...
% 	[0;0;0;0;0;iHR/xi;0],...
% 	[0;0;0;0;0;-iHR/xi;0],...
% 	];

sys = getMergedSSMdl('stvisdsr', paramsVhcl, vx, LAD, paramsSteer);
UserData = sys.UserData;

intsys = ss(0,1,1,0,...
	'InputName','latOff(LAD)',...
	'StateName','Int{latOff}',...
	'OutputName','Int{latOff}');
sys = connect(sys, intsys, sys.u, [sys.y; intsys.y]);

end%fcn

function [sys,UserData] = STM_LT_DSR_yL2int(paramsVhcl,vx,LAD,paramsSteer)
% State-space model of combined single track modell, lane tracking model,
% steering model & double integral action on lateral offset.
%   SYS = singleTrack_lt_DSR_yL2int(paramFile,vx,lad)
% 
% Source int. Model: 'Dynamic Controller for Lane Keeping and Obstacle 
% Avoidance Assistance System'.
% 

% % load parameter
% [csv,csh,lv,lh,m,Iz] = getParams_STM(paramsVhcl);
% [J,V,alph,drack,drot,iHR,mL,mR,mr,xi] = getParams_CarMakerDSR(paramsSteer);
% % set state space elements
% A = [-(csh+csv)/(m*vx), (csh*lh-csv*lv)/(m*vx) - vx, 0, 0;...
%     (csh*lh-csv*lv)/(Iz*vx), -(csh*lh^2+csv*lv^2)/(Iz*vx), 0, 0;    
%     -1, -LAD, 0, vx;...
%     0 -1 0 0];
% B = [csv/m; csv*lv/Iz; 0; 0];
% 
% A = [A,B/alph,[0;0;0;0],zeros(4,2);...
%     0 0 0 0 0 1 0 0;...
%     0,0,0,0,0,-1/xi*(drot*iHR^2+drack),0,0;...
%     0 0 0 0 0 0 0 1;...
%     0 0 1 0 0 0 0 0];
% B = [...
% 	[0;0;0;0;0;iHR^2*V/xi;0;0],...
% 	[0;0;0;vx;0;0;0;0],...
% 	[0;0;0;0;0;iHR/xi;0;0],...
% 	[0;0;0;0;0;-iHR/xi;0;0],...
% 	];


sys = getMergedSSMdl('stvisdsr_1int', paramsVhcl, vx, LAD, paramsSteer);
UserData = sys.UserData;

intsys = ss(0,1,1,0,...
	'InputName','Int{latOff}',...
	'StateName','Int2{latOff}',...
	'OutputName','Int2{latOff}');
% make sure outputs are in the desired order
sys = connect(sys, intsys, sys.u, [sys.y(1:6); intsys.y; sys.y(7)]);

% make sure states are in the desired order
sys = xperm(sys, [1:6,8,7]);

end%fcn


% --- Get parameters from structure ------------------------------------- %
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
