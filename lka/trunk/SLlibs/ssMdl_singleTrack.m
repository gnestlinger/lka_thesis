function [sys] = ssMdl_singleTrack(sw,paramsVhcl,vx)
% SSMDL_SINGLETRACK		Single-track state-space model.
%   SYS = SSMDL_SINGLETRACK(SW,PARAMSVHCL,VX,LAD,varargin)
%   ________________
%   Input arguments:
%   SW ............ string to choose state space model (st/stvis/. see below)
%   PARAMSVHCL .... Struct of vehicle parameters
%   VX ............ longitudinal velocity [m/s]
%   LAD ........... look-ahead distance [m]
%   PFILE_STEER ... opt. input arguments (e.g. additional parameter files)
%   ________________________________________________
% 	SW can take the following values:
%     - 'stm'
% 
% Source: see subfunctions

% Subject: Diplomarbeit - LKA
% $Author$
% $LastChangedDate$
% $Revision$


%%% handle input arguments
% tunable parameter
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
if ~isscalar(vx)
	error('Input argument VX must be scalar!');
end%if



% get state space data
switch lower(sw)
    
    case {'stm'} % Linear dynamic single track model
        [A,B,C,D,InDesc,StateDesc,OutDesc,UserData] = STM(paramsVhcl,vx);
		Name = 'Single Track Model';
	
    otherwise
        error('Unknown string SW');
        
end%switch


% create state space model
sys = ss(A,B,C,D);

% set input/state/output names
sys.InputName = InDesc(1,:);
sys.InputUnit = InDesc(2,:);
try
	% GENSS objects don't have properties STATENAME/STATEUNIT before R2017a!
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

end%fcn


function [A,B,C,D,InputDesc,StateDesc,OutputDesc,UD] = STM(params,vx)
% singleTrack     state-space model of single track model
%	[A,B,C,D,INPUTDESC,STATEDESC,OUTPUTDESC,UD] = STM(PARAMS,VX)
%	
% 	State variables:
%   x1 = Fahrzeugquergeschwindigkeit (v_y)
%   x2 = Giergeschwindigkeit (psiDot)
% 
% Source: A Comparative Study of Vision-Based Lateral Control Strategies 
% for Autonomous Highway Driving.
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
	'v_y','yaw rate';
	'm/s','rad/s'};
InputDesc = {...
	'front wheel angle';
	'rad'};
OutputDesc = StateDesc;

% info
UD.vx.about = 'longitudinal velocity';
UD.vx.value = vx;
UD.vx.unit	= 'm/s';
UD.params	= params;

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
