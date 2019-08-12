function sys = ssMdl_laneTracking(sw,vx,LAD)
% SSMDL_LANETRACKING	Lane tracking state-space model.
%	SYS = SSMDL_LANETRACKING(SW,VX,LAD) returns the lane tracking model
%	representation SW using longitudinal velocity VX and look-ahead
%	distance LAD as a state-space model SYS.
% 
%	Source: A Comparative Study of Vision-Based Lateral Control Strategies
%	for Autonomous Highway Driving. Kosecka, 1998.

% Subject: Master's Thesis - LKA
% $Author$
% $LastChangedDate$
% $Revision$


%%% check input arguments
narginchk(0,3);

% tunable parameter
if nargin < 3 || isempty(LAD)
	LAD = realp('LAD',0);
end%if
if nargin < 2 || isempty(vx)
	vx	= realp('vx',10);
end%if
if nargin < 1
	sw = 'LTM';
end%if


%%% get state space data
switch upper(sw)
	case 'LTM'
		 [A,B,C,D,InDesc,StateDesc,OutDesc,UD] = LTM(vx, LAD);
		 Name = 'Lane Tracking Model';
		
	otherwise
		error('Unknown case!');
end%fcn



%%% create state space model
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

% set additional info
sys.Name = Name;
sys.UserData = UD;

end%fcn


function [A,B,C,D,InputDesc,StateDesc,OutputDesc,UserData] = LTM(vx, LAD)

narginchk(2,2);

A = [0 vx; 0 0];
B = [-1 -LAD 0; 0 -1 vx];
C = eye(size(A));
D = zeros(size(B));

% set state/input/output names
StateDesc = {...
	'latOff(LAD)','angDev(LAD)';
	'm','rad'};
InputDesc = {...
	'v_y','yaw rate','curv(LAD)';
	'm/s','rad/s','1/m'};
OutputDesc = StateDesc;

% User Data
UserData.LAD.about = 'look-ahead distance';
UserData.LAD.value = LAD;
UserData.LAD.unit = 'm';

end%fcn
