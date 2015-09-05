function sys = ssMdl_SingleTrack(sw,paramFile,vx,lad,varargin)
% ssMdl_SingleTrack     returns single-track-based state-space models
%   _______
%   Syntax:
%   sys = ssMdl_SingleTrack(sw,paramFile,vx,lad,varargin)
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
% Subject: lka
% Author: georgnoname
% Date: 29.11.2012 - 19.03.2013


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
        sys = singleTrack_lt(paramFile,vx,lad);
        
    case {'stvis_2int'} % Einspurmdl. + lane tracking + 2fach int. bzgl. yL
        sys = singleTrack_lt_yL2int(paramFile,vx,lad);
        
    %%% mit Lenkmodell %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    case {'stdsr'} % Einspurmdl. + Lenkmodell CarMaker DSR
        sys = singleTrack_DSR(paramFile,vx,varargin{1});
    
    case {'stvisdsr'} % Einspurmdl. + lane tracking + Lenkmdl. CarMaker DSR
        sys = singleTrack_lt_DSR(paramFile,vx,lad,varargin{1});
        
    case {'stvisdsr_2int'}
        sys = singleTrack_lt_DSR_yL2int(paramFile,vx,lad,varargin{1});
    
    %%% otherwise %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
    otherwise
        error('Unknown string sw');
        
end%switch


end%fcn



function ret = singleTrack(paramFile,vx)
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
% Author: georgnoname
% Date: 24.10.2012 - 15.03.2013


% load parameter
eval(paramFile);

% init state space model
ret = ss([],[],[],[]);

% set state space elements
% ret.A = ...
%     [-(csh+csv)/(m*vx), (csh*lh-csv*lv)/(m*vx) - vx;...
%     (csh*lh-csv*lv)/(Iz*vx), -(csh*lh^2+csv*lv^2)/(Iz*vx)];
% ret.B = [csv/m; csv*lv/Iz];
ret.A = [0 1 0 0;...
    0, -(csh+csv)/(m*vx), 0, (csh*lh-csv*lv)/(m*vx) - vx;...
    0 0 0 1;...
    0, (csh*lh-csv*lv)/(Iz*vx), 0, -(csh*lh^2+csv*lv^2)/(Iz*vx)];
ret.B = [0; csv/m; 0; csv*lv/Iz];
ret.C = eye(size(ret.A));
ret.D = 0;

% set state/input/output names
ret.StateName = {'sy','vy','psi','psiDot'};
ret.Inputname = {'steer angle'};
% ret.OutputName = {};


% info
ret.userdata.about = 'Einspurmodell';
ret.userdata.vx.about = 'Längsgeschwindigkeit';
ret.userdata.vx.value = vx;
ret.userdata.vx.unit = 'm/s';

end%fcn



function ret = singleTrack_lt(paramFile,vx,lad)
% singleTrack_lt    state-space model of single track modell + 
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
% Author: georgnoname
% Date: 24.10.2012 - 15.03.2013


% load parameter
eval(paramFile);

% init state space model
ret = ss([],[],[],[]);

% set state space elements
ret.A = [-(csh+csv)/(m*vx), (csh*lh-csv*lv)/(m*vx) - vx, 0, 0;...
    (csh*lh-csv*lv)/(Iz*vx), -(csh*lh^2+csv*lv^2)/(Iz*vx), 0, 0;    
    -1, -lad, 0, vx;...
    0, -1, 0, 0];
ret.B = [csv/m; csv*lv/Iz; 0; 0];
% ret.C = [0 1 0 0; 0 0 1 0; 0 0 0 1];
ret.C = [0 0 1 0];
% ret.C = [0 1 0 0; 0 0 1 0];
ret.D = 0;

% set state/input/output names
ret.StateName = {'vy','psiDot','yL','epsL'};
ret.Inputname = {'steer angle'};
% ret.OutputName = {};


% info
ret.userdata.about = 'Einspurmodell + Relativposition';
ret.userdata.vx.about = 'Längsgeschwindigkeit';
ret.userdata.vx.value = vx;
ret.userdata.vx.unit = 'm/s';
ret.userdata.lad.about = 'look-ahead distance';
ret.userdata.lad.value = lad;
ret.userdata.lad.unit = 'm';

end%fcn



function ret = singleTrack_lt_yL2int(paramFile,vx,lad)
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
% Author: georgnoname
% Date: 29.11.2012 - 15.03.2013


% load parameter
eval(paramFile);

% init state space model
ret = ss([],[],[],[]);

% set state space elements
ret.A = [-(csh+csv)/(m*vx),(csh*lh-csv*lv)/(m*vx) - vx,0,0,0,0;...
    (csh*lh-csv*lv)/(Iz*vx),-(csh*lh^2+csv*lv^2)/(Iz*vx),0,0,0,0;    
    -1,-lad,0,vx,0,0;...
    0,-1,0,0,0,0;...
    0,0,0,0,0,1;...
    0,0,1,0,0,0];
ret.B = [csv/m; csv*lv/Iz; 0; 0; 0; 0];
ret.C = zeros(1,length(ret.A));
ret.C(3) = 1;
ret.D = 0;

% set state/input/output names
ret.StateName = {'vy','psiDot','yL','epsL','IntInt{yL}','Int{yL}'};
ret.Inputname = {'steer angle'};
% ret.OutputName = {};


% info
ret.userdata.about = 'Einspurmodell + Relativposition + 2fach int. bzgl yL';
ret.userdata.vx.about = 'Längsgeschwindigkeit';
ret.userdata.vx.value = vx;
ret.userdata.vx.unit = 'm/s';
ret.userdata.lad.about = 'look-ahead distance';
ret.userdata.lad.value = lad;
ret.userdata.lad.unit = 'm';

end%fcn



function ret = singleTrack_DSR(paramFile,vx,paramFileSteer)
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
% Author: georgnoname
% Date: 29.01.2013 - 19.03.2013


% check input arguments
if ~ischar(paramFileSteer); 
    error('Input argument steering-parameter-filename not of type char'); 
end


% load parameter
eval(paramFile);
eval(paramFileSteer);


% init state space model
ret = ss([],[],[],[]);

% set state space elements
A = [0 1 0 0;...
    0, -(csh+csv)/(m*vx), 0, (csh*lh-csv*lv)/(m*vx) - vx;...
    0 0 0 1;...
    0, (csh*lh-csv*lv)/(Iz*vx), 0, -(csh*lh^2+csv*lv^2)/(Iz*vx)];
B = [0; csv/m; 0; csv*lv/Iz];

ret.A = [A,B/alph,[0;0;0;0];...
    0,0,0,0,0,1;...
    0,0,0,0,0,-1/xi*(drot*iHR^2+drack)];
ret.B = [0; 0; 0; 0; 0; iHR^2*V/xi];
ret.C = zeros(1,length(ret.A));
ret.C(3) = 1;
ret.D = 0;

% set state/input/output names
ret.StateName = {'sy','vy','psi','psiDot','deltaH','deltaHDot'};
ret.Inputname = {'sw torque'};
% ret.OutputName = {};


% info
ret.userdata.about = 'Einspurmodell + Lenkmodell CarMaker DSR';
ret.userdata.vx.about = 'Längsgeschwindigkeit';
ret.userdata.vx.value = vx;
ret.userdata.vx.unit = 'm/s';

end%fcn



function ret = singleTrack_lt_DSR(paramFile,vx,lad,paramFileSteer)
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
% Author: georgnoname
% Date: 29.01.2013 - 19.03.2013


% check input arguments
if ~ischar(paramFileSteer); 
    error('Input argument steering-parameter-filename not of type char'); 
end


% load parameter
eval(paramFile);
eval(paramFileSteer);

% init state space model
ret = ss([],[],[],[]);

% set state space elements
A = [-(csh+csv)/(m*vx),(csh*lh-csv*lv)/(m*vx) - vx,0,0;...
    (csh*lh-csv*lv)/(Iz*vx),-(csh*lh^2+csv*lv^2)/(Iz*vx),0,0;    
    -1,-lad,0,vx;...
    0,-1,0,0];
B = [csv/m; csv*lv/Iz; 0; 0];

ret.A = [A,B/alph,[0;0;0;0];...
    0,0,0,0,0,1;...
    0,0,0,0,0,-1/xi*(drot*iHR^2+drack)];
ret.B = [0; 0; 0; 0; 0; iHR^2*V/xi];
ret.C = zeros(1,length(ret.A));
ret.C(3) = 1;
ret.D = 0;

% set state/input/output names
ret.StateName = {'vy','psiDot','yL','epsL','deltaH','deltaHDot'};
ret.Inputname = {'sw torque'};
% ret.OutputName = {};


% info
ret.userdata.about = 'Einspurmdl + Rel.pos. + Lenkmdl CarMaker DSR';
ret.userdata.vx.about = 'Längsgeschwindigkeit';
ret.userdata.vx.value = vx;
ret.userdata.vx.unit = 'm/s';
ret.userdata.lad.about = 'look-ahead distance';
ret.userdata.lad.value = lad;
ret.userdata.lad.unit = 'm';

end%fcn



function ret = singleTrack_lt_DSR_yL2int(paramFile,vx,lad,paramFileSteer)
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
% Author: georgnoname
% Date: 01.03.2013 - 19.03.2013


% check input arguments
if ~ischar(paramFileSteer); 
    error('Input argument steering-parameter-filename not of type char'); 
end


% load parameter
eval(paramFile);
eval(paramFileSteer);

% init state space model
ret = ss([],[],[],[]);

% set state space elements
A = [-(csh+csv)/(m*vx), (csh*lh-csv*lv)/(m*vx) - vx, 0, 0;...
    (csh*lh-csv*lv)/(Iz*vx), -(csh*lh^2+csv*lv^2)/(Iz*vx), 0, 0;    
    -1, -lad, 0, vx;...
    0 -1 0 0];
B = [csv/m; csv*lv/Iz; 0; 0];

ret.A = [A,B/alph,[0;0;0;0],zeros(4,2);...
    0 0 0 0 0 1 0 0;...
    0,0,0,0,0,-1/xi*(drot*iHR^2+drack),0,0;...
    0 0 0 0 0 0 0 1;...
    0 0 1 0 0 0 0 0];
ret.B = [0; 0; 0; 0; 0; iHR^2*V/xi; 0; 0];
ret.C = zeros(1,length(ret.A));
ret.C(3) = 1;
ret.D = 0;

% set state/input/output names
ret.StateName = {'vy','psiDot','yL','epsL','deltaH','deltaHDot',...
    'IntInt{yL}','Int{yL}'};
ret.Inputname = {'sw torque'};
% ret.OutputName = {};


% info
ret.userdata.about = 'Einspurmdl + Rel.pos. + Lenkmdl CarMaker DSR + 2fach int. bzgl. yL';
ret.userdata.vx.about = 'Längsgeschwindigkeit';
ret.userdata.vx.value = vx;
ret.userdata.vx.unit = 'm/s';
ret.userdata.lad.about = 'look-ahead distance';
ret.userdata.lad.value = lad;
ret.userdata.lad.unit = 'm';

end%fcn
