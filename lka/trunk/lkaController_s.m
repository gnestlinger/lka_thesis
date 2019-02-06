function ret = lkaController_s(paramFile,vx,lad)
% lkaController_s   LKA Reglerentwurf im Frequenzbereich
%   _______
%   Syntax:
%   ret = lkaController_s(paramFile,vx,lad)
%   ________________
%   Input arguments:
%   paramFile ... m-file filename containing vehicle parameters
%   vx .......... longitudinal velocity [m/s]
%   lad ......... look-ahead distance [m]
% 
% Subject: lka
% Author: georgnoname
% Date: 24.09.2012 - 01.03.2013


% check input arguments
% if ~ischar(paramFile); error('2nd input argument not of type char'); end
% if numel(vx) > 1; error('Dimension of input argument vx > 1'); end
% if numel(lad) > 1; error('Dimension of input argument lad > 1'); end


% load state space models

% Single Track + "lane tracking"
%   x1 = Fahrzeugquergeschwindigkeit (vy_V)
%   x2 = Giergeschwindigkeit (psiDot)
%   x3 = yL
%   x4 = epsL
sysSingleTrackVis = ssMdl_SingleTrack('stvis',paramFile2Struct(paramFile),vx,lad);

% single track modell + "lane tracking" + internal model
%   State variables: (Index V...Vehicle, 0...global) 
%   x1 = Fahrzeugquergeschwindigkeit (vy_V)
%   x2 = Giergeschwindigkeit (psiDot)
%   x3 = Querablage yL
%   x4 = Relativwinkel epsL
%   x5 = Int{Int{yL}}
%   x6 = Int{yL}
% sysSingleTrackInt = ssMdl_SingleTrack('stvis_2int',paramFile,vx,lad);



%%% Reglerentwurf

%%% load parameter %%%
eval(paramFile);
vehicleParam = paramFile2Struct(paramFile);


%% Lead-Lag-Glied

% G = tf(sysSingleTrackVis);
% sisotool('bode',G)
% load('contrFKL');

% quick info
ret.FKL.label = 'FKL';
ret.FKL.desc = 'Frequenzkennlinienverfahren';

% design parameter
ret.FKL.designParam.file = paramFile;
ret.FKL.designParam.vehicle = vehicleParam;
ret.FKL.designParam.vx = vx;
ret.FKL.designParam.lad = lad;

% controller values
ret.FKL.R = zpk(-1.7923,-25.2449,1.9217);
ret.FKL.u = [];


%% Lead-Lag-Glied

% G = tf(sysSingleTrackVis);
% sisotool('bode',G)
% load('contrFKL');

% quick info
ret.FKL2.label = 'FKL';
ret.FKL2.desc = 'Frequenzkennlinienverfahren';

% design parameter
ret.FKL2.designParam.file = paramFile;
ret.FKL2.designParam.vehicle = vehicleParam;
ret.FKL2.designParam.vx = vx;
ret.FKL2.designParam.lad = lad;

% controller values
ret.FKL2.R = zpk(-1.28,[-145.6 -145.6],-981.2128); %für BMW5
ret.FKL2.u = [];


%% Algebraische Synthese

% Tabelle (vgl. RT Horn/Dourdoumas S. 347)
Mp_tab = [1.00;1.02;1.05;1.10;1.16;1.25;1.37;1.53;1.73];
xi_tab = [0.9,0.8,0.7,0.6,0.5,0.4,0.3,0.2,0.1];
wt_tab = [2.85,2.55,2.28,2.04,1.84,1.66,1.50,1.37,1.26];

% Überschwingen (Mp) und Anstiegszeit (tr)
Mp = 1.0;
tr = .1;

% Tabelle: Interpolation
xi = interp1(Mp_tab,xi_tab,Mp,'spline','extrap');
wt = interp1(Mp_tab,wt_tab,Mp,'spline','extrap');

% Soll-Führungsübertragungsfunktion
s = tf('s');
Tsoll = (wt/tr)^2/(s^2 + 2*xi*wt/tr*s + (wt/tr)^2);
% eta1 = 100;
% eta2 = 200;
ws = s+30;
% Tsoll = Tsoll*eta1/(s+eta1)*eta2/(s+eta2);
Tsoll = Tsoll*(ws/ws)^(4+3-2-2+2);

% transfer function (plant)
% sys = ss(sysSingleTrackInt.a,sysSingleTrackInt.b,[0,0,0,0,0,1],sysSingleTrackInt.d);
% P = tf(sys);
P = tf(sysSingleTrackVis);

% quick info
ret.AlgSynth.label = 'AlgSynth';
ret.AlgSynth.desc = 'Algebraische Synthese';

% design parameter
ret.AlgSynth.designParam.file = paramFile;
ret.AlgSynth.designParam.P = P;
ret.AlgSynth.designParam.vehicle = vehicleParam;
ret.AlgSynth.designParam.vx = vx;
ret.AlgSynth.designParam.lad = lad;
ret.AlgSynth.designParam.Mp = Mp;
ret.AlgSynth.designParam.tr = tr;
ret.AlgSynth.designParam.Tsoll = Tsoll;

% controller values
[R,V] = algsynth(P(3,1),Tsoll,[0,0]);
ret.AlgSynth.R.tf = R;
ret.AlgSynth.R.ss = ssNormForm(R,'obs');
ret.AlgSynth.V.tf = V;
ret.AlgSynth.V.ss = ssNormForm(V,'obs');
ret.AlgSynth.u = [];


end%fcn

