% Simulation CarMaker-Lenkmodell: Simulink
% 
% Associated Simulink-Models: 
%   .) Lenkung_CarMaker_DSR.mdl
% 
% Subject: lka
% Author: georgnoname
% Date: 11.12.2012 - 17.03.2013

clc
clear all
% close all

% load Simulink model
mdlName = 'Lenkung_CarMaker_DSR';
load_system(mdlName);


%% steering model controller design

% parameter file
stringC_steer = 'paramFile_SteeringMdl_CarMaker_DSR';

% state-space model
sysSteeringModel = ssMdl_EPS('CarMakerDSR',stringC_steer);

% transfer function
P = tf(sysSteeringModel);

% Wunsch: Überschwingen (Mp) und Anstiegszeit (tr)
Mp = 1.0;
tr = 0.51;

% Tabelle (vgl. RT Horn/Dourdoumas S. 347)
xi = interp1([1.00;1.02;1.05;1.10;1.16;1.25;1.37;1.53;1.73],...
	[0.9,0.8,0.7,0.6,0.5,0.4,0.3,0.2,0.1],Mp,'spline','extrap');
wntr = interp1([1.00;1.02;1.05;1.10;1.16;1.25;1.37;1.53;1.73],...
    [2.85,2.55,2.28,2.04,1.84,1.66,1.50,1.37,1.26],Mp,'spline','extrap');

% Soll-Führungsübertragungsfunktion (System mit dominantem Polpaar)
s = tf('s');
Tsoll = (wntr/tr)^2/(s^2 + 2*xi*wntr/tr*s + (wntr/tr)^2);

% Ordnungserhöhung von Tsoll
spec = [0];
ws = s+10; % Hurwitz-Polynom!
Tsoll = Tsoll*(ws/ws)^(1+length(spec));

% steering model: controller design (algebraische Synthese)
pin.Controller.steerSys.desc = 'Algebraische Synthese';
[R,V] = algsynth(P,Tsoll,spec);
pin.Controller.steerSys.R.tf = R;
pin.Controller.steerSys.R.ss = ssNormForm(R,'obs');
pin.Controller.steerSys.V.tf = V;
pin.Controller.steerSys.V.ss = ssNormForm(V,'obs');
pin.Controller.steerSys.initialValue = ...
    zeros(length(pin.Controller.steerSys.R.ss.a),1);

Rpol = polvor(P,Tsoll.den{1},spec);


%% Simulation parameter

% simulation: steering parameter
stringSim = 'paramFile_SteeringMdl_CarMaker_DSR';
pin.SteeringModel.parameterFile = stringSim;
pin.SteeringModel.parameter = paramFile2Struct(stringSim);


