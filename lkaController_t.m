function ret = lkaController_t(paramFile,vx,lad)
% lkaController_t   LKA Reglerentwurf im Zeitbereich
%   _______
%   Syntax:
%   contr = lkaController_t(paramFile,vx,lad)
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
if ~ischar(paramFile); error('2nd input argument not of type char'); end
if numel(vx) > 1; error('Dimension of input argument vx > 1'); end
if numel(lad) > 1; error('Dimension of input argument lad > 1'); end


% load state space models

% Single Track + "lane tracking"
%   x1 = Fahrzeugquergeschwindigkeit (vy_V)
%   x2 = Giergeschwindigkeit (psiDot)
%   x3 = yL
%   x4 = epsL
sysSingleTrackVis = ssMdl_SingleTrack('stvis',paramFile,vx,lad);

% single track modell + "lane tracking" + internal model
%   State variables: (Index V...Vehicle, 0...global) 
%   x1 = Fahrzeugquergeschwindigkeit (vy_V)
%   x2 = Giergeschwindigkeit (psiDot)
%   x3 = Querablage yL
%   x4 = Relativwinkel epsL
%   x5 = Int{Int{yL}}
%   x6 = Int{yL}
sysSingleTrackVis_2int = ssMdl_SingleTrack('stvis_2int',paramFile,vx,lad);


%%% load parameter %%%
% eval(paramFile);
vehicleParam = paramFile2Struct(paramFile);

%% P-Glied: u = kP*yL

kP = 0.05;

% quick info
ret.P.label = 'P';
ret.P.desc = 'P-Glied';

% design parameter
ret.P.designParam.file = paramFile;
ret.P.designParam.vehicle = vehicleParam;
ret.P.designParam.vx = vx;
ret.P.designParam.lad = lad;

% controller values
ret.P.k = kP;
ret.P.u_fh = @(t,x,yL,epsL) (ret.P.k)*yL;



%% Zustandsregler durch EW-Vorgabe: u = -k'*x

% p = eig(sysSingleTrackVis.A);
% xi = 0.707; wn = 0.989; % 0<xi<1
% sigma = -xi*wn; omega = wn*sqrt(1-xi^2);
% eigCC = -xi*wn + 1i*wn*sqrt(1-xi^2);
% pneu = [p(ind); eigCC; conj(eigCC)];
% liefert gleiches Simulationsergebnis wie EA-ZR mit EW bei [-4,-4] für
% vx=30m/s und lad=6m
pneu = [-2.10347+1i*3.80752,-2.10347-1i*3.80752,-4,-4];
k = acker(sysSingleTrackVis.A,sysSingleTrackVis.b,pneu);
% eig(sysSingleTrackVis.A - sysSingleTrackVis.b*kZR);
% figure
% step(ss((sysSingleTrackVis.A-sysSingleTrackVis.b*kZR),sysSingleTrackVis.B,sysSingleTrackVis.C,0));

% quick info
ret.ZR.label = 'ZR';
ret.ZR.desc = 'Zustandsregler durch Eigenwertvorgabe';

% design parameter
ret.ZR.designParam.file = paramFile;
ret.ZR.designParam.vehicle = vehicleParam;
ret.ZR.designParam.vx = vx;
ret.ZR.designParam.lad = lad;
ret.ZR.designParam.EW = pneu;

% controller values
ret.ZR.k = k;
ret.ZR.u_fh = @(t,x,yL,epsL) (-ret.ZR.k)*[x(2);x(4);yL;epsL];



%% Zustandsregler durch EW-Vorgabe (2fach Int. bzgl yL): u = -k'*x

pneu = [-2.10347+1i*3.80752,-2.10347-1i*3.80752,-4,-4,-1,-1];
k = acker(sysSingleTrackVis_2int.A,sysSingleTrackVis_2int.b,pneu);

% quick info
ret.ZR_2int.label = 'ZR int. Mdl.';
ret.ZR_2int.desc = 'Zustandsregler durch Eigenwertvorgabe';

% design parameter
ret.ZR_2int.designParam.file = paramFile;
ret.ZR_2int.designParam.vehicle = vehicleParam;
ret.ZR_2int.designParam.vx = vx;
ret.ZR_2int.designParam.lad = lad;
ret.ZR_2int.designParam.EW = pneu;

% controller values
ret.ZR_2int.k = k;
ret.ZR_2int.u_fh = @(t,x,yL,epsL) (-ret.ZR_2int.k)*[x(2);x(4);yL;epsL;x(8);x(9)];



%% LQR: u = -k'*x

Q = diag([0 0 1 0]);
R = 10;
[k,~,~] = lqr(sysSingleTrackVis,Q,R);

% quick info
ret.LQR.label = 'LQR';
ret.LQR.desc = 'Linear Quadratic Regulator';

% design parameter
ret.LQR.designParam.file = paramFile;
ret.LQR.designParam.vehicle = vehicleParam;
ret.LQR.designParam.vx = vx;
ret.LQR.designParam.lad = lad;
ret.LQR.designParam.Q = Q;
ret.LQR.designParam.R = R;

% controller values
ret.LQR.k = k;
ret.LQR.u_fh = @(t,x,yL,epsL) (-ret.LQR.k)*[x(2);x(4);yL;epsL];



%% LQR (2fach Int. bzgl yL): u = -k'*x

% Q = diag([0 0 1 1 1 1]);
Q = [0 0 0 0 0 0;...
     0 0 0 0 0 0;...
     0 0 1 0 0 0;...
     0 0 0 0 0 0;...
     0 0 0 0 1 0;...
     0 0 0 0 0 1];
R = 10;
[k,~,~] = lqr(sysSingleTrackVis_2int,Q,R);

% quick info
ret.LQR_2int.label = 'LQR int. Mdl.';
ret.LQR_2int.desc = 'Linear Quadratic Regulator with int. Model';

% design parameter
ret.LQR_2int.designParam.file = paramFile;
ret.LQR_2int.designParam.vehicle = vehicleParam;
ret.LQR_2int.designParam.vx = vx;
ret.LQR_2int.designParam.lad = lad;
ret.LQR_2int.designParam.Q = Q;
ret.LQR_2int.designParam.R = R;

% controller values
ret.LQR_2int.k = k;
ret.LQR_2int.u_fh = @(t,x,yL,epsL)...
    (-ret.LQR_2int.k)*[x(2);x(4);yL;epsL;x(8);x(9)];


end%fcn

