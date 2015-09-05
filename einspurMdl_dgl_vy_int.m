function XDot = einspurMdl_dgl_vy_int(t,x,uIn,pIn)
% einspurMdl_dgl_vy_int     DGL-System des linearen Einspurmodells + 
% interes Modell der Störungen bzgl. yL.
%   _______
%   Syntax:
%   XDot = einspurMdl_dgl_vy_int(t,x,uIn,pIn)
%   ________________
%   Input arguments:
%   t ..... time
%   x ..... state
%   uIn ... control input (function handle/anonymus function)
%   pIn ... structure of input parameters
%   ________________________________________________
%   State variables: (Index V...Vehicle, g...global)
%   x1 = Weg in Fahrzeugquerrichtung (sy_V)
%   x2 = Fahrzeugquergeschwindigkeit (vy_V)
%   x3 = Gierwinkel (psi)
%   x4 = Giergeschwindigkeit (psiDot)
%   x5 = glob. Position in x_g-Richtung = Int{vx*cos(psi) - vy*sin(psi)}
%   x6 = glob. Position in y_g-Richtung = Int{vx*sin(psi) + vy*cos(psi)}
%   x7 = glob. zurueckgelegte Wegstrecke = Int{sqrt(vx^2+vy^2)} (s_g)
%   x8 = Int{Int{yL}}
%   x9 = Int{yL}
% 
% Source: LV 331.094 Modellbildung und Simulation in der Fahrzeugtechnik.
% Source Int. Model: "Dynamic Controller for Lane Keeping and Obstacle 
% Avoidance Assistance System".

% Subject: lka
% Author: georgnoname
% Date: ??.09.2012 - 13.02.2013


% vehicle parameters
F = pIn.VehicleModel.singleTrack.parameter.values;

% system matrices
A = [0 1 0 0;...
    0, -(F.csh+F.csv)/(F.m*pIn.vx),...
        0, (F.csh*F.lh-F.csv*F.lv)/(F.m*pIn.vx)-pIn.vx;...
    0 0 0 1;...
    0, (F.csh*F.lh-F.csv*F.lv)/(F.Iz*pIn.vx),...
        0, -(F.csh*F.lh^2+F.csv*F.lv^2)/(F.Iz*pIn.vx)];
B = [0; F.csv/F.m; 0; F.csv*F.lv/F.Iz];

% lane tracking
% [~,yL,epsL,rInv] = visionSystem([x(5),x(6),x(3)],pIn.traj,pIn.lad);
    
% dynamic system
[u,yL] = uIn(t,x);
xDot = A*x(1:4) + B*u;

% erw. Zustandsvektor um globale Fahrzeugposition/zurueckgelegte Wegstrecke
x_gDot = pIn.vx*cos(x(3)) - x(2)*sin(x(3));
y_gDot = pIn.vx*sin(x(3)) + x(2)*cos(x(3));
s_gDot = sqrt( pIn.vx^2 + x(2)^2 );

% erw. Zustandsvektor um internes Modell der Störungen (2-fach
% integrierendes Verhalten bzgl. yL)
% x8Dot = x(end);
% x9Dot = yL;

% evaluated right side of differential equations
XDot = [xDot;x_gDot;y_gDot;s_gDot;x(end);yL];

end%fcn
