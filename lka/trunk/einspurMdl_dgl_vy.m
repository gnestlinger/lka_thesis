function XDot = einspurMdl_dgl_vy(t,x,uIn,pIn)
% einspurMdl_dgl_vy     DGL-System des linearen Einspurmodells.
%   _______
%   Syntax:
%   XDot = einspurMdl_dgl_vy(t,x,uIn,pIn)
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
% 
% Source: LV 331.094 Modellbildung und Simulation in der Fahrzeugtechnik.
%
% Subject: lka
% Author: georgnoname
% Date: ??.09.2012 - 13.02.2013


% vehicle parameters
F = pIn.VehicleModel.singleTrack.parameter;

% system matrices
A = [0 1 0 0;...
    0, -(F.cs_rear+F.cs_front)/(F.m*pIn.vx),...
        0, (F.cs_rear*F.l_rear-F.cs_front*F.l_front)/(F.m*pIn.vx)-pIn.vx;...
    0 0 0 1;...
    0, (F.cs_rear*F.l_rear-F.cs_front*F.l_front)/(F.Izz*pIn.vx),...
        0, -(F.cs_rear*F.l_rear^2+F.cs_front*F.l_front^2)/(F.Izz*pIn.vx)];
B = [0; F.cs_front/F.m; 0; F.cs_front*F.l_front/F.Izz];

% dynamic system
u = uIn(t,x);
xDot = A*x(1:4) + B*u;

% erw. Zustandsvektor um globale Fahrzeugposition/zurueckgelegte Wegstrecke
x_gDot = pIn.vx*cos(x(3)) - x(2)*sin(x(3));
y_gDot = pIn.vx*sin(x(3)) + x(2)*cos(x(3));
s_gDot = sqrt( pIn.vx^2 + x(2)^2 );

% evaluated right side of differential equations
XDot = [xDot;x_gDot;y_gDot;s_gDot];

end%fcn
