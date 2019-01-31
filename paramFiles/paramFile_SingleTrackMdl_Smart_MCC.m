% PARAMETER FILE
% 
% Single track model
% MCC-Smart
% 
% Source: VO Modellbildung und Simulation in der Fahrzeugtechnik 2012
% 

% Subject: Parameter file
% $Author$
% $LastChangedDate$
% $Revision$


about.vehicle = 'MCC-Smart';
about.source = 'FTG-LV';

m = 820;		% Fahrzeugmasse [kg]
Izz = 3000;		% Tr�gheitsmoment um die Hochachse [kg*m^2]

l_front = 1.142;		% Schwerpunktlage/Vorderachse [m]
l_rear = 0.670;		% Schwerpunktlage/Hinterachse [m]
l = lv+l_rear;		% Achsabstand [m]

cs_front = 70000;	% Schr�glaufsteifigkeit vorne [N/rad]
cs_rear = 90000;	% Schr�glaufsteifigkeit hinten [N/rad]

vmax = 150/3.6;	% maximale Geschwindigkeit [m/s]

% is = 15.3;		% Lenk�bersetzung
%
% b = 1.5;		% Spurweite [m]
% h = 0.56;		% Schwerpunkth�he [m]
% hsp = 0.46;
% h0 = 0.1;
%
% Ix = 730;		% Tr�gheitsmoment [kgm^2]
% cx = 70000;		% Wankfedersteifigkeit [Nm/rad]
% dx = 6000;		% Wankd�mpfungskonstante [Nms/rad]
