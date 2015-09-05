% PARAMETER FILE
% 
% Single track model
% MCC-Smart
% 
% Source: VO Modellbildung und Simulation in der Fahrzeugtechnik 2012
% 

% Subject: 
% $Author: georgnoname@gmail.com $
% $LastChangedDate: 2015-08-12 11:19:04 +0200 (Mi, 12 Aug 2015) $
% $Revision: 196 $


about.vehicle = 'MCC-Smart';
about.source = 'FTG-LV';

m = 820;        % Fahrzeugmasse [kg]
Iz = 3000;      % Tr�gheitsmoment um die Hochachse [kgm^2]

lv = 1.142;     % Schwerpunktlage [m]
lh = 0.670;     % Schwerpunktlage [m]
l = lv+lh;      % Achsabstand [m]
                
csv = 70000;    % Schr�glaufsteifigkeit vorne [N/rad]
csh = 90000;    % Schr�glaufsteifigkeit hinten [N/rad]

vmax = 150/3.6; % maximale Geschwindigkeit [m/s]

% is = 15.3;      % Lenk�bersetzung
% 
% b = 1.5;        % Spurweite [m]
% h = 0.56;       % Schwerpunkth�he [m]
% hsp = 0.46;
% h0 = 0.1;
% 
% Ix = 730;       % Tr�gheitsmoment [kgm^2]
% cx = 70000;     % Wankfedersteifigkeit [Nm/rad]
% dx = 6000;      % Wankd�mpfungskonstante [Nms/rad]
