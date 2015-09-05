% PARAMETER FILE
% 
% Single track model
% BMW 5
% 
% Source: IPG Carmaker 4.0
% 

% Subject: 
% $Author: georgnoname@gmail.com $
% $LastChangedDate: 2015-08-12 11:19:04 +0200 (Mi, 12 Aug 2015) $
% $Revision: 196 $


about.vehicle = 'BMW 5';
about.source = 'IPG Carmaker 4.0';

m = 1564;       % Fahrzeugmasse [kg]
Iz = 2230;      % Trägheitsmoment um die Hochachse [kgm^2]

lv = 1.268;     % Schwerpunktlage [m]
lh = 1.620;     % Schwerpunktlage [m]
l = lv+lh;      % Achsabstand [m]
                
csv = 2*70000;  % Schräglaufsteifigkeit vorne [N/rad]
csh = 2*70000;  % Schräglaufsteifigkeit hinten [N/rad]

