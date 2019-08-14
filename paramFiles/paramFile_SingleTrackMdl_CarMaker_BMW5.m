% PARAMETER FILE
% 
% Single track model
% BMW 5
% 
% Source: IPG Carmaker 4.0
% 

% Subject: Parameter file
% $Author$
% $LastChangedDate$
% $Revision$


about.vehicle = 'BMW 5';
about.source = 'IPG Carmaker 4.0';

m = 1564;		% Fahrzeugmasse [kg]
Izz = 2230;		% Trägheitsmoment um die Hochachse [kg*m^2]

l_front = 1.268;		% Schwerpunktlage/Vorderachse [m]
l_rear = 1.620;		% Schwerpunktlage/Hinterachse [m]
l = l_front+l_rear;		% Achsabstand [m]

cs_front = 2*70000;	% Schräglaufsteifigkeit vorne [N/rad]
cs_rear = 2*70000;	% Schräglaufsteifigkeit hinten [N/rad]

