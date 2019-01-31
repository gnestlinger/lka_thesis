% PARAMETER FILE
% 
% Single track model
% BMW 5 - Scenario 1
% 
% Source: IPG Carmaker 4.0
% 

% Subject: Parameter file
% $Author$
% $LastChangedDate$
% $Revision$


about.vehicle = 'BMW 5 - Variation';
about.source = 'IPG Carmaker 4.0';

m = 1564;		% Fahrzeugmasse [kg]
Izz = 2230;		% Trägheitsmoment um die Hochachse [kg*m^2]

l_front = 1.268;		% Schwerpunktlage/Vorderachse [m]
l_rear = 1.620;		% Schwerpunktlage/Hinterachse [m]
l = l_front+l_rear;		% Achsabstand [m]

cs_front = 2*70000;	% Schräglaufsteifigkeit vorne [N/rad]
cs_rear = 2*70000;	% Schräglaufsteifigkeit hinten [N/rad]


%%% Parametervariation

% Schwerpunkt (l=0)
% O------x--------O
% |<-----|------->|
% lv     0      -lh

%%% Zuladung 
mzu = [400,0];
% bei Position
lzu = [-l_rear,l_front];

% neue Gesamtmasse
m = m + sum(mzu);
% Schwerpunktverschiebung (vgl. de.wikipedia.org/wiki/Massenmittelpunkt)
S = 1/m*dot(mzu,lzu);

% neue Schwerpunktlage
l_front = l_front - S;
l_rear = l_rear + S;

% neues Trägheitsmoment
Izz = Izz + dot(mzu,lzu.^2);

