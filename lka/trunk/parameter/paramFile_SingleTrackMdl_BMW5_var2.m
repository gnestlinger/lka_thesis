% PARAMETER FILE
% 
% Single track model
% BMW 5 - Scenario 2
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
Iz = 2230;		% Tr�gheitsmoment um die Hochachse [kg*m^2]

lv = 1.268;		% Schwerpunktlage/Vorderachse [m]
lh = 1.620;		% Schwerpunktlage/Hinterachse [m]
l = lv+lh;		% Achsabstand [m]

csv = 2*70000;	% Schr�glaufsteifigkeit vorne [N/rad]
csh = 2*70000;	% Schr�glaufsteifigkeit hinten [N/rad]


%%% Parametervariation

% Schwerpunkt (l=0)
% O------x--------O
% |<-----|------->|
% lv     0      -lh

%%% Zuladung 
mzu = [0,0];
% bei Position
lzu = [-lh,lv];

% neue Gesamtmasse
m = m + sum(mzu);
% Schwerpunktverschiebung (vgl. de.wikipedia.org/wiki/Massenmittelpunkt)
S = 1/m*dot(mzu,lzu);

% neue Schwerpunktlage
lv = lv - S;
lh = lh + S;

% neues Tr�gheitsmoment
Iz = Iz + dot(mzu,lzu.^2);

%%% Schr�glaufsteifigkeiten

csv = 1.15*csv;
csh = 0.85*csh;

