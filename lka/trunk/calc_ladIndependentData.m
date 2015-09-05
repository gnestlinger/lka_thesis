function ret = calc_ladIndependentData(t,x,fh)
% calc_ladIndependentData   Calculate internal (independent of lad)
% values from solution of DEQ-System 'einspurMdl_fcn_vy'.
%   _______
%   Syntax:
%   ret = calc_ladIndependentData(solIn)
%   ________________
%   Input arguments:
%   t ....... time vector of solution returned by solver
%   state ... solution returned by solver (Matlab/Simulink)
%   fh ...... function handle to be evaluated
%   _________________
%   Output arguments:
%   ret.u ... values of control plant input (Stellgröße u = delta) [rad]
%   ret.s ... covered vehicle distance [m]
%   ________________________________________________
%   Structure of state
%   (x(1) ... Weg in Fahrzeugquerrichtung (sy_V))
%   (x(2) ... Fahrzeugquergeschwindigkeit (vy_V))
%   (x(3) ... Gierwinkel (psi))
%   (x(4) ... Giergeschwindigkeit (psiDot))
%   x(5) ... glob. Position in x_g-Richtung
%   x(6) ... glob. Position in y_g-Richtung
%   (x(7) ... glob. zurueckgelegte Wegstrecke)
%   ....     ... zusätzliche Zustände für function handle
% 
% Subject: lka
% Author: georgnoname
% Date: 19.10.2012 - 23.04.2013


% check input arguments
if nargin < 3; error('Not enough input arguments'); end%if
if nargin > 3; error('Too many input arguments'); end%if

% get number of time points
length_t = length(t);

% pre-allocation
ui(length_t,1) = 0;

% pre-calculation: Wegänderungen in x- und y-Richtung
diffx = diff(x(5,:));
diffy = diff(x(6,:));

% Länge der einzelnen Segmente
sSegment = sqrt(diffx.^2 + diffy.^2);

% values of control input
for i = 1:length_t
    
    % Stellgröße (anonymus function) (yL/epsL/kapL werden in Fkt.
    % einspurMdl_u berechnet, diese wird über function handle aufgerufen)
    ui(i) = fh(t(i),x(:,i));
            
end%for

% output arguments
ret.u.about = 'values of control plant input';
ret.u.value = ui;
ret.u.unit = 'rad';

ret.s.about = 'covered vehicle distance';
ret.s.value = [0;cumsum(sSegment)'];
ret.s.unit = 'm';

end%fcn
