function ret = calc_relPos(solIn,lad)
% calc_relPos   Calculate relative position of vehicle to intended
% trajectory based on vehicle position/orientation.
%   _______
%   Syntax:
%   ret = calc_relPos(solIn,lad)
%   ________________
%   Input arguments:
%   solIn ... solution returned by solver (Matlab/Simulink) including
%             intended trajectory in solIn.simIn.traj
%   lad ..... look-ahead distance(s) für die yL/epsL/kapL berechnet werden
%             sollen
%   _________________
%   Output arguments:
%   ret.lad_<lad(i)>.lad .... lad-value <lad(i)> used to compute yL/... [m]
%   ret.lad_<lad(i)>.yL ..... lateral error yL [m]
%   ret.lad_<lad(i)>.epsL ... relative angle epsL [rad]
%   ret.lad_<lad(i)>.kapL ... cuvature kapL [1/m]
%   ________________________________________________
%   Structure of solIn (only required fields listed):
%   | simProg
%   | simIn
%       - traj
%   | simOut
%       - t
%       - vehicleState: sy, vy, psi, psiDot, x_g, y_g, (yL, epsL)
% 
% Subject: lka
% Author: georgnoname
% Date: 19.10.2012 - 01.05.2013


% check input arguments
if nargin < 2; lad = 0; end%if

% sort lad and keep only unique entries
lad = unique(lad);

% get number of time points
length_t = length(solIn.simOut.t);

% pre-allocation
yL(length_t,1) = 0;
epsL(length_t,1) = 0;
kapL(length_t,1) = 0;

% run loop for all values in lad
for j = 1:length(lad)
    
    % Berechne yL/epsL/kapL für look-ahead-distance
    for i = 1:length_t
        [~,yL(i),epsL(i),kapL(i)] = lkaLaneTracking(...
            [...
            solIn.simOut.vehicleState.x_g(i),...
            solIn.simOut.vehicleState.y_g(i),...
            solIn.simOut.vehicleState.psi(i)],...
            solIn.simIn.traj_sd,lad(j));
    end%for

    % output arguments
    ladStr = num2str(lad(j));
    ret.(['lad_',ladStr]).lad.about = 'look-ahead distance lad';
    ret.(['lad_',ladStr]).lad.value = lad(j);
    ret.(['lad_',ladStr]).lad.unit = 'm';

    ret.(['lad_',ladStr]).yL.about = ['yL at look-ahead distance ',ladStr,'m'];
    ret.(['lad_',ladStr]).yL.value = yL;
    ret.(['lad_',ladStr]).yL.unit = 'm';

    ret.(['lad_',ladStr]).epsL.about = ['epsL at look-ahead distance ',ladStr,'m'];
    ret.(['lad_',ladStr]).epsL.value = epsL;
    ret.(['lad_',ladStr]).epsL.unit = 'rad';
    
    ret.(['lad_',ladStr]).kapL.about = ['kapL at look-ahead distance ',ladStr,'m'];
    ret.(['lad_',ladStr]).kapL.value = kapL;
    ret.(['lad_',ladStr]).kapL.unit = '1/m';

end%for
 
end%fcn
