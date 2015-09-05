function [traj,err] = lka_trajectory_straight(swPlot,delta)
% Soll-Trajektorie: straight segment
% 
% Subject: lka
% Author: georgnoname
% Date: 19.10.2012 - 28.04.2013


% check input arguments
if nargin < 2; delta = []; end
if nargin < 1; swPlot = 0; end

% Kurzbezeichnung der Trajektorie
ID.label = 'straight';
ID.nbr = [];


% segment 1: 200m straight road
p1 = createTraj('straight',[0,0],[500;0],delta);


% connect segments
[traj,err] = createTraj('connect',p1);
traj.ID = ID;

% % plot
% if swPlot == 1
%     figure;
%     plot(traj.x,traj.y,'ok',...
%         'MarkerSize',2,...
%         'MarkerFaceColor','k'); 
%     grid on; 
%     axis equal;
%     xlabel('x_0 [m]'); ylabel('y_0 [m]');
%     title(['Sollbahn: ',label]);
% end%if
        
plotTraj(swPlot,traj);

end%fcn
