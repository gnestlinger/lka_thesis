function [traj,err] = lka_trajectory_07(swPlot,delta)
% Soll-Trajektorie: sprungförmige Kurvenkrümmung
% 
% Subject: lka
% Author: georgnoname
% Date: 17.12.2012 - 24.04.2013


% check input arguments
if nargin < 2; delta = []; end
if nargin < 1; swPlot = 0; end

% Kurzbezeichnung der Trajektorie
ID.label = 'Circle\_R200';
ID.nbr = 7;


% segment 1: 100m straight road
p1 = createTraj('straight',[0,0],[50;0],delta);

% segment 2: curved road with curvature of 0.002 1/m
p2 = createTraj('curv',[p1.x(end),p1.y(end)],3/2*pi,4/2*pi,200,delta);


% connect segments
[traj,err] = createTraj('connect',p1,p2);
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
%     title(['Sollbahn Nr. 07: ',label]);
% end%if

plotTraj(swPlot,traj);

end%fcn
