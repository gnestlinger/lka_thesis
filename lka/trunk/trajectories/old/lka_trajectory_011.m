function [traj,err] = lka_trajectory_011(swPlot,delta)
% Soll-Trajektorie: sprungförmige Kurvenkrümmung
% 
% Daten aus: Vision-Based Lateral Control of Vehicles: Look-ahead and Delay
% Issues, Fig 10b.
% 
% Subject: lka
% Author: georgnoname
% Date: 24.09.2012 - 11.01.2013


% check input arguments
if nargin < 2; delta = []; end
if nargin < 1; swPlot = 0; end

% Kurzbezeichnung der Trajektorie
ID.label = 'Circle\_R500';
ID.nbr = 1;


% segment 1: 100m straight road
p1 = createTraj('straight',[0,0],[100;0],delta);

% segment 2: curved road with curvature of 0.002 1/m
p2 = createTraj('curv',[p1.x(end),p1.y(end)],3/2*pi,5/2*pi,1/0.002,delta);


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
%     title(['Sollbahn Nr. 01: ',label]);
% end%if

plotTraj(swPlot,traj)

end%fcn

