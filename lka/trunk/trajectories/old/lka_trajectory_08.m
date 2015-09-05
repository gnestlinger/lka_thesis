function [traj,err] = lka_trajectory_08(swPlot,delta)
% Soll-Trajektorie: rampenförmige Kurvenkrümmung
% 
% Subject: lka
% Author: georgnoname
% Date: 17.12.2012 - 24.04.2013


% check input arguments
if nargin < 2; delta = []; end
if nargin < 1; swPlot = 0; end

% Kurzbezeichnung der Trajektorie
ID.label = 'Cloth\_K0p005';
ID.nbr = 8;


% segment 1: straight
p1 = createTraj('straight',[0;0],[50;0],delta);

% segment 2: clothoid
xyStart = [p1.x(end),p1.y(end)];
kStart = 0;
kStop = 0.005;
A = 400;
slopeStart = 0;
% createTraj('cloth2',xyStart,kStart,kStop,A,slopeStart,delta)
p2 = createTraj('cloth2',xyStart,kStart,kStop,A,slopeStart,delta);

% zur Kontrolle
% p3 = createTraj('curv',[p2.x(end),p2.y(end)],p2.phi(end)-pi/2,5/2*pi,200,delta);

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
%     title(['Sollbahn Nr. 08: ',label]);
% end%if

plotTraj(swPlot,traj);

end%fcn
