function [traj,err] = lka_trajectory_06(swPlot,delta)
% Soll-Trajektorie: 'gewagte Autobahnabfahrt'
% 
% Subject: lka
% Author: georgnoname
% Date: 09.11.2012 - 11.01.2013


% check input arguments
if nargin < 2; delta = []; end
if nargin < 1; swPlot = 0; end

% Kurzbezeichnung der Trajektorie
ID.label = '';
ID.nbr = 6;


% segment 1: straight
p1 = createTraj('straight',[0;0],[100;0],delta);

% segment 2: clothoid
xyStart = [p1.x(end),p1.y(end)];
kStart = 0;
kStop = 0.005;
A = 400;
slopeStart = 0;
% createTraj('cloth2',xyStart,kStart,kStop,A,slopeStart,delta)
p2 = createTraj('cloth2',xyStart,kStart,kStop,A,slopeStart,delta);

% segment 3: clothoid
xyStart = [p2.x(end),p2.y(end)];
kStart = p2.k(end);
kStop = 0;
% A = -A;
slopeStart = p2.phi(end);
% createTraj('cloth2',xyStart,kStart,kStop,A,slopeStart,delta)
p3 = createTraj('cloth2',xyStart,kStart,kStop,-A,slopeStart,delta);

% segment 4: straight
[a,b] = pol2cart(p3.phi(end),100);
p4end = [p3.x(end),p3.y(end)] + [a,b];
p4 = createTraj('straight',[p3.x(end),p3.y(end)],p4end,delta);

% p5end = [p3.x(end),p3.y(end)] - [a,b];
% p5 = createTraj('straight',[p3.x(end),p3.y(end)],p5end,delta);


% connect segments
[traj,err] = createTraj('connect',p1,p2,p3,p4);
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
%     title(['Sollbahn Nr. 06: ',label]);
% end%if

plotTraj(swPlot,traj);

end%fcn
