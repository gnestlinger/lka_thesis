function [traj,err] = lka_trajectory_03(swPlot,delta)
% Soll-Trajektorie: S-Kurve aus Klothoiden
% 
% Subject: lka
% Author: georgnoname
% Date: 19.10.2012 - 11.01.2013


% check input arguments
if nargin < 2; delta = []; end
if nargin < 1; swPlot = 0; end

% Kurzbezeichnung der Trajektorie
ID.label = 'S\_Cloth';
ID.nbr = 3;


% segment 1: straight
p1 = createTraj('straight',[0;0],[50;0],delta);

% segment 2: clothoid
kStop = 0.005;
A = 400;
% createTraj('cloth2',xyStart,kStart,kStop,A,slopeStart,delta)
p2 = createTraj('cloth2',[p1.x(end),p1.y(end)],0,kStop,A,0,delta);

% segment 3: clothoid
p3 = createTraj('cloth2',[p2.x(end),p2.y(end)],p2.k(end),0,A,p2.phi(end),delta);

% segment 4: straight
p4 = createTraj('straight',[p3.x(end),p3.y(end)],[p3.x(end)+50,p3.y(end)],delta);


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
%     title(['Sollbahn Nr. 03: ',label]);
% end%if

plotTraj(swPlot,traj)

end%fcn
