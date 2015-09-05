function [traj,err] = lka_trajectory_02(swPlot,delta)
% Soll-Trajektorie: geschlossene Trajektorie mit konst. Kurvenkrümmung
% 
% Subject: lka
% Author: georgnoname
% Date: 05.10.2012 - 11.01.2013


% check input arguments
if nargin < 2; delta = []; end
if nargin < 1; swPlot = 0; end

% Kurzbezeichnung der Trajektorie
ID.label = 'Oval\_Circle\_R500';
ID.nbr = 2;


% segment 1: 200m straight road
p1 = createTraj('straight',[0,0],[200;0],delta);

% segment 2: curved road with curvature of 0.002 1/m
r2 = 1/0.002;
p2 = createTraj('curv',[p1.x(end),p1.y(end)],3/2*pi,5/2*pi,r2,delta);

% segment 3: 100m straight road
p3 = createTraj('straight',[p2.x(end),p2.y(end)],[0;1000],delta);

% segment 4: curved road with curvature of 0.002 1/m
r4 = r2;
p4 = createTraj('curv',[p3.x(end),p3.y(end)],pi/2,3/2*pi,r4,delta);


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
%     title(['Sollbahn Nr. 02: ',label]);
% end%if

plotTraj(swPlot,traj)

end%fcn
