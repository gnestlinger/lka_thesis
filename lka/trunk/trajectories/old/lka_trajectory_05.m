function [traj,err] = lka_trajectory_05(swPlot,delta)
% Soll-Trajektorie: Spurwechsel mit Kreisen
% 
% Subject: lka
% Author: georgnoname
% Date: 30.10.2012 - 11.01.2013


% check input arguments
if nargin < 2; delta = []; end
if nargin < 1; swPlot = 0; end

% Kurzbezeichnung der Trajektorie
ID.label = 'LaneChange\_Circ';
ID.nbr = 5;


% segment 1: straight
p1 = createTraj('straight',[0;0],[100;0],delta);

% segment 2: curvature
angStart1 = -1/2*pi;
angStop1 = -1/3*pi;
r2 = 1/0.002;
p2 = createTraj('curv',[p1.x(end),p1.y(end)],angStart1,angStop1,r2,delta);

% segment 3: curvature
angStart = pi/2 + pi/2-abs(angStop1);%pi/2+pi/6
angStop = -angStart1;
r3 = r2;
p3 = createTraj('curv',[p2.x(end),p2.y(end)],angStart,angStop,r3,delta);

% segment 4: straight
p4 = createTraj('straight',[p3.x(end),p3.y(end)],[p3.x(end)+100,p3.y(end)],delta);


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
%     title(['Sollbahn Nr. 05: ',label]);
% end%if

plotTraj(swPlot,traj);

end%fcn
