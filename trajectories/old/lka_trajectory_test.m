function [traj,err] = lka_trajectory_test(swPlot,delta)
% Soll-Trajektorie: Testfunktion
% 
% Subject: lka
% Author: georgnoname
% Date: 19.10.2012 - 11.01.2013


% check input arguments
if nargin < 2; delta = []; end
if nargin < 1; swPlot = 0; end

% Kurzbezeichnung der Trajektorie
ID.label = 'test';
ID.nbr = [];


% segment 1:
kStop = 0.003;
% createTraj('cloth2',xyStart,kStart,kStop,A,slopeStart,delta)
% p1 = createTraj('cloth2',[],0,kStop,335);
p1 = createTraj('curv',[0 0],pi,-1/2*pi,500,delta);

% segment 2:
% p2 = createTraj('curv',[p1.x(end),p1.y(end)],p1.phi(end),0,1/kStop,delta);


% connect segments
[traj,err] = createTraj('connect',p1);
traj.ID = ID;

% figure
% plot(p1.x,p1.y,'.',...
%     p2.x,p2.y,'-'); 
% grid on; axis equal;

% % plot
% if swPlot == 1
%     figure;
%     plot(traj.x,traj.y,'-',...
%         'MarkerSize',5,...
%         'MarkerFaceColor','blue');
%     hold on
%     % start point
%     plot(p1.x(1),p1.y(1),'.k','MarkerSize',15);
%     grid on; axis equal;
%     xlabel('x [m]'); ylabel('y [m]');
% end%if

plotTraj(swPlot,traj);

end%fcn
