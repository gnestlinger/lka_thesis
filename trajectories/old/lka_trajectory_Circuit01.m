function [traj,err] = lka_trajectory_Circuit01(swPlot,delta)
% Soll-Trajektorie: 
% 
% Subject: lka
% Author: georgnoname
% Date: 16.05.2013 - 19.05.2013


% check input arguments
if nargin < 2; delta = []; end
if nargin < 1; swPlot = 0; end

% Kurzbezeichnung der Trajektorie
ID.label = 'Circuit01';
ID.nbr = 01;


% segment 1: 200m straight road
p1 = createTraj('straight',[0,0],[200;0],delta);
% p1.s(end)


% segment 2: curved road with curvature of 1/200 1/m
p2 = createTraj('curv',[p1.x(end),p1.y(end)],3/2*pi,4/2*pi,200,delta);
% p2.s(end)


% segment 3: clothoid
xyStart = [p2.x(end),p2.y(end)];
kStart = 1/200;
kStop = 1/100;
A = sqrt(104.72^2/(2*(2*pi/360*45)))+61.155;
slopeStart = pi/2;
% createTraj('cloth2',xyStart,kStart,kStop,A,slopeStart,delta)
p3 = createTraj('cloth2',xyStart,kStart,kStop,A,slopeStart,delta);
% p3.s(end)


% segment 4: clothoid
xyStart = [p3.x(end),p3.y(end)];
kStart = 1/100;
kStop = 0; % 1/1e6
A = -sqrt(157.06^2/(2*(2*pi/360*45)));
slopeStart = 3/4*pi;
% createTraj('cloth2',xyStart,kStart,kStop,A,slopeStart,delta)
p4 = createTraj('cloth2',xyStart,kStart,kStop,A,slopeStart,delta);
% p4.s(end)


% segment 5: clothoid
xyStart = [p4.x(end),p4.y(end)];
kStart = 0; % 1/1e6
kStop = 1/300;
A = -sqrt(3140.6^2/(2*(2*pi/360*300)));
slopeStart = pi;
% createTraj('cloth2',xyStart,kStart,kStop,A,slopeStart,delta)
p5 = createTraj('cloth2',xyStart,kStart,kStop,A,slopeStart,delta);
% p5.s(end)


% segment 6: clothoid
xyStart = [p5.x(end),p5.y(end)];
kStart = 1/300;
kStop = 0; % 1/1e6
A = sqrt(209.4^2/(2*(2*pi/360*20)));
slopeStart = p5.phi(end);
% createTraj('cloth2',xyStart,kStart,kStop,A,slopeStart,delta)
p6 = createTraj('cloth2',xyStart,kStart,kStop,A,slopeStart,delta);
% p6.s(end)


% segment 7: 500m straight road
xend = p6.x(end) + 500*cos(2*pi/360*(180-320));
yend = p6.y(end) + 500*sin(2*pi/360*(180-320));
p7 = createTraj('straight',[p6.x(end),p6.y(end)],[xend;yend],delta);
% p7.s(end)


% segment 8: clothoid
xyStart = [p7.x(end),p7.y(end)];
kStart = 0; % 1e6
kStop = 1/162.2;
A = sqrt(792.5^2/(2*(2*pi/360*140)));
slopeStart = 2*pi/360*(-140);
% createTraj('cloth2',xyStart,kStart,kStop,A,slopeStart,delta)
p8 = createTraj('cloth2',xyStart,kStart,kStop,A,slopeStart,delta);
% p8.s(end)


% segment 9: 900m straight road
xend = p8.x(end) + 900*cos(0);
yend = p8.y(end) + 900*sin(0);
p9 = createTraj('straight',[p8.x(end),p8.y(end)],[xend;yend],delta);
% p9.s(end)


% connect segments
[traj,err] = createTraj('connect',p1,p2,p3,p4,p5,p6,p7,p8,p9);
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

plotTraj(swPlot,traj);


% export
filename = 'traj_Circuit01_segments';

% save as mat-file
% save(filename,'traj');

% % save to *.dat file
% addpath('C:\Users\Georg\Documents\MATLAB\LKA')
% saveFigure_txt(filename,...
%     200,...
%     'x1',p1.x,...
%     'y1',p1.y,...
%     's1',p1.s,...
%     'kappa1',p1.k,...
%     'x2',p2.x,...
%     'y2',p2.y,...
%     's2',p2.s,...
%     'kappa2',p2.k,...
%     'x3',p3.x,...
%     'y3',p3.y,...
%     's3',p3.s,...
%     'kappa3',p3.k,...
%     'x4',p4.x,...
%     'y4',p4.y,...
%     's4',p4.s,...
%     'kappa4',p4.k,...
%     'x5',p5.x,...
%     'y5',p5.y,...
%     's5',p5.s,...
%     'kappa5',p5.k,...
%     'x6',p6.x,...
%     'y6',p6.y,...
%     's6',p6.s,...
%     'kappa6',p6.k,...
%     'x7',p7.x,...
%     'y7',p7.y,...
%     's7',p7.s,...
%     'kappa7',p7.k,...
%     'x8',p8.x,...
%     'y8',p8.y,...
%     's8',p8.s,...
%     'kappa8',p8.k,...
%     'x9',p9.x,...
%     'y9',p9.y,...
%     's9',p9.s,...
%     'kappa9',p9.k);

end%fcn

