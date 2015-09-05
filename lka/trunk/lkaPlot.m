function [p1,p2] = lkaPlot(solIn,lad,chPlot)
% lkaPlot   plots solutions stored in 'solIn.*' of lka-simulation   
%   _______
%   Syntax:
%   [p1,p2] = lkaPlot(solIn,lad,chPlot)
%   ________________
%   Input arguments:
%   solIn ... Structure of solution solIn.<lbl1> = ode**() has to be as 
%             follows (<lbl1> is used as label in plot):
%             solIn.<lbl1>.simIn 
%               vx ..... longitudinal velocity [m/s]
%               lad .... look-ahead distance [m]
%               traj ... 
%             solIn.<lbl1>.simOut 
%               t .............. time-vector
%               vehicleState ... solution of vehicle-DEQ-system
%               lkaSensor ...... yL/epsL/kapL computed by visionSystem.m
%               simParameter ... Parameter used in simulation like vx, lad
%               u .............. anonymus funtion of control input
%             solIn.<lbl1>.simPost
%               
%   lad ..... look-ahead distance für die yL/epsL/kapL geplotet werden
%             sollen
%   chPlot .. wähle Plot: Trajektorie ('traj'), yL/u ('add') oder beide
%   _________________
%   Output arguments:
%   p1 ... handle to figure of trajectories
%   p2 ... handle to figure of yL/u
%   
% 
% Subject: lka
% Author: georgnoname
% Date: xx.10.2012 - 25.02.2013


% check input arguments
if nargin < 3; chPlot = 'undefined'; end%if
if nargin < 2; lad = []; end%if


% get field names of struct 'solIn'
fNamesSolIn = fieldnames(solIn);

% number of field names
N = length(fNamesSolIn);

% create colormap
% colmap = lines(N);
% colmap = brighten(colmap,0.85);

% set(0,'DefaultAxesColorOrder',[1 0 0;0 1 0;0 1 1;1 0 1]);
set(0,'DefaultAxesLineStyleOrder','-|--|:');

% create labels
lblList{N} = [];
for i = 1:N
lblList{i} = [regexprep(fNamesSolIn{i},'_','\\_'),...
    ' (',...
    'v_x=',num2str(solIn.(fNamesSolIn{i}).simIn.vx),'m/s, ',...
    'L=',num2str(solIn.(fNamesSolIn{i}).simIn.lad),'m',...
    ')'];
end%for

% switch to plot trajectory or yL+u or both
switch lower(chPlot)
    case 'traj'
        p1 = plotTrajectory(solIn,fNamesSolIn,lblList);
        p2 = [];
        
    case 'add'
        p1 = [];
        p2 = plotAddData(solIn,fNamesSolIn,lblList,lad);
        
    otherwise
        [p1,~] = lkaPlot(solIn,lad,'traj');
        [~,p2] = lkaPlot(solIn,lad,'add');

end%switch

end%fcn



function p = plotTrajectory(solIn,fNamesSolIn,lblList)

% init figures
fig = figure;
ax = gca;

% get length of longest t-vectors of solutions stored in solIn.*
nmax = maxLength(solIn,fNamesSolIn);

% preallocation of X and Y
X(1:nmax,1:length(fNamesSolIn)) = NaN;
Y(1:nmax,1:length(fNamesSolIn)) = NaN;

for i = 1:length(fNamesSolIn)
    
    % choose i-th solution corresponding to control-strategy solIn.i
    soli = solIn.(fNamesSolIn{i});
    
    % store x0 and y0-values in matrices X and Y
    rowEnd = length(soli.simOut.vehicleState.x_g);
    X(1:rowEnd,i) = soli.simOut.vehicleState.x_g;
    Y(1:rowEnd,i) = soli.simOut.vehicleState.y_g;  
    
end%for

% plot intended trajectory
hold(ax,'all');
p1 = plot(ax,...
    soli.simIn.traj.x,soli.simIn.traj.y,'-k');

% plot simulation results
p2 = plot(ax,X,Y); 
set(p1,'LineWidth',4);
set(p2,'LineWidth',1.5)

% set current axes
% set(0,'CurrentFigure',fig); 
% set(gcf,'CurrentAxes',ax);

% set labels
legend('Sollbahn',lblList{:},'Location','Best');
axis equal;
ylabel('y_{g} [m]'); xlabel('x_{g} [m]');
% title(['v_x = ',num2str(sol.simParameter.vx),' m/s','   ',...
%     'L = ',num2str(sol.simParameter.lad),' m']);
grid on;

p = [p1;p2];

end%fcn



function p = plotAddData(solIn,fNamesSolIn,lblList,lad)

% init figures
fig = figure;
ax1 = subplot(2,1,1);
ax2 = subplot(2,1,2);

% flag to decide: plot over x_g or s_g
flag_xs = find(diff(solIn.(fNamesSolIn{1}).simIn.traj.x) <= 0, 1);

% get length of longest t-vectors of solutions stored in solIn.*
nmax = maxLength(solIn,fNamesSolIn);

% preallocation of data to plot
X(1:nmax,1:length(fNamesSolIn)) = NaN;
YL(1:nmax,1:length(fNamesSolIn)) = NaN;
U(1:nmax,1:length(fNamesSolIn)) = NaN;

% für alle simulierten Regler ('fNamesSolIn') speichere Daten in Matrizen
for i = 1:length(fNamesSolIn)
    
    % choose i-th solution corresponding to control-strategy solIn.i
    soli = solIn.(fNamesSolIn{i});
    
    % get lad if undefined
    if isempty(lad); lad = soli.simIn.lad; end%if
    
    % überprüfe ob yL usw. für gewünschtes lad bereits vorliegen
    fNameslad = fieldnames(soli.simPost.relPos);
    ind = [];
    for j = 1:length(fNameslad)
        try
            ladExist = soli.simPost.relPos.(fNameslad{j}).lad.value;
            if ladExist == lad; ind = j; break; end%if
        catch expression
%             disp(expression.message);
        end%try
    end%for
    
    % if cond. is true: choose alredy computed addData from struct soli
    if ~isempty(ind)
        yL = soli.simPost.relPos.(fNameslad{ind}).yL.value;
        epsL = soli.simPost.relPos.(fNameslad{ind}).epsL.value;
    else
        ret = calc_ladDependentData(soli,lad);
        fNameret = fieldnames(ret);
        yL = ret.(fNameret{1}).yL.value;
        epsL = ret.(fNameret{1}).epsL.value;     
    end%if
       
    % store in matrix
    % if x_g is strictly increasing, plot over x, else plot over s
    if isempty(flag_xs)
        X(1:length(soli.simOut.vehicleState.x_g),i) = soli.simOut.vehicleState.x_g;
        xLbl = 'x_g [m]';
    else
        X(1:length(soli.simOut.vehicleState.s_g),i) = soli.simOut.vehicleState.s_g;
        xLbl = 's_g [m]';
    end%if
    
    
    % mit Matlab simuliert (u kann erst nach Simulation berechnet werden)
    switch soli.simProg
        case {'Matlab'}
            u = rad2deg(soli.simPost.controlInp.steeringAngle.value);
        case {'Simulink'}
            u = rad2deg(soli.simOut.controlInp.steeringAngle);
    end%switch
    
    % store in matrix
    YL(1:length(soli.simOut.vehicleState.x_g),i) = yL;
    U(1:length(soli.simOut.vehicleState.x_g),i) = u;
         
end%for

% plot yL
p(:,1) = plot(ax1,X,YL,'LineWidth',1.5);

% plot Stellgröße
p(:,2) = plot(ax2,X,U,'LineWidth',1.5);


% subfigure 1: set current axes
set(0,'CurrentFigure',fig);
set(gcf,'CurrentAxes',ax1);

% subfigure 1: set labels
title(['Querversatz y_L bei L = ',num2str(lad),'m']);
ylabel('y_{L} [m]'); %xlabel('x_{g} [m]');
% legend(['L = ',num2str(addData.yL.lad),addData.yL.unit]);
legend(lblList{:},'Location','Best');
grid on
box on

% subfigure 2: set current axes
set(gcf,'CurrentAxes',ax2);

% subfigure 2: set labels
title('Stellgröße u = \delta_f'),
ylabel('u [°]'); xlabel(xLbl);
grid on
box on

end%fcn



function nmax = maxLength(solIn,fNamesSolIn)
% get length of longest t-vectors of solutions stored in solIn.*

nmax = 0;
for i = 1:length(fNamesSolIn)
    n = length(solIn.(fNamesSolIn{i}).simOut.t);
    if n > nmax; nmax = n; end
end%for

end%fcn
