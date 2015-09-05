function [] = plotTraj(swPlot,trajIn)
% plotTraj   Visualisierung von Sollbahnen
%   _______
%   Syntax:
%   [] = plotTraj(swPlot,trajIn,ID)
%   ________________
%   Input arguments:
%   swPlot ... Switch ob plotten (swPlot=1) oder nicht (sonst) 
%   trajIn ... Sollbahn-Structure basierend auf createTraj.m
% 
% Subject: lka
% Author: georgnoname
% Date: 11.01.2013


% plot or not
if swPlot ~= 1; return; end%if

% create vetor of interesting indices to separate segments of different
% type (straight, curve,...)
ind = 0;
ind = [ind,find(diff(trajIn.type))'];
ind = [ind,length(trajIn.x)];

% detect steps in sign of curvature
vorz = sign(trajIn.k);
ind = [ind,find(diff(vorz) < -1)];
ind = [ind,find(diff(vorz) > +1)];

% sort
ind = sort(ind);

% preallocation of X and Y
X(1:max(diff(ind)),1:length(ind)-1) = NaN;
Y(1:max(diff(ind)),1:length(ind)-1) = NaN;

% store distinct segments in columns of X/Y
for i = 1:length(ind)-1
    x = col(trajIn.x(ind(i)+1:ind(i+1)));
    X(1:length(x),i) = x;
    
    y = col(trajIn.y(ind(i)+1:ind(i+1)));
    Y(1:length(y),i) = y;
end%for


% plot
% figure;

% odd columns
plot(X(:,1:2:end),Y(:,1:2:end),'o',...
        'MarkerSize',4);

% enable distinct colours
hold all;

% even columns
plot(X(:,2:2:end),Y(:,2:2:end),'.',...
        'MarkerSize',5);
    
% release plot
hold off

% labels, ...
grid on; 
axis equal;
xlabel('x_0 [m]'); ylabel('y_0 [m]');
title(['Sollbahn Nr. ',num2str(trajIn.ID.nbr),': ',trajIn.ID.label]);

end%fcn



function ret = col(in)
% COL column vector
% check if 'in' is of dimension 1*x or x*1 (vector, not matrix)

[a,b] = size(in);

if (a>1) && (b>1)
    error('Input Argument should be vector but is matrix.')
end

ret = in(:);

end%fcn
