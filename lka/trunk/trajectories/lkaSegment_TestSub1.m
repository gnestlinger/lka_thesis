function lkaSegment_TestSub1(t,p,fig,varargin)

% Subject: lka
% $Author$
% $LastChangedDate$
% $Revision$

if nargin < 4
    varargin = {'bo','MarkerSize',4};
end%if

sbplt = @(nbr) subplot(3,2,nbr);

figure(fig);


% ax(1) = sbplt(1);
% plotTraj(1,t);
% 
% ax(2) = sbplt(2);
% plot(p);
% 
% linkaxes(ax);


plotStyle_old = {'k-','LineWidth',1};

% x/y
sbplt([1 3 5]);
plot(t.x,t.y,plotStyle_old{:});
hold all
plot(p,varargin{:});
legend('old','oop');


% curve length
ax(1) = sbplt([2]);
plot((0:length(t.s)-1)/length(t.s),t.s,plotStyle_old{:});
hold all
plotdiff_(p.segmentData,[],'s');
hold off
title('')

% curvature
ax(2) = sbplt([4]);
plot((0:length(t.k)-1)/length(t.k),t.k,plotStyle_old{:});
hold all
plotdiff_(p.segmentData,[],'k');
hold off
title('')

% tangent angle
ax(3) = sbplt([6]);
plotdiff_(p.segmentData,[],'phi');
title('')

linkaxes(ax,'x');
end%fcn

