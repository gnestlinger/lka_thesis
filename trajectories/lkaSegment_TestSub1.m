function lkaSegment_TestSub1(t,p,fig,varargin)

% Subject: lka
% $Author$
% $LastChangedDate$
% $Revision$

if nargin < 4
    varargin = {'bo','MarkerSize',4};
end%if

sbplt = @(nbr) subplot(5,2,nbr);

figure(fig);


ax(1) = sbplt(1);
plotTraj(1,t);

ax(2) = sbplt(2);
plot(p);

linkaxes(ax);


plotStyle_old = {'r-','LineWidth',1};

% x/y
sbplt([3 4]);
plot(t.x,t.y,plotStyle_old{:});
hold all
plot(p,varargin{:});
legend('old','oop');


% curve length
sbplt([5 6]);
plot((1:length(t.s))/length(t.s),t.s,plotStyle_old{:});
hold all
plotdiff_(p.segmentData,[],'s');
title('')

% curvature
sbplt([7 8]);
plot(t.s,t.k,plotStyle_old{:});
hold all
plotdiff_(p.segmentData,[],'k');
title('')

% tangent angle
sbplt([9 10]);
plotdiff_(p.segmentData,[],'phi');
title('')

end%fcn

