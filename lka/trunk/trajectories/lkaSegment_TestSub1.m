function lkaSegment_TestSub1(t,p,fig,varargin)

if nargin < 4
    varargin = {'ro','MarkerSize',4,'MarkerFaceColor','r'};
end%if

sbplt = @(nbr) subplot(2,2,nbr);

figure(fig);


ax(1) = sbplt(1);
plotTraj(1,t);

ax(2) = sbplt(2);
plot(p);

linkaxes(ax);


sbplt([3 4]);
plot(t.x,t.y,'bo','MarkerSize',8,'MarkerFaceColor','b');
hold all
plot(p,varargin{:});
legend('old','oop');

end%fcn

