function lkaSegmentClothoid_TestSub1(a,b,c,d,objprop,fig,ylbl,varargin)


figure(fig);

toplot = @(obj) plot(...
    1:length(obj.segmentData.(objprop)),...
    obj.segmentData.(objprop),varargin{:});
sbplt = @(nbr) subplot(2,2,nbr);

xlbl = 'Index';

sbplt(1);
toplot(a);
grid on;
title(['a: ',objprop]);
xlabel(xlbl);
ylabel(ylbl);

sbplt(3); 
toplot(b);
grid on;
title(['b: ',objprop]);
xlabel(xlbl);
ylabel(ylbl);

sbplt(2);
toplot(c);
grid on;
title(['c: ',objprop]);
xlabel(xlbl);
ylabel(ylbl);

sbplt(4);
toplot(d);
grid on;
title(['d: ',objprop]);
xlabel(xlbl);
ylabel(ylbl);


end%fcn

