% Testfile: lkaSegmentClothoid
% 
% Tests for proper results with the four main settings (curvStart >(<)
% curvStop and A >(<) 0).

clear all;
close all
clc;

sbplt = @(nbr) subplot(2,2,nbr);


%% the clothoid itself

figure(1)

delta = 0.5;
curvStart = 0;
curvStop = 0.04;
slopeStart = 0;
A = 100;

a = lkaSegmentClothoid(delta,+curvStart,+curvStop,slopeStart,A);
b = lkaSegmentClothoid(delta,-curvStart,-curvStop,slopeStart,A);
c = lkaSegmentClothoid(delta,+curvStop,+curvStart,slopeStart,A);
d = lkaSegmentClothoid(delta,-curvStop,-curvStart,slopeStart,A);

ind = [1,101,201,401,701];
col = {'r','g'};

sbplt(1);
plottangent(a,ind,col{:});
title('a: should turn counter-clockwise');

sbplt(3); 
plottangent(b,ind,col{:});
title('b: should turn clockwise');

sbplt(2);
plottangent(c,ind,col{:});
title('c: should turn counter-clockwise');

sbplt(4);
plottangent(d,ind,col{:});
title('d: should turn clockwise');


%% the curve length

lkaSegmentClothoid_TestSub1(a,b,c,d,'s',2,'Curve length s');


%% the curvature

lkaSegmentClothoid_TestSub1(a,b,c,d,'k',3,'Curvature k');

%% the tangent angle

lkaSegmentClothoid_TestSub1(a,b,c,d,'phi',4,'Tangent angle phi');

