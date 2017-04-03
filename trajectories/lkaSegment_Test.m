% Testfile: lkaSegment.m
% 
% Compares the results of the old script with the new OOP approach.
% 
% Subject: lka
% $Author$
% $LastChangedDate$
% $Revision$

clear;
close all;
clc;

sbplt = @(nbr) subplot(2,1,nbr);

addpath('old')

%% trajectory_01

t = lka_trajectory_01;
deltaSet = 10;
a = lkaSegmentStraight(deltaSet,100,0);
b = lkaSegmentCircle(deltaSet,-pi/2,0,500);
p = a + b;

lkaSegment_TestSub1(t,p,1);


%% trajectory_02
t = lka_trajectory_02;
deltaSet = 50;
a = lkaSegmentStraight(deltaSet,200,0);
b = lkaSegmentCircle(deltaSet,-pi/2,pi/2,500);
c = lkaSegmentStraight(deltaSet,200,pi);
d = lkaSegmentCircle(deltaSet,pi/2,3*pi/2,500);
p = a + b + c + d;

lkaSegment_TestSub1(t,p,2);


%% trajectory_03
t = lka_trajectory_03;
deltaSet = 20;
a = lkaSegmentStraight(deltaSet,50,0);
b = lkaSegmentClothoid(deltaSet,0,0.005,0,400);
c = lkaSegmentClothoid(deltaSet,-b.segmentData.k(end),0,b.segmentData.phi(end),400);
d = lkaSegmentStraight(deltaSet,50,0);
p = a + b + c + d;

lkaSegment_TestSub1(t,p,3);


%% trajectory_04
t = lka_trajectory_04;
deltaSet = 10;
a = lkaSegmentStraight(deltaSet,100,0);
b = lkaSegmentClothoid(deltaSet,0,0.004,0,300);
c = lkaSegmentClothoid(deltaSet,-b.segmentData.k(end),0,b.segmentData.phi(end),300);
d = lkaSegmentStraight(deltaSet,100,0);
p = a + b + c + d;

lkaSegment_TestSub1(t,p,4);


%% trajectory_05
t = lka_trajectory_05;
deltaSet = 10;
a = lkaSegmentStraight(deltaSet,100,0);
b = lkaSegmentCircle(deltaSet,-pi/2,-pi/3,500);
c = lkaSegmentCircle(deltaSet,pi-pi/3,pi/2,500);
d = lkaSegmentStraight(deltaSet,100,0);
p = a + b + c + d;

lkaSegment_TestSub1(t,p,5);


%% trajectory_06
t = lka_trajectory_06;
deltaSet = 20;
a = lkaSegmentStraight(deltaSet,100,0);
b = lkaSegmentClothoid(deltaSet,0,0.005,0,400);
c = lkaSegmentClothoid(deltaSet,b.segmentData.k(end),0,b.segmentData.phi(end),400);
d = lkaSegmentStraight(deltaSet,100,c.segmentData.phi(end));
p = a + b + c + d;

lkaSegment_TestSub1(t,p,6);


%% trajectory_07
t = lka_trajectory_07;
deltaSet = 5;
a = lkaSegmentStraight(deltaSet,50,0);
b = lkaSegmentCircle(deltaSet,-pi/2,0,200);
p = a + b;

lkaSegment_TestSub1(t,p,7);


%% trajectory_08

tic
t = lka_trajectory_08;
toc

tic
deltaSet = 20;
a = lkaSegmentStraight(deltaSet,50,0);
b = lkaSegmentClothoid(deltaSet,0,0.005,0,400);
p = a + b;
toc

lkaSegment_TestSub1(t,p,8);


%% trajectory_Circuit01

dispString = '%s %.2f %s\n\n';

tic
t = lka_trajectory_Circuit01;
time = toc;
fprintf(dispString,'Funktion:',time,'sec');

deltaSet = 50;

a = lkaSegmentStraight(deltaSet,200,0);
b = lkaSegmentCircle(deltaSet,-pi/2,0,200);

A = sqrt(104.72^2/(2*(2*pi/360*45)))+61.155;
c = lkaSegmentClothoid(deltaSet,1/200,1/100,pi/2,A);

A = sqrt(157.06^2/(2*(2*pi/360*45)));
d = lkaSegmentClothoid(deltaSet,1/100,0,3/4*pi,A);

A = sqrt(3140.6^2/(2*(2*pi/360*300)));
e = lkaSegmentClothoid(deltaSet,0,-1/300,pi,A);

A = sqrt(209.4^2/(2*(2*pi/360*20)));
f = lkaSegmentClothoid(deltaSet,-1/300,0,-2.0911,A);

g = lkaSegmentStraight(deltaSet,500,2*pi/360*(180-320));

A = sqrt(792.5^2/(2*(2*pi/360*140)));
h = lkaSegmentClothoid(deltaSet,0,1/162.2,2*pi/360*(-140),A);

i = lkaSegmentStraight(deltaSet,900,0);



tic
p = a + b + c + d + e + f + g + h + i;
time = toc;
fprintf(dispString,'OOP:',time,'sec');

lkaSegment_TestSub1(t,p,09);

