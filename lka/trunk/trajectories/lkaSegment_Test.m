% Testfile: lkaSegment.m
% 
% Compares the results of the old script with the new oop approach.

clear all;
close all;
clc;

sbplt = @(nbr) subplot(2,1,nbr);

addpath('old')

%% trajectory_01

t = lka_trajectory_01;

a = lkaSegmentStraight([],100,0);
b = lkaSegmentCircle([],-pi/2,0,500);
p = a + b;

lkaSegment_TestSub1(t,p,1);


%% trajectory_02
t = lka_trajectory_02;

a = lkaSegmentStraight([],200,0);
b = lkaSegmentCircle([],-pi/2,pi/2,500);
c = lkaSegmentStraight([],200,pi);
d = lkaSegmentCircle([],pi/2,3*pi/2,500);
p = a + b + c + d;

lkaSegment_TestSub1(t,p,2);


%% trajectory_03
t = lka_trajectory_03;

a = lkaSegmentStraight([],50,0);
b = lkaSegmentClothoid([],0,0.005,0,400);
c = lkaSegmentClothoid([],-b.segmentData.k(end),0,b.segmentData.phi(end),400);
d = lkaSegmentStraight([],50,0);
p = a + b + c + d;

lkaSegment_TestSub1(t,p,3);


%% trajectory_04
t = lka_trajectory_04;

a = lkaSegmentStraight([],100,0);
b = lkaSegmentClothoid([],0,0.004,0,300);
c = lkaSegmentClothoid([],-b.segmentData.k(end),0,b.segmentData.phi(end),300);
d = lkaSegmentStraight([],100,0);
p = a + b + c + d;

lkaSegment_TestSub1(t,p,4);


%% trajectory_05
t = lka_trajectory_05;

a = lkaSegmentStraight([],100,0);
b = lkaSegmentCircle([],-pi/2,-pi/3,500);
c = lkaSegmentCircle([],pi-pi/3,pi/2,500);
d = lkaSegmentStraight([],100,0);
p = a + b + c + d;

lkaSegment_TestSub1(t,p,5);


%% trajectory_06
t = lka_trajectory_06;

a = lkaSegmentStraight([],100,0);
b = lkaSegmentClothoid([],0,0.005,0,400);
c = lkaSegmentClothoid([],b.segmentData.k(end),0,b.segmentData.phi(end),400);
d = lkaSegmentStraight([],100,c.segmentData.phi(end));
p = a + b + c + d;

lkaSegment_TestSub1(t,p,6);


%% trajectory_07
t = lka_trajectory_07;

a = lkaSegmentStraight([],50,0);
b = lkaSegmentCircle([],-pi/2,0,200);
p = a + b;

lkaSegment_TestSub1(t,p,7);


%% trajectory_08

tic
t = lka_trajectory_08;
toc

tic
a = lkaSegmentStraight([],50,0);
b = lkaSegmentClothoid([],0,0.005,0,400);
p = a + b;
toc

lkaSegment_TestSub1(t,p,8);


%% trajectory_Circuit01

dispString = '%s %.2f %s\n\n';

tic
t = lka_trajectory_Circuit01;
time = toc;
fprintf(dispString,'Funktion:',time,'sec');


a = lkaSegmentStraight([],200,0);
b = lkaSegmentCircle([],3/2*pi,4/2*pi,200);

A = sqrt(104.72^2/(2*(2*pi/360*45)))+61.155;
c = lkaSegmentClothoid([],1/200,1/100,pi/2,A);

A = sqrt(157.06^2/(2*(2*pi/360*45)));
d = lkaSegmentClothoid([],1/100,0,3/4*pi,A);

A = sqrt(3140.6^2/(2*(2*pi/360*300)));
e = lkaSegmentClothoid([],0,-1/300,pi,A);

A = sqrt(209.4^2/(2*(2*pi/360*20)));
f = lkaSegmentClothoid([],-1/300,0,-2.0911,A);

g = lkaSegmentStraight([],500,2*pi/360*(180-320));

A = sqrt(792.5^2/(2*(2*pi/360*140)));
h = lkaSegmentClothoid([],0,1/162.2,2*pi/360*(-140),A);

i = lkaSegmentStraight([],900,0);



tic
p = a + b + c + d + e + f + g + h + i;
time = toc;
fprintf(dispString,'OOP:',time,'sec');

lkaSegment_TestSub1(t,p,09);

