clc
clear all

initFile_

addpath(genpath('dataStructure'));



%%

%%% vehicle model
p = paramFile2Struct('paramFile_SingleTrackMdl_BMW5');
SubMdl_Vehicle = VhclMdl_SingleTrack(...
	VhclMdl_SingleTrack_parameters(p),...
	VhclMdl_SingleTrack_initValues(0,1,2,3),...
	VhclMdl_SingleTrack_stateNames());

%%% steering model
p = paramFile2Struct('paramFile_SteeringMdl_CarMaker_DSR_');
SubMdl_Steering = SteerMdl_CarMaker_DSR(...
	SteerMdl_CarMaker_DSR_parameters(p),...
	SteerMdl_CarMaker_DSR_initValues(0,1),...
	SteerMdl_CarMaker_DSR_stateNames(),...
	'some info');

v1 = Vehicle(SubMdl_Vehicle,SubMdl_Steering);


