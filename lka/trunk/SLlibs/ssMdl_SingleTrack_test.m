

clc

vx = 10;
LAD = 0;
% vx = realp('vx',vx);
% LAD = realp('LAD',LAD);
% vx = [];
% LAD = [];

paramsVhcl = paramFile2Struct('paramFile_SingleTrackMdl_CarMaker_BMW5');
paramsSteer = paramFile2Struct('paramFile_SteeringMdl_CarMaker_DSR');

sys_a = ssMdl_singleTrack('stm', paramsVhcl, vx);

sys_b = getMergedSSMdl('stvis', paramsVhcl, vx, LAD);

sys_c = getMergedSSMdl('stvis_1int', paramsVhcl, vx, LAD, paramsSteer);

sys_d = getMergedSSMdl('stvis_2int', paramsVhcl, vx, LAD, paramsSteer);

sys_e = getMergedSSMdl('stdsr', paramsVhcl, vx, LAD, paramsSteer);

sys_f = getMergedSSMdl('stvisdsr', paramsVhcl, vx, LAD, paramsSteer);

sys_g = getMergedSSMdl('stvisdsr_1int', paramsVhcl, vx, LAD, paramsSteer);

sys_h = getMergedSSMdl('stvisdsr_2int', paramsVhcl, vx, LAD, paramsSteer);
