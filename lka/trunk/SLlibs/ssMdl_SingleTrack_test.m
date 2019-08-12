

clc

vx = 10;
LAD = 0;
% vx = realp('vx',vx);
% LAD = realp('LAD',LAD);
vx = [];
LAD = [];

paramsVhcl = paramFile2Struct('paramFile_SingleTrackMdl_CarMaker_BMW5');

sys_a = ssMdl_SingleTrack('st',paramsVhcl,vx);

sys_b = getMergedSSMdl('stvis',paramsVhcl,vx,LAD);

sys_c = getMergedSSMdl('stvis_2int',paramsVhcl,vx,LAD,...
	'paramFile_SteeringMdl_CarMaker_DSR');

sys_d = getMergedSSMdl('stdsr',paramsVhcl,vx,LAD,...
	'paramFile_SteeringMdl_CarMaker_DSR');

sys_e = getMergedSSMdl('stvisdsr',paramsVhcl,vx,LAD,...
	'paramFile_SteeringMdl_CarMaker_DSR');

sys_f = getMergedSSMdl('stvisdsr_2int',paramsVhcl,vx,LAD,...
	'paramFile_SteeringMdl_CarMaker_DSR');
