



vx = 10;
LAD = 0;
% vx = realp('vx',vx);
% LAD = realp('LAD',LAD);
vx = [];
LAD = [];

paramsVhcl = paramFile2Struct('paramFile_SingleTrackMdl_CarMaker_BMW5');

sys_a = ssMdl_SingleTrack('st',paramsVhcl,vx,LAD);

sys_b = ssMdl_SingleTrack('stvis',paramsVhcl,vx,LAD);

sys_c = ssMdl_SingleTrack('stvis_2int',paramsVhcl,vx,LAD,...
	'paramFile_SteeringMdl_CarMaker_DSR');

sys_d = ssMdl_SingleTrack('stdsr',paramsVhcl,vx,LAD,...
	'paramFile_SteeringMdl_CarMaker_DSR');

sys_e = ssMdl_SingleTrack('stvisdsr',paramsVhcl,vx,LAD,...
	'paramFile_SteeringMdl_CarMaker_DSR');

sys_f = ssMdl_SingleTrack('stvisdsr_2int',paramsVhcl,vx,LAD,...
	'paramFile_SteeringMdl_CarMaker_DSR');
