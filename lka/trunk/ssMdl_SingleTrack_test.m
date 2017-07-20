



initFile_

vx = 10;
LAD = 0;
% vx = realp('vx',vx);
% LAD = realp('LAD',LAD);
vx = [];
LAD = [];

sys_a = ssMdl_SingleTrack('st','paramFile_SingleTrackMdl_BMW5',vx,LAD);

sys_b = ssMdl_SingleTrack('stvis','paramFile_SingleTrackMdl_BMW5',vx,LAD);

sys_c = ssMdl_SingleTrack('stvis_2int','paramFile_SingleTrackMdl_BMW5',vx,LAD,...
	'paramFile_SteeringMdl_CarMaker_DSR');

sys_d = ssMdl_SingleTrack('stdsr','paramFile_SingleTrackMdl_BMW5',vx,LAD,...
	'paramFile_SteeringMdl_CarMaker_DSR');

sys_e = ssMdl_SingleTrack('stvisdsr','paramFile_SingleTrackMdl_BMW5',vx,LAD,...
	'paramFile_SteeringMdl_CarMaker_DSR');

sys_f = ssMdl_SingleTrack('stvisdsr_2int','paramFile_SingleTrackMdl_BMW5',vx,LAD,...
	'paramFile_SteeringMdl_CarMaker_DSR');
