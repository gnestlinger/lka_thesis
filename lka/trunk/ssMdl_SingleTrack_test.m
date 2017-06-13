



initFile_


sys_a = ssMdl_SingleTrack('st','paramFile_SingleTrackMdl_BMW5',20,10);

sys_b = ssMdl_SingleTrack('stvis','paramFile_SingleTrackMdl_BMW5',20,10);

sys_c = ssMdl_SingleTrack('stvis_2int','paramFile_SingleTrackMdl_BMW5',20,10,...
	'paramFile_SteeringMdl_CarMaker_DSR');

sys_d = ssMdl_SingleTrack('stdsr','paramFile_SingleTrackMdl_BMW5',20,10,...
	'paramFile_SteeringMdl_CarMaker_DSR');

sys_e = ssMdl_SingleTrack('stvisdsr','paramFile_SingleTrackMdl_BMW5',20,10,...
	'paramFile_SteeringMdl_CarMaker_DSR');

sys_f = ssMdl_SingleTrack('stvisdsr_2int','paramFile_SingleTrackMdl_BMW5',20,10,...
	'paramFile_SteeringMdl_CarMaker_DSR');
