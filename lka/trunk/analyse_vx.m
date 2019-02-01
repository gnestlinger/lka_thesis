% LKA-Analyse: vx variabel, lad konstant
% 
% .) des um Relativposition erw. Einspurmodells
% 
% Subject: 
% Author: georgnoname
% Date: 07.03.2013

clc
clear;

figure;
hold all;

lad = 5;
i = 0;
w = logspace(-1,2,100);

params = paramFile2Struct('paramFile_SingleTrackMdl_BMW5');
VX = 5:5:50;
for vx = VX
    i = i+1;
    
    sys = ssMdl_SingleTrack('stvis',params,vx,lad);
    sys = sys(3,1);
	
    bode(sys,w);
    [mag,phase] = bode(sys,w);
    magdb(:,i) = 20*log10(mag(:));
    Phase(:,i) = phase(:);
    
    
end%for

grid on
legend(cellfun(@num2str,num2cell(VX),'UniformOutput',false));

% saveFigure_txt('EinspurMdlVis_Bode_BMW5_vxVar',[],...
%     'w_rad/s',w,...
%     'abs_vx5',magdb(:,1),'arg_vx5',Phase(:,1),...
%     'abs_vx10',magdb(:,2),'arg_vx10',Phase(:,2),...
%     'abs_vx15',magdb(:,3),'arg_vx15',Phase(:,3),...
%     'abs_vx20',magdb(:,4),'arg_vx20',Phase(:,4),...
%     'abs_vx25',magdb(:,5),'arg_vx25',Phase(:,5),...
%     'abs_vx30',magdb(:,6),'arg_vx30',Phase(:,6),...
%     'abs_vx35',magdb(:,7),'arg_vx35',Phase(:,7),...
%     'abs_vx40',magdb(:,8),'arg_vx40',Phase(:,8),...
%     'abs_vx45',magdb(:,9),'arg_vx45',Phase(:,9),...
%     'abs_vx50',magdb(:,10),'arg_vx50',Phase(:,10))


