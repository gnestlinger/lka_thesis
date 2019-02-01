% LKA-Analyse: vx konstant, LAD variabel
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

vx = 20;
i = 0;
w = logspace(-1,2,100);

params = paramFile2Struct('paramFile_SingleTrackMdl_BMW5');
LAD = [10,15,20,30];
for lad = LAD
    i = i+1;
    
    sys = ssMdl_SingleTrack('stvis',params,vx,lad);
	sys = sys(3,1);
    
    bode(sys,w);
    [mag,phase] = bode(sys,w);
    magdb(:,i) = 20*log10(mag(:));
    Phase(:,i) = phase(:);
    
end%for

grid on
legend(cellfun(@num2str,num2cell(LAD),'UniformOutput',false));


% i = 0;
% string = [];
% 
% for lad = 0:5:50
%     i = i+1;
%     
%     string = [string,...
%         '''abs_lad',num2str(lad),''',magdb(:,',num2str(i),'),',...
%         '''arg_lad',num2str(lad),''',Phase(:,',num2str(i),'),'];
%     
% end

% saveFigure_txt('EinspurMdlVis_Bode_BMW5_ladVar',[],...
%     'w_rad/s',w,...
%     'abs_lad0',magdb(:,01),'arg_lad0',Phase(:,01),...
%     'abs_lad5',magdb(:,02),'arg_lad5',Phase(:,02),...
%     'abs_lad10',magdb(:,03),'arg_lad10',Phase(:,03),...
%     'abs_lad15',magdb(:,04),'arg_lad15',Phase(:,04),...
%     'abs_lad20',magdb(:,05),'arg_lad20',Phase(:,05),...
%     'abs_lad25',magdb(:,06),'arg_lad25',Phase(:,06),...
%     'abs_lad30',magdb(:,07),'arg_lad30',Phase(:,07),...
%     'abs_lad35',magdb(:,08),'arg_lad35',Phase(:,08),...
%     'abs_lad40',magdb(:,09),'arg_lad40',Phase(:,09),...
%     'abs_lad45',magdb(:,10),'arg_lad45',Phase(:,10),...
%     'abs_lad50',magdb(:,11),'arg_lad50',Phase(:,11))

