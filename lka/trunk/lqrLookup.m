% LQR-Regler für Lookup-Tables
% 
% Subject: lka
% Author: georgnoname
% Date: 04.12.2012 - 12.03.2013

clear all

% vehicle parameter
vehicleParameter = 'paramFile_SingleTrackMdl_BMW5';

% range of longitudinal velocity vx
vxStart = 5;
vxStop = 50;
vxDelta = 5;
vxRow = vxStart:vxDelta:vxStop;

% range of look-ahead distance lad
ladStart = 5;
ladStop = 30;
ladDelta = 5;
ladCol = ladStart:ladDelta:ladStop;


for i = 1:length(ladCol)
    for j = 1:length(vxRow)

        contr.t = lkaController_t(vehicleParameter,vxRow(j),ladCol(i));
        k = contr.t.LQR_2int.k;
        
        K1(j,i) = k(1);
        K2(j,i) = k(2);
        K3(j,i) = k(3);
        K4(j,i) = k(4);
        K5(j,i) = k(5);
        K6(j,i) = k(6);
        
    end%for
end%for

