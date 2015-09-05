


%%% dateinamen anlegen
filename = 'lkaSim_IPG_mitLenkung';

%%% Sollbahn
% filename = [filename,'_SprungR200_LQR2int_ZRMdl'];
% filename = [filename,'_RampeR200_LQR2int_ZRMdl'];
filename = [filename,'_Circuit01_LQR2int_ZRMdl'];

%%% Regler verwendet vy oder nicht
filename = [filename,'_vy0'];
% filename = [filename,'_vy1'];

%%% Parametervariation
% filename = [filename,'_mu0p4'];
% filename = [filename,'_paramVar'];


% Daten des verwendeten Reglers
solIPG.simIn.Controller.lka = lkaContr.t.LQR_DSR_2int.fix;

% Sim-Output IPG
solIPG.simOut.t = ipgSimout_v.time;
solIPG.simOut.v = ipgSimout_v.signals.values;
solIPG.simOut.ay = ipgSimout_ay.signals.values;
solIPG.simOut.beta = ipgSimout_beta.signals.values;
solIPG.simOut.lad = ipgSimout_lad.signals.values;
solIPG.simOut.xg = ipgSimout_xg.signals.values;
solIPG.simOut.yg = ipgSimout_yg.signals.values;
solIPG.simOut.sg = ipgSimout_sg.signals.values;
solIPG.simOut.yL = ipgSimout_yL.signals.values;
solIPG.simOut.epsL = ipgSimout_epsL.signals.values;
solIPG.simOut.kapL = ipgSimout_kapL.signals.values;
solIPG.simOut.TorqueH = ipgSimout_TorqueH.signals.values;
solIPG.simOut.vy = ipgSimout_vy.signals.values;
solIPG.simOut.vy_contrIn = ipgSimout_vy_contrIn.signals.values;



% speichere als mat-File
save(filename,'solIPG');

% speichere als ASCII-File
saveFigure_txt(filename,500,...
    't',solIPG.simOut.t,...
    'v',solIPG.simOut.v,...
    'ay',solIPG.simOut.ay,...
    'beta',solIPG.simOut.beta,...
    'LAD',solIPG.simOut.lad,...
    'xg',solIPG.simOut.xg,...
    'yg',solIPG.simOut.yg,...
    'sg',solIPG.simOut.sg,...
    'yL_SensorIPG',solIPG.simOut.yL,...
    'epsL_SensorIPG',solIPG.simOut.epsL,...
    'kapL_SensorIPG',solIPG.simOut.kapL,...
    'TorqueH',solIPG.simOut.TorqueH,...
    'vy',solIPG.simOut.vy,...
    'vy_contrIn',solIPG.simOut.vy_contrIn);
