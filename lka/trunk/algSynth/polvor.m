% POLVOR   Polvorgabe mit Hilfe des Standardregelkreises
%
% ===============================================================
% M. Horn, N. Dourdoumas: Regelungstechnik, Pearson-Studium, 2004
%                         Kapitel 18.4 und 18.6
% ===============================================================
%
% [R,T]=polvor(P,nut)  ermittelt bei gegebener Strecke mit
% der �bertragungsfunktion P den gesuchten Regler R so, dass die
% F�hrungs�bertragungsfunktion T das Nennerpolynom nut besitzt
% 
%
% P ..... Strecken�bertragungsfunktion (Datentyp tf)
% nut ... gew�nschtes Nennerpolynom von T
%
% R ..... gesuchte Regler�bertragungsfunktion R (Datentyp tf)
% T ..... F�hrungs�bertragungsfunktion T=minreal(R*P/(1+R*P))
%
% [R,T]=polvor(P,nut,arg,val) Polvorgabe mit Zusatzw�nschen 
% f�r das Reglernennerpolynom a, d.h.
%
% a(arg(1))=val(1),  a(arg(2))=val(2),  a(arg(3))=val(3) ...
%
% [R,T]=polvor(P,nut,arg) hier wird val=[0 0 0 ...] gesetzt, d.h.
%
% a(arg(1))=0,  a(arg(2))=0,  a(arg(3))=0 ...
%
% siehe auch: ALGSYNTH
%
% (c) 2003 M.Horn

function [R,T]=polvor(p1,p2,p3,p4)


% Aufrufmodus abfragen
switch nargin
    case 2,       
        % [R,T]=polvor(P,nut)        
        P=p1; nut=p2;     
        [R,T]=polvor_c(P,nut);
    case 3,    
        % [R,T]=polvor(P,nut,arg)                
        P=p1; nut=p2; arg=p3; val=zeros(1,length(p3));
        [R,T]=polvor_c(P,nut,val);
    case 4,
        % [R,T]=polvor(P,nut,arg,val)          
        P=p1; nut=p2; arg=p3; val=p4;      
        [R,T]=polvor_c(P,nut,arg,val);
    otherwise
        error('Bitte Eingabe �berpr�fen - siehe  HELP POLVOR')
end