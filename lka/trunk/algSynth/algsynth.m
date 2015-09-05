% ALGSYNTH   algebraischer Entwurf für die erweiterte Regelkreisstruktur 
%
% ===============================================================
% M. Horn, N. Dourdoumas: Regelungstechnik, Pearson-Studium, 2004
%                         Kapitel 18.5 und 19.2
% ===============================================================
%
% [R,V]=algsynth(P,T)    ermittelt bei gegebener Strecke mit der
% Übertragungsfunktion P den gesuchten Regler [R,V] so, dass der
% Regelkreis die Führungsübertragungsfunktion T besitzt.
%
% P ... Streckenübertragungsfunktion (Datentyp tf)
% T ... gewünschte Führungsübertragungsfunktion T (Datentyp tf)
%
% R,V ... gesuchte Reglerübertragungsfunktionen (Datentyp tf)
%                                
% [R,V]=algsynth(P,T,arg,val) algebraische Synthese mit
% Zusatzwünschen für das Reglernennerpolynom a, d.h.
%
% a(arg(1))=val(1),  a(arg(2))=val(2),  a(arg(3))=val(3) ...
%
% [R,V]=algsynth(P,T,arg) hier wird val=[0 0 0 ...] gesetzt, d.h.
%
% a(arg(1))=0,  a(arg(2))=0,  a(arg(3))=0 ...
%
% siehe auch: POLVOR
%
% (c) 3.12.2003 M.Horn

function [R,V]=algsynth(p1,p2,p3,p4)

% Init
R=[]; V=[];

% Aufrufmodus abfragen
switch nargin
    case 2,       
        % [R,V]=algsynth(P,T)     
        P=p1; T=p2;
        [R,V]=algsynth_c(P,T); 
    case 3,    
        % [R,V]=algsynth(P,T,arg)            
        P=p1; T=p2; arg=p3; val=zeros(1,length(p3));
        [R,V]=algsynth_c(P,T,arg);  
    case 4,
        % [err,R,V]=algsynth(P,T,arg,val)       
        P=p1; T=p2; arg=p3; val=p4;      
        [R,V]=algsynth_c(P,T,arg,val);  
    otherwise
        error('Bitte Eingabe überprüfen - siehe  HELP ALGSYNTH')
end