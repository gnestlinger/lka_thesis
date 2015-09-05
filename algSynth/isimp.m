% ISIMP  Überprüfung auf Implementierbarkeit
%
% ===============================================================
% M. Horn, N. Dourdoumas: Regelungstechnik, Pearson-Studium, 2004
%                         Kapitel 18 und 19
% ===============================================================
%
% flag=isimp(P,T)  überprüft, ob bei gegebener Strecke mit der
% Übertragungsfunktion P die Führungsübertragungsfunktion T 
% implementierbar ist.
%
% P ... Streckenübertragungsfunktion (Datentyp tf)
% T ... gewünschte Führungsübertragungsfunktion T (Datentyp tf)
%
% flag=1 ... T ist implementierbar
% flag=0 ... T ist nicht implementierbar
%                                
% Hinweis: im Falle von z-Übertragungsfunktionen müssen P und T die
%          gleiche Diskretisierungszeit besitzen, sonst ist flag=0.
%
% (c) 2003 M.Horn

function flag=isimp(P,T)


flag=isimp_c(P,T);
