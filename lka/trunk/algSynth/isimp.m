% ISIMP  �berpr�fung auf Implementierbarkeit
%
% ===============================================================
% M. Horn, N. Dourdoumas: Regelungstechnik, Pearson-Studium, 2004
%                         Kapitel 18 und 19
% ===============================================================
%
% flag=isimp(P,T)  �berpr�ft, ob bei gegebener Strecke mit der
% �bertragungsfunktion P die F�hrungs�bertragungsfunktion T 
% implementierbar ist.
%
% P ... Strecken�bertragungsfunktion (Datentyp tf)
% T ... gew�nschte F�hrungs�bertragungsfunktion T (Datentyp tf)
%
% flag=1 ... T ist implementierbar
% flag=0 ... T ist nicht implementierbar
%                                
% Hinweis: im Falle von z-�bertragungsfunktionen m�ssen P und T die
%          gleiche Diskretisierungszeit besitzen, sonst ist flag=0.
%
% (c) 2003 M.Horn

function flag=isimp(P,T)


flag=isimp_c(P,T);
