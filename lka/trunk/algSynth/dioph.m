% DIOPH   Diophantische Gleichung als lineares Gleichungssystem anschreiben
%
% ===============================================================
% M. Horn, N. Dourdoumas: Regelungstechnik, Pearson-Studium, 2004
%                         Kapitel 18
% ===============================================================
%
% [R,f]=dioph(mu,nu,nut) liefert die Resultante R und die rechte Seite f
% des linearen Gleichungssystems R*x=f, das sich aus
%
%        nu*a+mu*b=nut
%
% ergibt. Für x gilt dabei: x=[a0 a1 ... b0 b1].
%
% nu,mu .... vorgegebene Polynome mit max{Grad(nu),Grad(mu)}=n
% nut ...... vorgegebenes Polynom mit Grad(nut)=m 
%
% siehe auch: ALGSYNTH, POLVOR

% (c) 2003 M.Horn

function [R,f]=dioph(mu,nu,nut)

% führende Nullen eliminieren
while(mu(1)==0)
    mu=mu(2:end);
end
while(nu(1)==0)
    nu=nu(2:end);
end
while(nut(1)==0)
    nut=nut(2:end);
end

% stürzen für den Aufbau der Resultante
mu=mu(:);   mu=flipud(mu);
nu=nu(:);   nu=flipud(nu);
nut=nut(:); nut=flipud(nut); f=nut;

% n und m ermitteln
n=max((length(nu)-1),(length(mu)-1));
m=length(nut)-1;
if (m<n) 
    error('Problem prinzipiell nicht lösbar.')    
end

% Ordnung von a und b ermitteln
r=m-n;

% leere Resultante 
R=zeros(m+1,2*r+2);

% Gleichungssystem aufstellen

% Block für Koeffizienten von nu
col=zeros(m+1,1);
col(1:length(nu))=nu;
row=zeros(1,r+1);
row(1)=col(1);
T1=toeplitz(col,row);

% Block für Koeffizienten von mu
col=zeros(m+1,1);
col(1:length(mu))=mu;
row=zeros(1,r+1);
row(1)=col(1);
T2=toeplitz(col,row);

% Resultante
R=[T1 T2];