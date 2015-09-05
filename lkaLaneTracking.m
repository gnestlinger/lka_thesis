function [out,yL,epsL,kapL] = lkaLaneTracking(in,traj,lad)
% lkaLaneTracking  Berechnung von yL, epsL, kapL
%   _______
%   Syntax: [out,yL,epsL,kapL] = lkaLaneTracking(in,traj,lad)
%   ________________
%   Input arguments:
%   in ..... vector of global vehicle position and orientation
%   traj ... structure of intended trajectory
%   lad .... look-ahead distance [m]
%   ________________
%   Output arguments:
%   out .... [yL;epsL;kapL;]
%   yL ..... Querversatz in der Entfernung lad vor dem Fahrzeug [m]
%   epsL ... Relativwinkel in der Entfernung lad vor dem Fahrzeug [rad]
%   kapL ... Sollbahnkrümmung in der Entfernung lad vor dem Fahrzeug [1/m]
% 
% Subject: lka
% Author: georgnoname
% Date: 26.09.2012 - 30.04.2013

% Methode: Koordinatentransformation
% struct g ... globales Koordinatensystem (Index 0 in Funktion doPlot)
% struct V ... glob. Koord.system rotiert sodass Fzg. parallel zu x-Achse

% globale Fzg.position (CG) und Orientierung
% g.CG.x = in(1);
% g.CG.y = in(2);
% g.psi = in(3);

% Trajektorie (global)
% g.traj = traj;

% Koordinatentransformation: Drehmatrix mit in(3) = psi
M = [cos(-in(3)) -sin(-in(3));sin(-in(3)) cos(-in(3))];

% CG (transformiert)
V.CG.x = M(1,:)*[in(1); in(2)];
V.CG.y = M(2,:)*[in(1); in(2)];

% Solltrajektorie (transformiert)
trajt = [traj.x, traj.y];
V.traj.x = M(1,:)*trajt';
V.traj.y = M(2,:)*trajt';

% Koordinaten des Punkts at look-ahead distance (transformiert)
V.L.x = V.CG.x + lad;
V.L.y = V.CG.y;

% maximaler Abstand zwischen x-Werten
deltax = max(abs(diff(V.traj.x)));

% binärwertige Indizes für Elemente aus 'V.traj.x' im Bereich von 'V.L.x'
logIndx1 = V.L.x - deltax/2 <= V.traj.x;
logIndx2 = V.L.x + deltax/2 >= V.traj.x;
logIndx = logIndx1 & logIndx2;

% error if no element of 'logIndx' is logical 1
if ~any(logIndx)
    doPlot(in,traj,lad,V,M);
    error(['visionSystem.m: Keine Elemente der Sollbahn im Bereich der',... 
        ' aktuellen Fahrzeugposition gefunden'])
end%if

% explizite Indizes für Elemente aus 'V.traj.x' im Bereich von 
% V.L.x - deltax/2 < V.traj.x < V.L.x + deltax/2
indx = find(logIndx); % indx = traj.ind(logIndx);

% Inter- bzw. Extrapoliere y-Werte der transf. Solltrajektorie
yTraj = zeros(length(indx),1);
try
    for i = 1:length(indx)
        [indl,indu] = interpIndex(indx(i),1,traj.ind(end));
%         yTraj(i) = interp1(V.traj.x(indl:indu),V.traj.y(indl:indu),...
%             V.L.x,'spline');
        % same result like interp1(..,'spline') but faster
        yTraj(i) = spline(V.traj.x(indl:indu),V.traj.y(indl:indu),V.L.x);
    end%for
catch exception
    disp(exception.message);
    doPlot(in,traj,lad,V,M);
    hold on
    plot([V.L.x - deltax/2,V.L.x - deltax/2],[V.L.y-10,V.L.y+30],'-k');
    plot([V.L.x + deltax/2,V.L.x + deltax/2],[V.L.y-10,V.L.y+30],'-k');
end%try

% Querabstand von interpolierten y-Werten der Trajektorie zu
% Fahrzeuglaengsachse
yL = yTraj - V.CG.y;

% wähle "wahrscheinlichsten" (betragsmäßig kleinsten) Wert aus, falls > 1
% y-Wert zu einem x-Wert existiert (zb. bei geschlossener Solltrajektorie)
if length(indx) < 2
    minInd = 1;
else
    [~,minInd] = min(abs(yL));
end%if

% output argument yL
yL = yTraj(minInd) - V.CG.y;

% refresh interpolating-indices
[indl,indu] = interpIndex(indx(minInd),1,traj.ind(end));

% Tangentenvektor
tangent = [V.traj.x(indu)-V.traj.x(indl); V.traj.y(indu)-V.traj.y(indl)];

% output argument epsL
epsL = atan2(tangent(2),tangent(1));

% output argument rInv: read appropriate element from 'traj.k'
% kapL = traj.k(indx(minInd));
% interpoliere für glatte Verläufe bei Regelung
kapL = spline(V.traj.x(indl:indu),traj.k(indl:indu),V.L.x);

% collect all output arguments (to be used in simulink)
out = [yL;epsL;kapL];

end%fcn



function [indl,indu] = interpIndex(ind,indMin,indMax)
% return the indices used to interpolate, basically ind-1 and ind+1 but do
% some error checking (if-blocks)

% lower interpolating-index
indl = ind-1; 
if indl < indMin; 
    indl = indl+1; 
end%if

% upper interpolating-index
indu = ind+1; 
if indu > indMax; 
    indu = indu-1; 
end%if

end%fcn



function doPlot(in,traj,lad,V,M)

x0CG = in(1);
y0CG = in(2);
psi = in(3);

% Solltrajektorie (global)
plot(traj.x,traj.y,'o','MarkerSize',3,'MarkerFaceColor','b'); 
hold on

% Fahrzeug
plot(x0CG,y0CG,'ob','MarkerSize',10,'MarkerFaceColor','b');
Fzglachse = [x0CG+lad*cos(psi); y0CG+lad*sin(psi)];
plot([x0CG,Fzglachse(1)],[y0CG,Fzglachse(2)],'b','LineWidth',2)

% Solltrajektorie (transformiert)
plot(V.traj.x,V.traj.y,'ro','MarkerSize',3,'MarkerFaceColor','r');

% Fahrzeug (transformiert)
plot(V.CG.x,V.CG.y,'or','MarkerSize',10,'MarkerFaceColor','r');
Fzglachset = M*Fzglachse;
plot([V.CG.x,Fzglachset(1)],[V.CG.y,Fzglachset(2)],'r','LineWidth',2)

hold off
grid on
axis equal

end%fcn
