function [u,yL,epsL,kapL] = einspurMdl_u(t,x,pIn,fh)
% einspurMdl_u     Eingangsgröße für das Lineare Einspurmodell
%   _______
%   Syntax:
%   [u,yL,epsL,kapL] = einspurMdl_u(t,x,pIn,fh)
%   ________________
%   Input arguments:
%   t ..... time
%   x ..... state
%   pIn ... structure of vehicle parameters stored in pIn.vehicle.values
%   fh .... control input (function handle/anonymus function)
% 
% Subject: lka
% Author: georgnoname
% Date: 01.12.2012 - 22.02.2013


% lane tracking
[~,yL,epsL,kapL] = lkaLaneTracking([x(5),x(6),x(3)],pIn.traj,pIn.lad);

% control input
u = fh(t,x,yL,epsL);

end%fcn

