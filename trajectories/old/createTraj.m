function [out,varargout] = createTraj(type,varargin)
% createTraj    create or connect trajectory segment(s)
%   _______
%   Syntax:
%   [out,varargout] = createTraj(type,varargin)
%   ________________
%   Input arguments:
%   type ....... one of the strings: straight/curv/cloth/cloth2/connect
%   varargin ... depends on input argument type (see corresponding subfcn.)
%   _________________
%   Output arguments:
%   out ......... structure of trajectory segment
%   varargout ... error index in case of type = 'connect'
% 
% 
%   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   createTrajSegm('straight',xyStart,xyStop,delta) %%%%%%%%%%%%%%%%%%%%%%%
%   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Draws a straight segment from 'xyStart' to 'xyStop'.
% 
%   xyStart ... Startpunkt [m]        
%   xyStop .... Endpunkt [m]
%   delta ..... Abstand zwischen zwei Punkten (opt.) [m]
% 
% 
%   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   createTrajSegm('curv',xyStart,angleStart,angleStop,radius,delta) %%%%%%
%   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   draws a circle of radius 'radius' 
%     .) clockwise if 'angleStart' > 'angleStop' with curvature < 0
%     .) counter-clockwise if 'angleStart' < 'angleStop' with curvature > 0
% 
%   xyStart ...... Startpunkt [m]
%   angleStart ... Startwinkel im Punkt xyStart [rad]
%   angleStop .... Endwinkel [rad]    
%   radius ....... Kurvenradius (>0) [m]
%   delta ........ Soll-Bogenlänge zwischen zwei Punkten (opt.) [m]
% 
% 
%   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   createTrajSegm('cloth',xyStart,slopeStart,sStop,A,delta) %%%%%%%%%%%%%%
%   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   xyStart ...... Startpunkt [m]
%   slopeStart ... Steigung der Klothoide im Startpunkt xyStart [rad]
%   A ............ Klothoidenparameter r = A^2/s [1]
%                  > 0 ... Krümmung gegen Uhrzeigersinn
%                  < 0 ... Krümmung mit dem Uhrzeigersinn
%   sStop ........ Bogenlänge der Klothoide [m]
%                  > 0 ... Verlauf von Krümmung 0 -> 1/r = sStop/A^2
%                  < 0 ... Verlauf von Krümmung 1/r = sStop/A^2 -> 0
%   delta ........ Soll-Bogenlänge zwischen zwei Punkten (opt.) [m]
%
% 
%   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   createTrajSegm('cloth2',xyStart,curvStart,curvStop,A,slopeStart,delta) 
%   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   draws a clothoid with start-curvature 'curvStart' and end-curvature
%   'curvStop'
%     .) clockwise if 'curvStart' > 'curvStop' with curvature < 0
%     .) counter-clockwise if 'curvStart' < 'curvStop' with curvature > 0
%   if A > 0 and 
%     .) clockwise if 'curvStart' < 'curvStop' with curvature < 0
%     .) counter-clockwise if 'curvStart' > 'curvStop' with curvature > 0
%   if A < 0.
% 
%   xyStart ..... Startpunkt [m]
%   curvStart ... Krümmung der Klothoide im Startpunkt xyStart [1/m]
%   curvStop .... Krümmung der Klothoide im Endpunkt [1/m]
%   A ........... Klothoidenparameter r = A^2/s [1]
%                 > 0 ... Krümmung gegen Uhrzeigersinn
%                 < 0 ... Krümmung mit dem Uhrzeigersinn
%   slopeStart .. Steigung der Klothoide im Startpunkt xyStart [rad]
%   delta ....... Bogenlänge zwischen zeii Punkten [m]
%
% 
% Subject: lka
% Author: georgnoname
% Date: 24.09.2012 - 19.05.2013


% check input arguments
if ~ischar(type); error('inputArg type should be string'); end%if

% Standardwert für Abstand zwischen zwei benachbarten Punkten einer
% Trajektorie
deltaStd = 0.2;

% switch to subfunction
switch lower(type)
    
    % straight segment
    case 'straight'
        
        % check delta
        indDelta = 3;
        if size(varargin,2) < indDelta; % set delta if unspecified
            varargin{indDelta} = deltaStd;
        else
            if isempty(varargin{end})	% set delta if empty
                varargin{indDelta} = deltaStd;
            end%if
        end%if
        
        % check number of input arguments
        if (size(varargin,2) < 2) || (size(varargin,2) > indDelta);
            error('Nbr. of input arguments does not match function type.');
        end%if
        ret = straight(varargin{1},varargin{2},varargin{3});
        
        
    % circular segment
    case 'curv'
        
        % check delta
        indDelta = 5;
        if size(varargin,2) < indDelta; % set delta if unspecified
            varargin{indDelta} = deltaStd;
        else
            if isempty(varargin{end})	% set delta if empty
                varargin{indDelta} = deltaStd;
            end%if
        end%if
        
        % check number of input arguments
        if (size(varargin,2) < 3) || (size(varargin,2) > indDelta);
            error('Input arguments do not match function type.');
        end%if
        ret = curvature(varargin{1},varargin{2},varargin{3},varargin{4},...
            varargin{5});
        
        
    % clothoid segment type 1
    case 'cloth'
        
        % check delta
        indDelta = 5;
        if size(varargin,2) < indDelta; % set delta if unspecified
            varargin{indDelta} = deltaStd;
        else
            if isempty(varargin{end})	% set delta if empty
                varargin{indDelta} = deltaStd;
            end%if
        end%if
        
        % check number of input arguments
        if (size(varargin,2) < 3) || (size(varargin,2) > indDelta);
            error('Input arguments do not match function type.');
        end%if
        ret = clothoid(varargin{1},varargin{2},varargin{3},varargin{4},...
            varargin{5});
        
        
    % clothoid segment type 2
    case 'cloth2'
        
        % check delta
        indDelta = 6;
        if size(varargin,2) < indDelta; % set delta if unspecified
            varargin{indDelta} = deltaStd;
        else
            if isempty(varargin{end})	% set delta if empty
                varargin{indDelta} = deltaStd;
            end%if
        end%if
        
        % check number of input arguments
        if (size(varargin,2) < 5) || (size(varargin,2) > indDelta);
            error('Input arguments do not match function type.');
        end%if
        ret = clothoid2(varargin{1},varargin{2},varargin{3},varargin{4},...
            varargin{5},varargin{6});
        
        
    % connect segments
    case 'connect'
        [ret,err] = connectTraj(varargin);
        varargout{1} = err;
        
    otherwise
        error('Unknown type')
        
end%switch


out = ret;

end%fcn


function ret = straight(xyStart,xyStop,delta)
% straight trajectory (type 0)
% 
% xyStart ..... Startpunkt [m]
% xyStop ...... Endpunkt [m]
% delta ....... Abstand zwischen zwei Punkten (opt.) [m]
% 
% Subject: lka
% Author: georgnoname
% Date: 24.09.2012 - 09.01.2013

    % check input arguments
    if isempty(xyStart); xyStart = [0;0]; end
    if numel(xyStart) ~= 2; error('Dimension of xyStart'); end
    if numel(xyStop) ~= 2; error('Dimension of xyStop'); end
    if nargin < 3; delta = 0.2; end % obsolet
    
    % ensure column-vectors
    xyStart = col(xyStart);
    xyStop = col(xyStop);
    
    % compute number of segments
    l = norm(xyStop-xyStart);
    n = ceil(l/delta) + 1;
    
    % output arguments
    ret.x = col(linspace(xyStart(1),xyStop(1),n));
    ret.y = col(linspace(xyStart(2),xyStop(2),n));
    ret.s = sqrt((ret.x-xyStart(1)).^2 + (ret.y-xyStart(2)).^2);
    ret.k = zeros(size(ret.x));
    ret.type = 0*ones(size(ret.x));
    
    % test
%     plot(x,y,'o',...
%         'MarkerSize',2,...
%         'MarkerFaceColor','blue');
%     grid on;

end%fcn


function ret = curvature(xyStart,angleStart,angleStop,radius,delta)
% circular trajectory
% 
% draws a circle of radius 'radius' 
%   .) clockwise if 'angleStart' > 'angleStop' with curvature < 0
%   .) counter-clockwise if 'angleStart' < 'angleStop' with curvature > 0
% 
% xyStart ..... Startpunkt [m]
% angleStart .. Startwinkel im Punkt xyStart [rad]
% angleStop ... Endwinkel [rad]    
% radius ...... Kurvenradius (>0) [m]
% delta ....... Soll-Bogenlänge zwischen zwei Punkten (opt.) [m]
% 
% Subject: lka
% Author: georgnoname
% Date: 24.09.2012 - 09.01.2013

    % check input arguments
    if isempty(xyStart); xyStart = [0;0]; end
    if numel(xyStart) ~= 2; error('Dimension of xyStart'); end
    if radius <= 0; error('Input argument radius <= 0'); end
    if nargin < 5; delta = 0.2; end % obsolet
    
    % ensure column-vectors
    xyStart = xyStart(:);
    
    % get sign of curvature k !!!!!!! liefert vlt. FALSCHE Ergebnisse
    if angleStop < angleStart; signk = -1; else signk = 1; end
    
    % create curvature trajectory
    alpha = delta/radius;    % Soll-Winkel des Bogensegments [rad]
    k = ceil(abs((angleStop-angleStart)/alpha));  % Anzahl der Bogensegmente
    phi = linspace(angleStart,angleStop,k+1)';
    x = radius*cos(phi);
    y = radius*sin(phi);
    
    % shift whole trajectory ([x(1);y(1)] matches xyStart)
    xShift = xyStart(1) - x(1);
    yShift = xyStart(2) - y(1);
    
    % output arguments x/y
    ret.x = x + xShift;
    ret.y = y + yShift;
    
    % output argument s/kappa/type
    ret.s = sort(abs((phi-phi(1))*radius),'ascend');
    ret.k = signk/radius*ones(size(x));
    ret.type = ones(size(ret.x));
    
    % test
%     plot(ret.x,ret.y,'o',...
%         'MarkerSize',2,...
%         'MarkerFaceColor','blue');
%     grid on;

end%fcn


function ret = clothoid(xyStart,slopeStart,sStop,A,delta)
% clothoid trajectory (type 1)
% 
% xyStart ..... Startpunkt [m]
% slopeStart .. Steigung der Klothoide im Startpunkt xyStart [rad]
% A ........... Klothoidenparameter r = A^2/s [1]
%               > 0 ... Krümmung gegen Uhrzeigersinn
%               < 0 ... Krümmung mit dem Uhrzeigersinn
% sStop ....... Bogenlänge der Klothoide [m]
%               > 0 ... Verlauf von Krümmung 0 -> 1/r = sStop/A^2
%               < 0 ... Verlauf von Krümmung 1/r = sStop/A^2 -> 0
% delta ....... Soll-Bogenlänge zwischen zwei Punkten (opt.) [m]
% 
% Subject: lka
% Author: georgnoname
% Date: 24.09.2012 - 09.01.2013

    % check input arguments
    if isempty(xyStart); xyStart = [0;0]; end
    if numel(xyStart) ~= 2; error('Dimension of xyStart'); end
    if nargin < 4; delta = 0.2; end % obsolet
    
    % get sign of clothoid parameter A
    signA = sign(A); 
    A = abs(A);
    
    % get sign of clothoid arc length 'sStop'
    signsStop = sign(sStop); 
%     sStop = abs(sStop);
    
    % clothoid integrand (parameterized form)
    intx = @(t) cos(pi.*t.^2./2);
    inty = @(t) sin(pi.*t.^2./2);
    
    % Constant
    sqpi = sqrt(pi);
    
    % num. integration
    s = col(0:signsStop*delta:sStop);
    l = s/(A*sqpi);
    for i = 1:length(l)
        cloth.x(i,1) = A*sqpi*quad(intx,0,l(i));
        cloth.y(i,1) = signA*A*sqpi*quad(inty,0,l(i));
    end%for
    
    % derivative to compute tangent vector (
    % http://mathworld.wolfram.com/TangentVector.html
    xD = col(A*sqpi*cos(pi*l.^2/2));
    yD = col(A*sqpi*sin(pi*l.^2/2));
    tang.x = xD./sqrt(xD.^2+yD.^2);
    tang.y = yD./sqrt(xD.^2+yD.^2);
    
    cloth.k = s/A^2*sign(A);
    if signsStop < 0
        cloth.x = flipud(cloth.x);
        cloth.y = flipud(cloth.y);
        cloth.k = flipud(cloth.k);
        tang.x = flipud(tang.x);
        tang.y = flipud(tang.y);
    end%if
    
    % rotate 
    M = [cos(slopeStart) -sin(slopeStart);sin(slopeStart) cos(slopeStart)];
    xyRotated = M*[cloth.x';cloth.y'];
    cloth.rot.x = col(xyRotated(1,:));
    cloth.rot.y = col(xyRotated(2,:));
    
    xyTangRotated = M*[tang.x';tang.y'];
    tang.rot.x = col(xyTangRotated(1,:));
    tang.rot.y = col(xyTangRotated(2,:));
        
    % shift whole trajectory (start point to xyStart)
    xShift = xyStart(1)-cloth.rot.x(1);
    yShift = xyStart(2)-cloth.rot.y(1);
    
    % output arguments
    ret.x = cloth.rot.x + xShift;
    ret.y = cloth.rot.y + yShift;
    ret.s = sort(abs(s),'ascend');
%     ret.k = s/A^2*sign(A);
    ret.k = cloth.k;
    ret.type = 2*ones(size(ret.x));
    ret.phi = angle(tang.rot.x + 1i*tang.rot.y) - pi/2;
    
    
%     % test
%     plot(x,y,'o',...
%         'MarkerSize',1.5,...
%         'MarkerFaceColor','blue');
%     hold on
%     t.x = tang.rot.x;
%     t.y = tang.rot.y;
%     quiver(ret.x,ret.y,t.x,t.y,0,'ob','MarkerSize',4);
%     % sStop -> inf
%     plot(A*sqpi/2,A*sqpi/2,'r.');
%     grid on; axis equal;

end%fcn


function ret = clothoid2(xyStart,curvStart,curvStop,A,slopeStart,delta)
% clothoid trajectory (type 2)
% 
% draws a clothoid with start-curvature 'curvStart' and end-curvature
% 'curvStop'
%   .) clockwise if 'curvStart' > 'curvStop' with curvature < 0
%   .) counter-clockwise if 'curvStart' < 'curvStop' with curvature > 0
% if A > 0 and 
%   .) clockwise if 'curvStart' < 'curvStop' with curvature < 0
%   .) counter-clockwise if 'curvStart' > 'curvStop' with curvature > 0
% if A < 0.
% 
% The curvature is given as
%   curv = s/A^2
% where 's' is the length of the curve and 'A' is a constant.
% 
% xyStart ..... Startpunkt [m]
% curvStart ... Krümmung der Klothoide im Startpunkt xyStart [1/m]
% curvStop .... Krümmung der Klothoide im Endpunkt [1/m]
% A ........... Klothoidenparameter r = A^2/s [1]
%               > 0 ... Krümmung gegen Uhrzeigersinn
%               < 0 ... Krümmung mit dem Uhrzeigersinn
% slopeStart .. Steigung der Klothoide im Startpunkt xyStart [rad]
% delta ....... Bogenlänge zwischen zeii Punkten [m]
% 
% Subject: lka
% Author: georgnoname
% Date: 24.09.2012 - 19.05.2013


    % check input arguments
    if isempty(xyStart); xyStart = [0;0]; end
    if numel(xyStart) ~= 2; error('Dimension of xyStart'); end
    if nargin < 6; delta = 0.2; end % obsolet
    if nargin < 5; slopeStart = 0; end
    if isempty(slopeStart); slopeStart = 0; end
    
    % get sign of clothoid parameter A
    signA = sign(A);
    A = abs(A);
    
    % 
    if curvStart > curvStop; signk = -1; else signk = 1; end
%     signk = signk*signA;
    
    % clothoid integrand (parameterized form)
    intx = @(t) cos(pi.*t.^2./2);
    inty = @(t) sin(pi.*t.^2./2);
    
    % Constant
    sqpi = sqrt(pi);
    
    % calc sStart/sStop
    sStart = A^2*curvStart;
    sStop = A^2*curvStop;
    
    % calculate number of segments within [sStart,sStop]
    n = abs(sStop-sStart)/delta;
        
    % pre-calculation
    s = linspace(sStart,sStop,ceil(n)+1);
    l = s/(A*sqpi)*signk;
    
    % pre-allocation
    length_l = length(l);
    cloth.x(length_l,1) = 0;
    cloth.y(length_l,1) = 0;
    
    % num. integration
    cloth.x(1) = A*sqpi*quad(intx,0,l(1));
    cloth.y(1) = signA*A*sqpi*quad(inty,0,l(1));
    for i = 2:length_l
        cloth.x(i) = cloth.x(i-1) + A*sqpi*quad(intx,l(i-1),l(i));
        cloth.y(i) = cloth.y(i-1) + signA*A*sqpi*quad(inty,l(i-1),l(i));
    end%for
    
    % derivative to compute tangent vector
    % http://mathworld.wolfram.com/TangentVector.html
    xD = col(A*sqpi*cos(pi*l.^2/2));
    yD = col(signA*A*sqpi*sin(pi*l.^2/2));
    tang.x = xD./sqrt(xD.^2+yD.^2);
    tang.y = yD./sqrt(xD.^2+yD.^2);
    
    % rotate clothoid
    diffAng = angle(tang.x(1) + 1i*tang.y(1)); % falls curvStart != 0
    slopeStart = slopeStart - diffAng;
    M = [cos(slopeStart) -sin(slopeStart);sin(slopeStart) cos(slopeStart)];
    xyRotated = M*[cloth.x';cloth.y'];
    cloth.rot.x = col(xyRotated(1,:));
    cloth.rot.y = col(xyRotated(2,:));
    
    % rotate tangent
    xyTangRotated = M*[tang.x';tang.y'];
    tang.rot.x = col(xyTangRotated(1,:));
    tang.rot.y = col(xyTangRotated(2,:));
        
    % shift whole trajectory ([x(1);y(1)] matches xyStart)
    xShift = xyStart(1)-cloth.rot.x(1);
    yShift = xyStart(2)-cloth.rot.y(1);
    
    % output arguments
    ret.x = cloth.rot.x + xShift;
    ret.y = cloth.rot.y + yShift;
    ret.sCloth = sort(abs(s));
    ret.s = ret.sCloth - ret.sCloth(1);
    ret.k = s/A^2*signk*signA;
    ret.type = 2*ones(size(ret.x));
    ret.phi = angle(tang.rot.x + 1i*tang.rot.y);
    
    
%     % test
%     plot(ret.x,ret.y,'o',...
%         'MarkerSize',1.5,...
%         'MarkerFaceColor','blue');
%     hold on
%     t.x = tang.rot.x;
%     t.y = tang.rot.y;
%     quiver(ret.x,ret.y,t.x,t.y,0,'ob','MarkerSize',4);
%     % sStop -> inf
%     plot(A*sqpi/2,A*sqpi/2,'r.');
%     grid on; axis equal;

end%fcn


function [ret,err] = connectTraj(varargin)
% connectTraj   connect trajectory segments
%   _______
%   Syntax:
%   [ret,err] = connectTraj(varargin)
%   ________________
%   Input arguments:
%   varargin ... segments created by subfunctions straight/curvature/..
%   like connectTraj(segm1,segm2,segm3,...)
%   _________________
%   Output arguments:
%   ret ... whole trajectory (connected segments) of type struct
%
% 05.10.2012 - 19.05.2013

varargin = varargin{1};
% err = [];

% get fieldnames of struct varargin{1}
fn = fieldnames(varargin{1});

% remove 'phi' and 'sCloth' from 'fn' if present (just at
% clothoid-segments)
keepIndex = true(size(fn));
removeIndex = [...
    strmatch('phi',fn,'exact'),...
    strmatch('sCloth',fn,'exact')];
keepIndex(removeIndex) = 0;
fn = fn(keepIndex);

% create first segment
% traj.x = col(varargin{1}.x);    % x-position
% traj.y = col(varargin{1}.y);    % y-position
% traj.s = col(varargin{1}.s);    % arc length
% traj.k = col(varargin{1}.k);    % curvature
for i = 1:length(fn);
    traj.(fn{i}) = col(varargin{1}.(fn{i}));
end

% append segments 2:end
for i = 2:size(varargin,2)
    traj.x = [traj.x; col(varargin{i}.x)];  % x-position
    traj.y = [traj.y; col(varargin{i}.y)];  % y-position
    traj.s = [traj.s; col(varargin{i}.s)+traj.s(end)];  % arc length
    traj.k = [traj.k; col(varargin{i}.k)];  % curvature
    traj.type = [traj.type; col(varargin{i}.type)]; % segment type
end%for

% finde x-Elemente mit Steigung 0 (führt zu Problemen mit
% interp1 @ visionSystem)
% logische 0 an Stelle wo Steigung von traj.x = 0 (x_werte identisch)
logIndx = logical(diff(traj.x));
logIndx = [logIndx;true(1)]; % Länge = Länge von traj.x/y/s/k

% finde y-Elemente mit Steigung 0 (führt zu Problemen mit
% interp1 @ visionSystem)
% logische 0 an Stelle wo Steigung von traj.x = 0 (x_werte identisch)
logIndy = logical(diff(traj.y));
logIndy = [logIndy;true(1)]; % Länge = Länge von traj.x/y/s/k

% Behalte Elemente, die Steigungn ~= 0 in x- ODER y-Richtung haben
logInd = logIndx | logIndy;
% logInd = logIndx;

% remove entries where diff(x) == 0
fn = fieldnames(traj);
for i = 1:length(fn)
    traj.(fn{i}) = traj.(fn{i})(logInd);
end

% output argument
traj.ind = col(1:length(traj.x));
ret = traj;

% check covered vehicle distance
indxErr = find(diff(traj.s)<0);
if isempty(indxErr)
    err = 0;
else 
    err = indxErr;
end%if

end%fcn


function ret = col(in)
% COL column vector
% checks if 'in' is of dimension 1*x or x*1 (vector, not matrix) and
% returns column vector of in.

[a,b] = size(in);

if (a>1) && (b>1)
    error('Input Argument should be vector but is matrix.')
end

ret = in(:);

end%fcn