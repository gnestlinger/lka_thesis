classdef segDat
%SEGDAT 	Street segment data.
%	
%	OBJ = SEGDAT(X,Y,S,K,PHI,TYPE,NBR) creates the street segment object
%	OBJ of class SEGDAT representing the geometric data of a street segment
%	using the properties listed below.
%	
%	A curvature K > 0 (K < 0) indicates a left (right) turn by moving from
%	[X(1),Y(1)] -> [X(end),Y(end)].
%	
%	NBR is usually 1, but becomes interesting when street segments are
%	connected together. Then NBR represents the number of the street
%	segment in the path of connected segments.
%	
%	SEGDAT Properties:
%	 X		- x-coordinate [m].
%	 Y		- y-coordinate [m].
%	 S		- Length measured from the starting point [x(1);y(1)] [m].
%	 K		- Curvature [1/m].
%	 PHI	- Angle beetween the positive x-axis and the tangent [rad].
%	 TYPE	- Integer indicating the street segment type.
%	 NBR	- Street segment number of connected segments.
%	
%	SEGDAT Methods:
%	 changeSignOfCurvature - Change street segments curvature sign.
%	 plus	- Connect street segments using '+'.
%	 reverseDirection - Reverse street segment direction.
%	 rotate - Rotate street segment.
%	 shift	- Shift street segment.
%	 plot	- Plot the street segment.
%	 plotdiff - Plot the street segment with specific appearance.
%	 plottangent - Plot the street segment and specified tangents.
%	 
%	
%	See also LKASEGMENT, LKASEGMENTSTRAIGHT, LKASEGMENTCIRCLE,
%	LKASEGMENTCLOTHOID.
%

% DEVELOPMENT NOTES:
%	(1) Add methods for plot labeling (x/y labels, title, ...). DONE
%	(2) Show the clothoid direction in plots.

% Subject: lka
% Author: $Author$
% Date: $LastChangedDate$
% Revision: $Revision$



	properties (SetAccess = private)
		% x-coordinate [m].
		x
		
		% y-coordinate [m].
		y
		
		% Covered distance [m].
		s
		
		% Curvature [1/m].
		k
		
		% Tangent angle [rad].
		phi
		
		% Street segment type [-]:
		%	0 .. straight
		%	1 .. circular
		%	2 .. clothoid
		type
		
		% Street segment number of connected segments [-].
		nbr
	end%properties
	
	
	properties (Hidden)
		
		%%% plot specifications
		
		% Plot color
		plotColor = {'b','g','r'};
		
		% Plot marker symbol/size
		plotMarker = {...
			'o','diamond';... % marker symbol
			5,4}; % marker size
		
	end%properties
	
	
	
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	%%% METHODS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	
	%%% CONSTRUCTOR
	methods
		
		function obj = segDat(x,y,s,k,phi,type,nbr)
		%SEGDAT 	Create an instance of the SEGDAT object.
		
		
			% all inputs have to be non-empty column vectors, 
			[m,n] = size(x);
			
			% expand the lenght of TYPE/NBR if they are scalar
			if isscalar(type) && isscalar(nbr)
				type = type*ones(m,n);
				nbr = nbr*ones(m,n);
			end%if
			
			% just check X and check the size of all others against X
			if (m<1) || (n>1)
				error('SEGDAT:segDat:nonemptyColumnVectors',...
					'Inputs have to be non-empty column vectors')
			end%if
			if ~isequal(size(x),size(y),size(s),size(k),size(phi),size(type),size(nbr))
				error('SEGDAT:segDat:unequalInputArgumentSizes',...
					'Inputs must have the same size.')
			end%if
			
			obj.x = x;
			obj.y = y;
			obj.s = s;
			obj.k = k;
			obj.phi = phi;
			obj.type = type;
			obj.nbr = nbr;
			
		end%fcn
		
	end%CONSTRUCTOR-methods
	
	
	%%% User-facing methods
	methods
		
		function obj = changeSignOfCurvature(obj)
		% CHANGESIGNOFCURVATURE		Change street segments curvature sign.
		%	OBJ = CHANGESIGNOFCURVATURE(OBJ) changes the curvature sign of
		%	street segment OBJ while maintaining the initial point.
		%	
		%	The other properties are manipulated accordingly.
		
			% get starting point/angle
			P0 = [obj.x(1); obj.y(1)];
			phi0 = obj.phi(1);
			
			% shift to origin and rotate so initial slope is zero
			obj = shift(obj);
			obj = rotate(obj,-phi0);
			
			% invert the curvature
			obj = segDat(...
				+obj.x,...
				-obj.y,...
				+obj.s,...
				-obj.k,...
				-obj.phi,...
				+obj.type,...
				+obj.nbr);
			
			% undo the shift/rotate procedure
			obj = rotate(obj,phi0);
			obj = shift(obj,P0);
			
		end%fcn
		
		
		function obj = plus(obj1,obj2)
		%+ Plus.
		%	OBJ12 = OBJ1 + OBJ2 adds the street segment data OBJ2 with its
		%	starting point to the end point of street segment data OBJ1
		%	resulting in the street segment data OBJ12.
		%	
		%	Note that here plus (+) is a non-commutative operation!
		
		
			obj = segDat(...
				[obj1.x; obj1.x(end) + obj2.x(2:end) - obj2.x(1)],...
				[obj1.y; obj1.y(end) + obj2.y(2:end) - obj2.y(1)],...
				[obj1.s; obj1.s(end) + obj2.s(2:end)],...
				[obj1.k; obj2.k(2:end)],...
				[obj1.phi; obj2.phi(2:end)],...
				[obj1.type; obj2.type(2:end)],...
				[obj1.nbr; obj2.nbr(2:end)+obj1.nbr(end)]);
			
		end%fcn
		
		
		function obj = reverseDirection(obj)
		% REVERSEDIRECTION	Reverse street segment direction.
		%	OBJ = REVERSEDIRECTION(OBJ) reverses the direction of street
		%	segment OBJ so [x(end),y(end)] becomes [x(1),y(1)] and so on.
		%	
		%	The other properties are manipulated accordingly.
			
			
			%%% handle input arguments
			narginchk(1,1);
			
			%%% reverse direction of segment
			obj = segDat(...
				+flipud(obj.x),...
				+flipud(obj.y),...
				+obj.s(end)-flipud(obj.s),... % in case of non-equally distributed x/y
				-obj.k,...
				+flipud(obj.phi) + pi,... 
				+obj.type,...
				+obj.nbr); 
			
		end%fcn
		
		
		function obj = rotate(obj,phi)
		% ROTATE	Rotate street segment.
		%	OBJ = ROTATE(OBJ,PHI) rotates the street segment OBJ by an
		%	angle PHI in radians.
			
			
			%%% handle input arguments
			narginchk(2,2);
			
			if numel(phi) ~= 1 || ~isnumeric(phi)
				error(['Method ROTATE requires a numeric input',...
					' argument with one element.']);
			end%if
			
			
			%%% rotate segment
			% rotation matrix in R^2.
			rotMat = [...
				cos(phi) -sin(phi);...
				sin(phi) +cos(phi)];
			
			xy_new = rotMat*[obj.x';obj.y'];
			
			obj = segDat(...
				xy_new(1,:)',...
				xy_new(2,:)',...
				obj.s,...
				obj.k,...
				obj.phi+phi,...
				obj.type,...
				obj.nbr); 
			
		end%fcn
		
		
		function obj = shift(obj,P)
		% SHIFT		Shift street segment.
		%	OBJ = SHIFT(OBJ,P) shifts the street segment OBJ so that its
		%	starting point [OBJ.x(1) OBJ.y(1)] matches P.
		%	
		%	OBJ = SHIFT(OBJ) applies the default value [0 0] for P.
			
			
			%%% handle input arguments
			narginchk(1,2);
			
			if nargin < 2
				P = [0 0];
			end%if
			
			if numel(P) ~= 2 || ~isnumeric(P)
				error(['Method SHIFT requires a numeric input',...
					' argument with two elements.']);
			end%if
			
			
			%%% shift segment
			obj = segDat(...
				obj.x - obj.x(1) + P(1),...
				obj.y - obj.y(1) + P(2),...
				obj.s,...
				obj.k,...
				obj.phi,...
				obj.type,...
				obj.nbr); 
			
		end%fcn
		
	end%methods
	
	
	%%% User-facing methods (plot related)
	methods
		
		function h = plot(obj,varargin)
		%PLOT	Plot the street segment.
		%	PLOT(OBJ) plots OBJ.y over OBJ.x
		%	
		%	PLOT(OBJ,S) additionally applies the line specification
		%	S.
		%	
		%	H = PLOT(...) returns the handle H to lineseries objects.
		%	
		%	The line specification S is a character string supported by the
		%	standard PLOT command. For example
		%		PLOT(OBJ,'LineWidth',2,'Color',[.6 0 0]) 
		%	will create a plot with a dark red line width of 2 points.
		%	
		%	See also PLOT, segDat/PLOTTANGENT, segDat/PLOTDIFF.
		
		
			% apply plot options if unspecified
			if isempty(varargin)
%				plotOpts = {'o','MarkerSize',2,'MarkerFaceColor','blue'};
				plotOpts = {'-b','LineWidth',2};
			else
				plotOpts = varargin;
			end%if
			
			% check for class
			if ~isa(obj,'segDat')
				error('OBJ has to be of class ''segDat''!')
			end%if
			
			% plot street segment
			h = plot(obj.x,obj.y,plotOpts{:});
			hold on; 
			plot(obj.x(1),obj.y(1),plotOpts{:},'Marker','o');
			hold off
			grid on;
			axis equal;
			title(getLegendCellString(obj));
			ylabel('y [m]');
			xlabel('x [m]');
		
		end%fcn
		
		
		function h = plotdiff(obj,fh)
		%PLOTDIFF	Plot the street segment with specific appearance.
		%	PLOTDIFF(OBJ) plots each street segment type of OBJ using the
		%	according pre-defined type-color.
		%	
		%	PLOTDIFF(OBJ,FH) similar to PLOTDIFF(OBJ) but marker symbols
		%	replace the solid plot line and the marker symbols are switched
		%	periodically at indices specified by the function handle FH.
		%	Some usefull values of FH might be
		%		.) @(OBJ) diff(OBJ.type) 
		%		.) @(OBJ) diff(OBJ.type) | diff(sign(OBJ.k))
		%		.) @(OBJ) diff(OBJ.nbr)
		%	where the first one is used by PLOTDIFF(OBJ). You can also
		%	specify the indices where marker symbols change manually by
		%	using something like
		%		.) @(OBJ) [0 0 1 0 1 0 0 0 1 ...].
		%	Internally FH is evaluated to FH(OBJ) ~= 0.
		%	
		%	H = PLOTDIFF(...) returns a column vector of handles to
		%	lineseries objects, one handle per plotted street segment.
		%	
		%	To modify the plot colors set the hidden property plotColour.
		%	The default colors are blue (straight), green (circle) and red
		%	(clothoid).
		%	
		%	To modify the plot marker set the hidden property plotMarker,
		%	where row one applies to the marker symbol and row two to the
		%	marker size. The default marker symbols are 'o' and 'diamond'
		%	repeated periodically.
		%	
		%	See also segDat/PLOT, segDat/TANGENT.
		
		
			%%% handle input arguments
			if nargin < 2; fh = @(arg) diff(arg.type); end%if
			
			if ~isa(fh,'function_handle')
				error('class')
			end%if
			
			
			%%% index numbering of number of elements of OBJ
			indBase = 1:length(obj.x);
			
			
			%%% get the segment-grouping indices and the number of indices
			flag = fh(obj) ~= 0;
			if length(flag) > length(indBase)
				error(['The length of your evaluated function handle ',...
					'FH exceeds the length of the street segment: ',...
					'%u > %u.'],length(flag),length(indBase));
			end%if
			ind = [0,indBase(flag),indBase(end)];
			nbrInd = length(ind);
			
			
			%%% extend the plot settings to the number of segments to plot
			n = ceil(nbrInd/size(obj.plotMarker,2));
			plotMarker_ = repmat(obj.plotMarker,1,n);
			
			
			%%% plot the segments with according plot options
			h = zeros(nbrInd-1,1);
			for i = 1:nbrInd-1
				
				% plot range
				indRange = ind(i)+1:ind(i+1);
				
				% create the segment data to plot
				sd = segDat(...
					obj.x(indRange),...
					obj.y(indRange),...
					obj.s(indRange),...
					obj.k(indRange),...
					obj.phi(indRange),...
					obj.type(indRange),...
					ones(size(indRange'))...
					);
				
				hold on
				h(i) = plot(sd,...
					'LineStyle','none',...
					'Color',obj.plotColor{obj.type(indRange(1))+1},...
					'Marker',plotMarker_{1,i},...
					'MarkerSize',plotMarker_{2,i});
				hold off
				
			end%for
			
			% override single segment legend with legend of connected OBJ
			title(getLegendCellString(obj));
			
			% unsure about usefulness
			% line style plotting if no additional ....
% 			if nargin < 2
% 				set(h,'LineStyle','-','LineWidth',2,'Marker','none');
% 			end%if
			
		end%fcn
		
		
		function h = plottangent(obj,ind,varargin) 
		%PLOTTANGENT	Plot the street segment and specified tangents.
		%	PLOTTANGENT(OBJ,IND) plots OBJ.y over OBJ.x, the tangents of
		%	the elements of indices IND and highlights the elements IND by
		%	a marker '*'.
		%	
		%	PLOTTANGENT(OBJ,IND,C) additionally applies the color C to the
		%	tangent and the marker.
		%	
		%	PLOTTANGENT(OBJ,IND,C1,...,CN) applies the colors C1,...,CN to
		%	IND(1),...,IND(N). If length(IND) > N, the pattern C1,...,CN
		%	will be extended to match length(IND). If length(IND) < N, the
		%	extra colors will be ignored.
		%	
		%	H = PLOTTANGENT(...) returns a handle H to lineseries objects.
		%	H(1,1) is the street segment-handle, H(i+1,1) and H(i+1,2) the
		%	marker- and the tangent-handle of IND(i) respectively.
		%	
		%	See also segDat/PLOT, segDat/PLOTDIFF.
		
		
			% check dimension of input ind
			if ~isvector(ind)
				error('Input argument IND has to be of size 1xN or Nx1.');
			end%if
			
			% apply plot options if unspecified
			if nargin < 3
				tangentColor = {'r'};
			else
				tangentColor = varargin;
			end%if
			
			% extend the color specification to the length of ind
			if length(ind) > length(tangentColor)
				tangentColor = repmat(tangentColor,1,ceil(length(ind)/length(tangentColor)));
			end
			
			% init handle to lineseries objects
			h = zeros(1+length(ind),2);
			
			% plot the segment and get corresponding axis limtis
			h(1,1) = plot(obj);
			h(1,2) = NaN;
			xLimits = xlim;
			yLimits = ylim;
			
			% check if some indexes are out of range
			if any(ind > length(obj.x))
				warning('MATLAB:lkaSegment:plottangent:indxOutOfRange',...
					['Some indexes exceed the number of elements ',...
					'of segment data which is ',num2str(length(obj.x)),'. ',...
					'Therefore their tangents will not be plotted.'])
				ind(ind > length(obj.x)) = [];
			end%if
			
			% plot the tangents and the corresponding segment-elements
			hold on
			for i = 1:length(ind)
				iind = ind(i);
				
				% marker of tangent point
				h(i+1,1) = plot(obj.x(iind),obj.y(iind),...
					'Marker','*','MarkerSize',7,'Color',tangentColor{i});
				
				% length of tangent
				[r1,r2] = obj.scaleTangentToAxis(xLimits,yLimits,...
					[obj.x(iind) obj.y(iind)],obj.phi(iind));
				
				% start/end point of tangent
				Pstart = [...
					obj.x(iind)+r2*cos(obj.phi(iind));...
					obj.y(iind)+r2*sin(obj.phi(iind))];
				Pstop = [...
					obj.x(iind)+r1*cos(obj.phi(iind));...
					obj.y(iind)+r1*sin(obj.phi(iind))];
				
				% draw the tangent using N arrows sticked together
				N = 25;
				xq = linspace(Pstart(1),Pstop(1),N);
				yq = linspace(Pstart(2),Pstop(2),N);
				% quiver would draw one arrow at every point of xq/yq of
				% length uq/vq; since the last element of xq/yq is the end
				% point, no arrow has to be drawn there
				xq = xq(1:end-1);
				yq = yq(1:end-1);
				uq = ones(size(xq))*(Pstop(1)-Pstart(1))/N;
				vq = ones(size(yq))*(Pstop(2)-Pstart(2))/N;
				scale = 0;
				h(i+1,2) = quiver(xq,yq,uq,vq,scale,'Color',tangentColor{i});
				
			end%for
			hold off
			
			% set the axis limtis corresponding to segment data
			axis([xLimits,yLimits]);
			
		end%fcn
		
	end%methods
	
	
	methods (Access = private)
		
		function cellStr = getLegendCellString(obj)
			
			ts = {'straight','circle','clothoid'};
			if all(obj.type(1) == obj.type)
				gettype = ts{obj.type(1)+1};
			else
				gettype = 'connected';
			end%if
			lengthStr = [sprintf('%.2f',obj.s(end)),' m'];
			pointsStr = sprintf('%.0d',length(obj.s));
			startStr = sprintf('(%g;%g)',obj.x(1),obj.y(1));
			endStr = sprintf('(%g;%g)',obj.x(end),obj.y(end));
			cellStr = {...
				['Street segment of type ''',gettype,''''],...
				['length: ',lengthStr,' using ',pointsStr,' points'],...
				['moves from ',startStr,' \rightarrow ', endStr]};
			
		end%fcn
		
	end%methods
	
	
	%%% SET-Methods
	methods
		
		function obj = set.plotColor(obj,value)
			
			[m,n] = size(value);
			if min(m,n) > 1
				error('Color specification seems to be matrix, but should be vector.')
			end%if
			
			if max(m,n) < 3
				error(['You have to specify at least three colors. ',...
					'One for straight/circle/clohoid segments.'])
			end%if
			
		end%fcn
		
	end%SET-Methods
	
	
	%%% Static-Methods
	methods (Static)
		
		function [r1 r2] = scaleTangentToAxis(xLimits,yLimits,xy,phi)
		%SCALETANGENTTOAXIS		Scale length of tangent to axis limits.
		%   [R1 R2] = SCALETANGENTTOAXIS(XLIMITS,YLIMITS,XY,PHI) calculates
		%   the lengths R1 and R2 for the tangent at point XY with angle
		%   PHI, so that the tangent does not exceed given limits XLIMITS =
		%   [xMin xMax] and YLIMITS = [yMin yMax].
		%	Length R1 counts in the direction of PHI, whereas R2 counts in
		%	the opposite direction.
		%
		%	The intended usage is for plotting something like tangents in
		%	existing plot figures maintaining the current axis limits.
			
			
			% assign inputs to meaningful variables
			xMin = xLimits(1);
			xMax = xLimits(2);
			yMin = yLimits(1);
			yMax = yLimits(2);
			xT = xy(1);
			yT = xy(2);
			
			%{ 
			this approach fails since some angles alphX do not exist if the
			tangents point is located in one of the plot corners.
			
			% map phi to [0,2*pi)
			phi = mod(phi,2*pi);
			
			% get the angles from tangent point to XLIMITS/YLIMITS corners
			alph1 = mod(atan2(yMax-yT,xMax-xT)+2*pi,2*pi);
			alph2 = mod(atan2(yMax-yT,xMin-xT)+2*pi,2*pi);
			alph3 = mod(atan2(yMin-yT,xMin-xT)+2*pi,2*pi);
			alph4 = mod(atan2(yMin-yT,xMax-xT)+2*pi,2*pi);
			
			% calc r1
			if 0 <= phi && phi <= alph1
				r1 = (xMax-xT)/cos(phi);
				
			elseif alph1 < phi && phi <= alph2
				r1 = (yMax-yT)/sin(phi);
				
			elseif alph2 < phi && phi <= alph3
				r1 = (xMin-xT)/cos(phi);
				
			elseif alph3 < phi && phi <= alph4
				r1 = (yMin-yT)/sin(phi);
				
			elseif alph4 < phi && phi <= 2*pi
				r1 = (xMax-xT)/cos(phi);
				
			else
				error('Fatal error at calculating r1!');
			end%if
			
			% calc r2
			phi_ = mod(phi + pi,2*pi);
			if 0 <= phi_ && phi_ <= alph1
				r2 = (xMax-xT)/cos(phi);
				
			elseif alph1 < phi_ && phi_ <= alph2
				r2 = (yMax-yT)/sin(phi);
				
			elseif alph2 < phi_ && phi_ <= alph3
				r2 = (xMin-xT)/cos(phi);
				
			elseif alph3 < phi_ && phi_ <= alph4
				r2 = (yMin-yT)/sin(phi);
				
			elseif alph4 < phi_ && phi_ <= 2*pi
				r2 = (xMax-xT)/cos(phi);
				
			else
				error('Fatal error at calculating r2!');
			end%if
			%}
			
			
			%%% map phi to [-pi,+pi) (plus/minus pi)
			if phi >= pi
				phi_pmPi = phi - 2*pi;
			else
				phi_pmPi = phi;
			end%if
			
			
			%%% reduce problem to angles of [-pi/2,+pi/2]
			if phi_pmPi <= pi/2 && phi_pmPi >= -pi/2
				% keep the original value
				phi_pmPi_ = phi_pmPi;
				isPhiReversed = false;
			else
				% consider the opposite direction of phi which is whithin
				% [-pi/2,+pi/2]
				phi_pmPi_ = mod(phi_pmPi - pi,pi); % ??? consider the sign of x in y for mod(x,y)
				isPhiReversed = true;
			end%if
			
			
			%%% calculate the equation of a line at [xT yT] for phi
			% y(xT) = yT -> k*xT + d = yT -> d = yT - k*xT
			k1 = tan(phi_pmPi_); % tangents slope
			d = yT - k1*xT; % tangents offset
			y_x = @(x) k1*x + d; % tangents equation of a line
			
			
			%%% calculate the plot-range [xL,xR] depending on the sign of
			%%% the tangents slope
			if k1 >= 0
				
				if y_x(xMax) <= yMax
					xR = xMax;
				else
					xR = (yMax-d)/k1;
				end%if
				
				if y_x(xMin) >= yMin
					xL = xMin;
				else
					xL = (yMin-d)/k1;
				end%if
				
			else
				
				if y_x(xMax) >= yMin
					xR = xMax;
				else
					xR = (yMin-d)/k1;
				end%if
				
				if y_x(xMin) <= yMax
					xL = xMin;
				else
					xL = (yMax-d)/k1;
				end%if
				
			end%if
			
			
			%%% calculate the tangent lengths
			r1 = +sqrt((xR-xT)^2 + (y_x(xR)-y_x(xT))^2);
			r2 = -sqrt((xL-xT)^2 + (y_x(xL)-y_x(xT))^2);
			
			% fix tangent lengths if the slope is +-inf
			if phi_pmPi_ == +pi/2
				r1 = yMax - yT;
				r2 = yMin - yT;
			elseif phi_pmPi == -pi/2
				r1 = -(yMin - yT);
				r2 = -(yMax - yT);
			end%if
			
			
			%%% switch R1/R2 and their signs if the opposite problem was
			%%% considered
			if isPhiReversed
				dummy = r1;
				r1 = -r2;
				r2 = -dummy;
			end%if
			
		end%fcn
		
	end%methods
	
	
	%%% Test-Methods
	methods (Static, Hidden)
		
		function test_changeSignOfCurvature(obj)
			
			if nargin < 1;
				b = lkaSegmentCircle([],5/3*pi,4/2*pi,50);
				b = shift(b,[20,10]);
				obj = b.segmentData;
			end%if
			
			fig = figure;
			
			ind = [1,length(obj.x)];
			h = plottangent(obj,ind);
			set(h(1,1),'Color','b','LineWidth',1);
			set(h(2:end,1),'Color','c','Marker','o','MarkerFaceColor','c');
			set(h(2:end,2),'Color','c');
			
			obj_ = changeSignOfCurvature(obj);
			hold on
			h_ = plottangent(obj_,ind);
			set(h_(1,1),'Color','r','LineWidth',1);
			set(h_(2:end,1),'Color','m','Marker','h','MarkerFaceColor','m');
			set(h_(2:end,2),'Color','m');
			
			axis auto
			legend([h(1,1),h_(1,1)],'original','curvature*(-1)','location','Best')
			pause
			close(fig)
			
		end%fcn
		
		
		function test_rotate(P0,phi)
			
			if nargin < 2; phi = pi/2; end%if
			if nargin < 1; P0 = [20 0]; end%if
			
			fig = figure;
			b = lkaSegmentCircle([],3/2*pi,4/2*pi,50);
			
			sd_b = b.segmentData;
			
			sd_b = shift(sd_b,P0);
			
			plotdiff(sd_b);
			
			plotdiff(rotate(sd_b,phi));
			
			pause
			close(fig)
			
		end%fcn
		
		
		function test_reverseDirection(obj)
			
			fig = figure;
			
			obj_ = reverseDirection(obj);
			
			ind = 1:10:21;
			h = plottangent(obj,ind,'b');
			
			hold on
			h = plottangent(obj_,ind,'r');
			set(h(1,1),'LineStyle',':','Color','r');
			
			pause
			close(fig)
			
		end%fcn
		
		
		function indFailed = test_scaleTangentToAxis(nbr)
			
			
			% define some test cases and their desired results
			tc = [...
				%solution	xLimits		yLimits		point		angle
				%
				% tangent point at lower left corner
				{[+10 0],	[0 10],		[0 10],		[0 0],		0};...
				{[+10 0],	[0 10],		[0 10],		[0 0],		pi/2};...
				{[0 -10],	[0 10],		[0 10],		[0 0],		pi};...
				{[0 -10],	[0 10],		[0 10],		[0 0],		3*pi/2};...
				%
				% tangent point at lower right corner
				{[0 -10],	[0 10],		[0 10],		[10 0],		0};...
				{[+10 0],	[0 10],		[0 10],		[10 0],		pi/2};...
				{[+10 0],	[0 10],		[0 10],		[10 0],		pi};...
				{[0 -10],	[0 10],		[0 10],		[10 0],		3*pi/2};...
				%
				% tangent point at upper right corner
				{[0 -10],	[0 10],		[0 10],		[10 10],	0};...
				{[0 -10],	[0 10],		[0 10],		[10 10],	pi/2};...
				{[+10 0],	[0 10],		[0 10],		[10 10],	pi};...
				{[+10 0],	[0 10],		[0 10],		[10 10],	3*pi/2};...
				%
				% tangent point at upper left corner
				{[+10 0],	[0 10],		[0 10],		[0 10],		0};...
				{[0 -10],	[0 10],		[0 10],		[0 10],		pi/2};...
				{[0 -10],	[0 10],		[0 10],		[0 10],		pi};...
				{[+10 0],	[0 10],		[0 10],		[0 10],		3*pi/2};...
				%
				% tangent point at upper left corner
				% k = tan(phi)
				% d = k*x - y(x) = 5 - 5*tan(phi)
				% r1 = sqrt(5^2 + (5-d)^2)
				% r2 = -r1
				{[sqrt(25+25*tan(pi/8)^2) -sqrt(25+25*tan(pi/8)^2)],...
							[0 10],		[0 10],		[5 5],		pi/8};
			];
		
% 			cell2struct(tc,{'sol','xLimits','yLimits','point','angle'},2);
			
			% handle input arguments
			if nargin < 1
				nbr = 1:size(tc,1);
			end%if
			
			% select test cases to perform
			tc = tc(nbr,:);
			nbrOfTestCases = size(tc,1);
			
			% pre-allocation
			isTestSuccesfull = false(nbrOfTestCases,1);
			
			% run test cases
			for i = 1:nbrOfTestCases
				
				[r1 r2] = segDat.scaleTangentToAxis(...
					tc{i,2},...
					tc{i,3},...
					tc{i,4},...
					tc{i,5});
				
				isTestSuccesfull(i) = all(tc{i,1} == [r1 r2]);
				
			end%for
			
			% evaluation and command line output
			nbrSuccesfull = sum(isTestSuccesfull);
			nbrFailed = sum(~isTestSuccesfull);
			fprintf('   Test results: \n')
			fprintf('     succesfull: %3d of %d.\n',nbrSuccesfull,nbrOfTestCases)
			fprintf('     failed:     %3d of %d.\n',nbrFailed,nbrOfTestCases)
			
			ind = 1:length(tc);
			indFailed = ind(~isTestSuccesfull);
			formatString = repmat('%d, ',1,nbrFailed);
			formatString = [formatString(1:end-2),'.'];
			if ~isempty(indFailed)
				fprintf(['     test cases failed: ',formatString,'\n'],...
					indFailed);
			end%if
			fprintf('\n');
			
		end%fcn
		
		
		function test_shift()
			
			fig = figure;
			a = lkaSegmentStraight([],100,pi/8);
			b = lkaSegmentCircle([],3/2*pi,4/2*pi,50);
			c = lkaSegmentClothoid([],1/100,0,0,100);
			
			sd_a = a.segmentData;
			sd_b = b.segmentData;
			sd_c = c.segmentData;
			
			sd_a = shift(sd_a,[10 30]);
			sd_b = shift(sd_b,[20 20]);
			sd_c = shift(sd_c,[30 10]);
			
			plotdiff(sd_a);
			plotdiff(sd_b);
			plotdiff(sd_c);
			
			plotdiff(shift(sd_a));
			plotdiff(shift(sd_b));
			plotdiff(shift(sd_c));
			
			pause
			close(fig)
			
		end%fcn
		
	end%methods
	
end%class