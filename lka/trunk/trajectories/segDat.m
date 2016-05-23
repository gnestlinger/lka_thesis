classdef segDat
%SEGDAT 	Street segment data.
%	
%	OBJ = SEGDAT(X,Y,S,K,PHI,TYPE,NBR) creates the street segment object
%	OBJ of class SEGDAT representing the geometric data of a street segment
%	using the properties listed below.
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
%	A curvature K > 0 (K < 0) indicates a left (right) turn by moving from
%	[X(1),Y(1)] -> [X(end),Y(end)].
%	
%	TYPE = 0 indicates a straight, 1 a circular and 2 a clothoidal street
%	segment.
%	
%	NBR is usually 1, but becomes interesting when some street segments are
%	connected together. Then NBR represents the number of the street
%	segment in the path of connected segments.
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
		% x-coordinate.
		x
		
		% y-coordinate.
		y
		
		% Covered distance.
		s
		
		% Curvature.
		k
		
		% Tangent angle.
		phi
		
		% Street segment type.
		type
		
		% Street segment number of connected segments.
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
		
		
			% all inputs have to be non-empty column vectors, just check x
			% and check the size of all others against x
			[m,n] = size(x);
			if (m<1) || (n>1)
				error('Inputs have to be non-empty column vectors')
			end%if
			if ~isequal(size(x),size(y),size(s),size(k),size(phi),size(type),size(nbr))
				error('Inputs must have the same size.')
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
		
		%%% connect street segment data
		function obj = plus(obj1,obj2)
		%+ Plus.
		%	SEGDAT = SEGDAT1 + SEGDAT2 adds the street segment data SEGDAT2
		%	with its starting point to the end point of street segment data
		%	SEGDAT1 resulting in the street segment data SEGDAT.
		%	
		%	Note that here + is a non-commutative operation!
		
		
			obj = segDat(...
				[obj1.x; obj1.x(end) + obj2.x(2:end) - obj2.x(1)],...
				[obj1.y; obj1.y(end) + obj2.y(2:end) - obj2.y(1)],...
				[obj1.s; obj1.s(end) + obj2.s(2:end)],...
				[obj1.k; obj2.k(2:end)],...
				[obj1.phi; obj2.phi(2:end)],...
				[obj1.type; obj2.type(2:end)],...
				[obj1.nbr; obj2.nbr(2:end)+obj1.nbr(end)]);
			
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
		
		
			
		%%% plot of street segment
		function h = plot(obj,varargin)
		%PLOT	Plots the street segment.
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
			grid on;
			axis equal;
			title(getLegendCellString(obj));
			ylabel('y [m]');
			xlabel('x [m]');
		
		end%fcn
		
		
		%%% tangent plot of street segment
		function h = plottangent(obj,ind,varargin) 
		%PLOTTANGENT	Plots the street segment and specified tangents.
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
			
			% plot the tangets and the corresponding segment-elements
			hold on
			for i = 1:length(ind)
				iind = ind(i);
				h(i+1,1) = plot(obj.x(iind),obj.y(iind),...
					'Marker','*','MarkerSize',7,'Color',tangentColor{i});
% 				r1 = -1e3; r2 = +1e3;
				[r1,r2] = obj.scaleTangentToAxis(xLimits,yLimits,...
					[obj.x(iind) obj.y(iind)],obj.phi(iind));
% 				h(i+1,2) = plot(...
% 					obj.x(iind)+[r1*cos(obj.phi(iind));r2*cos(obj.phi(iind))],...
% 					obj.y(iind)+[r1*sin(obj.phi(iind));r2*sin(obj.phi(iind))],...
% 					'Color',tangentColor{i});
				
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
				uq = ones(size(xq))*(Pstop(1)-Pstart(1))/N;
				vq = ones(size(yq))*(Pstop(2)-Pstart(2))/N;
				scale = 0;
				h(i+1,2) = quiver(xq,yq,uq,vq,scale,'Color',tangentColor{i});
				
			end%for
			hold off
			
			% set the axis limtis corresponding to segment data
			axis([xLimits,yLimits]);
			
		end%fcn
		
		
		%%% plot of street segment using appearance variations
		function h = plotdiff(obj,fh)
		%PLOTDIFF	Plots the street segment with specific appearance.
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
			
			% unsure about usefulness
			% line style plotting if no additional ....
			if nargin < 2
				set(h,'LineStyle','-','LineWidth',2,'Marker','none');
			end%if
			
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
			
		end%fcn
		
	end%methods
	
	
	%%% Test-Methods
	methods (Static, Hidden)
		
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
		
	end%methods
	
end%class
