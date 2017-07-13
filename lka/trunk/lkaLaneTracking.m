function [out,latOff_LAD,angDev_LAD,curvat_LAD,isValid_LAD] = ...
		lkaLaneTracking(obj,xyCG_global,yawAngle_global,LAD,mode)
% LANETRACKING  Calculate lane tracking pose.
%	
%	[~,LATOFF_LAD,ANGDEV_LAD,CURVAT_LAD,ISVALID_LAD] =
%	LANETRACKING(OBJ,XYCG_GLOBAL,YAWANGLE_GLOBAL,LAD) calculates the
%	lateral offset LATOFF_LAD, the angular deviation ANGDEV_LAD and the
%	curvature CURVAT_LAD at the look-ahead distances LAD in front of the
%	vehicles center of gravity XYCG_GLOBAL along its longitudinal axis
%	oriented with the angle YAWANGLE_GLOBAL and with respect to the street
%	segment OBJ. The logical vector ISVALID_LAD indicates if the
%	corresponding entry of all other output arguments is valid (true) or
%	not (false).
%	
%	OUT = LANETRACKING(...) is the syntax to be used from simulink, where
%	OUT = [LAD;ISVALID_LAD;LATOFF_LAD;ANGDEV_LAD;CURVAT_LAD].
%	
%	The number of columns of all output arguments matches the length of
%	input argument LAD, moreover the i-th column of any output argument
%	corresponds to the i-th element of look-ahead distance LAD.

% Subject: lka
% $Author$
% $LastChangedDate$
% $Revision$

	% Methode: Koordinatentransformation
	if nargin <5
		mode = 0;
	end%if
	
	% ensure row format
	LAD = LAD(:)';
	
	%%% shift origin to vehicles CG/rotate so vehicle is oriented along
	%%% x-axis
	% CG (transformed)
	xyCG_T = [0;0]; % xyCG_global - xyCG_global
	
	% Sollbahn (transformed)
% 	obj_T = shift(obj, [obj.x(1);obj.y(1)]-xyCG_global);
	P = [obj.x(1);obj.y(1)] - xyCG_global;
	obj_T.x = obj.x - obj.x(1) + P(1);
	obj_T.y	= obj.y - obj.y(1) + P(2);
	obj_T.s	= obj.s;
	obj_T.k = obj.k;
	obj_T.phi = obj.phi;
	obj_T.type = obj.type;
	
% 	obj_T = rotate(obj_T, -yawAngle_global);
	rotMat = [...
		cos(-yawAngle_global) -sin(-yawAngle_global);...
		sin(-yawAngle_global) +cos(-yawAngle_global)];
	
	% much faster than rotMat*[obj.x;obj.y]
	xy_new = [obj_T.x' obj_T.y']*rotMat';
	obj_T.x = xy_new(:,1)';
	obj_T.y = xy_new(:,2)';
	obj_T.s = obj_T.s;
	obj_T.k = obj_T.k;
	obj_T.phi = obj_T.phi-yawAngle_global;
	obj_T.type = obj_T.type;
	
	% Koordinaten des Punkts bei look-ahead distance (transformed)
	xyLAD_T = [xyCG_T(1) + LAD; xyCG_T(2)+zeros(size(LAD))];
	
	
% 	%%% shift origin to point at LAD
% 	xyCG_T	= xyCG_T - [lad;0];
% 	obj_T	= shift(obj_T,[obj_T.x(1);obj_T.y(1)] - [lad;0]);
%	xyLAD_T	= xyLAD_T - [lad;0];
	
	
	%%% get indices of potential elements of desired path
	% maximaler Abstand zwischen x-Werten der Sollbahn
	deltaX_max = max(abs(diff(obj_T.x)));
	
	% Depending on the shape of the desired path (e.g. closed path), there
	% might be multiple elements of the desired path whose x-coordinates
	% are within the range of LAD +- DELTAX_MAX/2.
	% 
	% Get logical indices where OBJ_T.X is in the range of LADs
	% x-coordinate
	%  rows: LAD 
	%  columns: street segments x-coordinate
	logIndx1 = bsxfun(@le,xyLAD_T(1,:)' - deltaX_max/2,obj_T.x);
	logIndx2 = bsxfun(@ge,xyLAD_T(1,:)' + deltaX_max/2,obj_T.x);
	logIndx = logIndx1 & logIndx2;
	
	% error if no element of LOGINDX is true
	if ~any(any(logIndx))
		plotLaneTracking(obj,xyCG_global,yawAngle_global,LAD,[],obj_T,xyCG_T);
		error('segDat:laneTracking',...
			['Keine Elemente der Sollbahn im Bereich der',... 
			' aktuellen Fahrzeugposition gefunden'])
	end%if
	
	% Numerical indices according to LOGINDX: NUMINDROW points to the
	% according LAD, NUMINDCOL points to the according x-coordinate (see
	% LOGINDX calculation above)
	[numIndRow,numIndCol] = find(logIndx);
	
	if mode
		%%% reduce number of potential elements
		% since there might be multiple elements within LAD +-
		% DELTAX_MAX/2, get the closest in x-direction per LAD AND per
		% contigous indices
		
		% split NUMINDCOL into segments of contigous indices and LAD values
		splitByLAD = find(diff(numIndRow)>0);
		splitByInd = find(diff(numIndCol)>1);	
		splitInd = [0;unique([splitByLAD;splitByInd]);length(numIndRow)];
		
		% get index of element closest to LAD, do it for each index-segment
		% defined by SPLITIND
		numIndRow_red = zeros(length(splitInd)-1,1);
		numIndCol_red = zeros(length(splitInd)-1,1);
		for i = 1:length(splitInd)-1
			numIndRow_i = numIndRow(splitInd(i)+1:splitInd(i+1));
			numIndCol_i = numIndCol(splitInd(i)+1:splitInd(i+1));
			
			% since we grouped by contigous indices and LAD, NUMINDWOR_I
			% must not contain different values
			if ~all(numIndRow_i(1) == numIndRow_i)
				error('This should not have happened!');
			end%if
			
			%
			if length(numIndRow_i) > 2
				[~,winInd] = min(abs(obj_T.x(numIndCol_i) - xyLAD_T(1,numIndRow_i(1))));
			else
				winInd = 1;
			end%if
			
			numIndRow_red(i) = numIndRow_i(winInd);
			numIndCol_red(i) = numIndCol_i(winInd);
		end%for
		numIndRow = numIndRow_red;
		numIndCol = numIndCol_red;
	end%if
	
% 	plotLaneTracking(obj,xyCG_global,yawAngle_global,LAD,numIndCol,obj_T,xyCG_T);
	
	
	%%% get lateral offset for potential elements
	% interpolation index-range: m>1 mainly increase the calculation time,
	% lateral offset and angular deviation are rarely affected.
	m = 1;
	[indl,indu] = interpIndexRange(numIndCol,[1,length(obj.x)],m);
	
	% preallocation of for-loop variable
	lanePose_LAD_candidates = zeros(3,length(numIndCol));
	try
		% Inter- bzw. Extrapoliere y-Werte der transf. Solltrajektorie
		for i = 1:length(numIndCol)
			% same result like interp1(..,'spline') but faster
			lanePose_LAD_candidates(:,i) = spline(...
				obj_T.x(indl(i):indu(i)),...
				[obj_T.y(indl(i):indu(i));...
				 obj_T.phi(indl(i):indu(i));...
				 obj_T.k(indl(i):indu(i))],...
				xyLAD_T(1,numIndRow(i)));
		end%for
	catch exception
		plotLaneTracking(obj,xyCG_global,yawAngle_global,LAD,numIndCol,obj_T,xyCG_T);
		error(exception.message);
	end%try
	latOff_LAD_candidates = lanePose_LAD_candidates(1,:);
	angDev_LAD_candidates = lanePose_LAD_candidates(2,:);
	curvat_LAD_candidates = lanePose_LAD_candidates(3,:);
	
	
	%%% select one of multiple lateral offsets per LAD-value
	% get most possible (by means of smallest) value per LAD-value
	latOff_LAD	= zeros(size(LAD));
	angDev_LAD	= zeros(size(LAD));
	curvat_LAD	= zeros(size(LAD));
	isValid_LAD = false(size(LAD));
	for i = 1:length(latOff_LAD_candidates)
		latOff_underTest = latOff_LAD_candidates(i);
		angDev_underTest = angDev_LAD_candidates(i);
		curvat_underTest = curvat_LAD_candidates(i);
		LAD_group = numIndRow(i);
		
		if abs(latOff_underTest) < abs(latOff_LAD(LAD_group)) || ...
				~isValid_LAD(LAD_group)
			latOff_LAD(LAD_group)	= latOff_underTest;
			angDev_LAD(LAD_group)	= angDev_underTest;
			curvat_LAD(LAD_group)	= curvat_underTest;
			isValid_LAD(LAD_group)	= true;
		end%if
		
	end%for
	
	% set multiple occurences of same LAD values to invalid
	[LAD_unique,ia] = unique(LAD,'stable');
	if numel(LAD_unique) ~= numel(LAD)
		% duplicate indices
		duplicate_ind = setdiff(1:numel(LAD),ia);
		
		% set to false
		isValid_LAD(duplicate_ind) = false;
	end%if
	
	% collect all output arguments (to be used in simulink)
	out = [LAD;isValid_LAD;latOff_LAD;angDev_LAD;curvat_LAD];
	
end%fcn


function [indl,indu,ind] = interpIndexRange(ind,indMinMax,m)
%INTERPINDEXRANGE	Index range for interpolation.
%	[INDL,INDU] = INTERPINDEXRANGE(IND,INDMINMAX,M) returns the lower and
%	upper indices INDL and INDU within index boundaries [INDMINMAX(1)
%	INDMINMAX(2)] using an index difference M for givenen indices IND.
%	
%	Indices INDL/INDU define range of interpolation, basically IND-M and
%	IND+M but ensure INDL>INDMINMAX(1) and INDU<INDMINMAX(2).


	%%% handle input arguments
	% input IND
	if  ~isvector(ind) % check size
		error('Size of IND must be vector!');
	else
		ind = ind(:); % ensure column orientation
	end%if

	% input INDMINMAX
	if any(size(indMinMax) ~= [1 2]) % check size
		error('Size of INDMINMAX must be 1-by-2!');
	else
		indMin = indMinMax(1);
		indMax = indMinMax(2);
	end%if

	% input M
	if nargin < 3 % apply default value if undefined
		m = 1;
	elseif ~isscalar(m) % check isze
		error('Size of M must be scalar!');
	elseif m < 1 % check value
		error('Value of m must be > 0!');
	end%if


	%%% check input arguments plausibility
	if any(ind<indMin | ind>indMax)
		error('laneTracking:interpIndexRange:IndexOutOfBounds',...
			'One or more indexes out of range [%i,%i].',indMin,indMax);
	end%if
	if indMax-indMin < 2*m
		error('laneTracking:interpIndexRange:IntervalExtOutOfRange',...
			'Index difference M=%i does not fit into given range [%i,%i].',...
			2*m,indMin,indMax);
	end%if


	%%% calculate interpolating indexes
	% lower/upper interpolating-index unbounded
	indl = ind-m;
	indu = ind+m;

	% bounded to INDMINMAX
	indLU = bsxfun(@plus,[indl,indu],...
		max(zeros(size(ind)),indMinMax(:,1)-indl) + ...
		min(zeros(size(ind)),indMinMax(:,2)-indu));
	indl = indLU(:,1);
	indu = indLU(:,2);

end%fcn
