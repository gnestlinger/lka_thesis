classdef lkaSegment
%LKASEGMENT 	Superclass to create a path of street segments.
%	
%	This class is not user-facing and cannot be instantiated!
%	User-facing subclasses include:
%	  * Straight segments:	LKASEGMENTSTRAIGHT
%	  * Circular segments:	LKASEGMENTCIRCLE
%	  * Clothoid segments:	LKASEGMENTCLOTHOID
%	  * Connected segments: LKASEGMENTCONNECT
%	
%	LKASEGMENT Properties:
%	 segmentType - Name of the segment type.
%	 deltaSet	 - Desired distance between nearby points of segment.
%	 deltaAct	 - Actual distance between nearby points of segment.
%	 nbrOfPoints - Number of points representing the segment.
%	 xyStart	 - Coordinate of segment data starting point.
%	 xyStop		 - Coordinate of segment data end point.
%	 segmentData - The segment data.
%	 length		 - Arc length of the segment.
%	
%	LKASEGMENT Methods:
%	 - DESIGN
%	 plus		 - Connect segments.
%	 resample	 - Apply a new value for DELTASET.
%	 rotate		 - Rotate segment.
%	 shift		 - Shift segment.
% 
%	 - ANALYSIS
%	 plot		 - Plot segment.
%	 plotdiff	 - Plot segment using segment-type specific plot styles.
%	 plottangent - Plot segment and tangents of specified segment points.
%	
%	--- (just used by subclasses) ---------------------------------------
%	OBJ = LKASEGMENT(SEGMENTTYPE,DELTASET,XYSTART) sets the segment type
%	SEGMENTTYPE, the desired spacing DELTASET between consecutive points of
%	the segment and the coordinate XYSTART of segment's starting point.
%	
%	OBJ = LKASEGMENT(SEGMENTTYPE,[],XYSTART) applies the default value
%	for DELTASET.
%	----------------------------------------------------------------------
%	
%	See also LKASEGMENTSTRAIGHT, LKASEGMENTCIRCLE, LKASEGMENTCLOTHOID,
%	LKASEGMENTCONNECT, SEGDAT.
% 

% Subject: lka
% Author: $Author$
% Date: $LastChangedDate$
% Revision: $Revision$


% DEVELOPMENT NOTES:
%	(1) get rid of property DESIGNPROPERTIES? Keep it, is used by
%	errorMsg_SetDependent.
%	
%	(2) Implement method GETNBROFPINTS_ABSTRACT in supercalls LKASEGMENT
%	and overload it for sublass LKASEGMENTCONNECT?
%	
%	(3) Check if properties ROTMAT* are used, ROTMAT() seems to fail. FIXED
%	in r228 @ 2018-04-10.
%	
%	(4) Fix value of property DELTAACT for class LKASEGMENTCONNECT.

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% PROPERTIES %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

	properties (Constant, Hidden, GetAccess = protected)
        
        % rotMatX - Rotation matrix x-component.
        %   The x-component of a vector p = [x;y] is rotated by an angle
        %   PHI counter-clockwise by rotMatX(PHI)*p.
        rotMatX = @(phi) [cos(phi) -sin(phi)];
        
        % rotMatY - Rotation matrix y-component.
        %   The y-component of a vector p = [x;y] is rotated by an angle
        %   PHI counter-clockwise by rotMatY(PHI)*p.
        rotMatY = @(phi) [sin(phi) +cos(phi)];
        
		% rotMat - Rotation matrix in R^2.
		%	A vector p = [x;y] is rotated by an angle PHI counter-clockwise
		%	by rotMat(PHI)*p.
		rotMat = @(phi) [cos(phi) -sin(phi); sin(phi) +cos(phi)];
		
	end%
	
	
    properties (Constant, Hidden, Abstract)
        
        % designProperties - User adjustable properties for segment design.
        %   Cell of strings defining user-adjustable properties used for
        %   segment design.
        designProperties;
        
    end%
	
	
    properties (SetAccess = private)
        
        % segmentType - The segment type.
        %   String indicating the type of the segment.
        segmentType(1,1) string; %segment info data
		
        % deltaSet - Desired distance between two consecutive points [m].
		%	The desired distance DELTASET between two consecutive points
		%	can not always be fullfilled exactly (depending on the other
		%	street segment design parameters) and is therefore treated as
		%	an upper bound.
		%	
		%	The default value is DELTASET = 1.
        deltaSet(1,:) double {mustBeFinite,mustBePositive} = 1; %segment design data
        
    end%
	
	
    properties (Dependent, SetAccess = private)
        
        % deltaAct - Actual distance between two nearby points [m].
		%	Currently used distance between two nearby segment points. In
		%	general, this value differs from the desired distance DELTASET.
        deltaAct(1,:) double; %segment info data
        
        % nbrOfPoints - Number of segment-points [-].
		%	This is the minimum and currently used number of segment points
		%	required, to fullfill DELTAACT < DELTASET.
        nbrOfPoints(1,:) uint64; %segment info data
        
	end%
	
	
    properties (SetAccess = private)
        
        % xyStart - Starting point in x/y-plane [m].
        xyStart(1,2) double {mustBeFinite}; %segment design data
        
    end%
	
	
    properties (Dependent, SetAccess = private)
        
        % xyStop - Endpoint in x/y-plane [m].  
        xyStop(1,2) double; %segment info data 
        
        % segmentData - Object of street segment data. 
        segmentData(1,1);
        
    end%
	
	
    properties (Abstract, SetAccess = private)
        
        % length - Arc length of the segment [m].
		%	Depending on the implementation of the subclass, this property
		%	is either just for info purpose or can be set by the user.
        length(1,1) double; %segment design or info data
        
    end%
    
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% METHODS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %%% CONSTRUCTOR
    methods
        
        function obj = lkaSegment(segmentType,deltaSet,xyStart)
            
            obj.segmentType = segmentType;
            
            % set only if specified, otherwise default values are used
            if ~isempty(deltaSet); obj.deltaSet = deltaSet; end%if
            
            obj.xyStart = xyStart;
             
        end%Constructor
        
    end%CONSTRUCTOR-methods
    
    
    %%% User-facing methods
	methods
		
        function obj = plus(obj1,obj2)
        %+ Plus.
        %   OBJ1 + OBJ2 works on the property 'segmentData' and connects
        %   the starting point of OBJ2 to the end point of OBJ1. The result
        %   is of subclass lkaSegmentConnect.
		%	
		%	Note that here + is a non-commutative operation!
        %
        %   See also LKASEGMENTCONNECT.
        
        
            obj_segDat = plus(obj1.segmentData,obj2.segmentData);
            obj = lkaSegmentConnect(obj_segDat);
                
        end%fcn
        
		
		function obj = resample(obj,deltaNew)
        %RESAMPLE   Apply a desired distance between points.
        %   OBJ = RESAMPLE(OBJ,DELTANEW) sets the maximum distance between
        %   two nearby points of street segment OBJ to DELTANEW.
        %
        
            obj.deltaSet = deltaNew;
            
        end%fcn
		
		
		function obj = rotate(obj,phi)
		%ROTATE		Rotate the street segment.
		%	OBJ = ROTATE(OBJ,PHI) rotates the street segment OBJ by an
		%	angle PHI counter-clockwise.
		%
		
			% call abstract method
			obj	= rotate_abstract(obj,phi);
			
		end%fcn
		
		
		function obj = shift(obj,point)
        %SHIFT  Shift the street segment.
        %   OBJ = SHIFT(OBJ,POINT) shifts the starting point of street
        %   segment OBJ to coordinate POINT.
        %
        %   OBJ = SHIFT(OBJ) shifts to the origin [0,0].
        %
        
            if nargin < 2
                point = [0,0];
            end%if
            
            obj.xyStart = point;
            
        end%fcn
		
	end%methods
	
	
	
	%%% User-facing methods (plot related)
	methods
		
        function h = plot(obj,varargin)
        %PLOT   Plot the street segment.
        %
        %   For the documentation see class SEGDAT.
        %
        %   See also SEGDAT/PLOT.
        
        
            h = plot(obj.segmentData,varargin{:});
            
%             % format specification for points
%             formatSpec = '(%g;%g)';
%             
%             % add some info to the plot title
%             title({...
%                 ['Segment type: ',obj.segmentType],...
%                 ['Segment length: ',num2str(obj.length,'%g'),' m'],...
%                 ['Segment points: ',num2str(obj.nbrOfPoints)],...
%                 ['Starting/end point: ',...
%                 sprintf(formatSpec,obj.xyStart(1),obj.xyStart(2)),' / ',...
%                 sprintf(formatSpec,obj.xyStop(1),obj.xyStop(2))],...
%                 });
        
        end%fcn
		
		
        function h = plotdiff(obj,fh)
        %PLOTDIFF   Plots the street segment with specific appearance.
        % 
        %   For the documentation see class SEGDAT.
        %
        %   See also SEGDAT/PLOTDIFF.
        
        
            %%% handle input arguments
			% call SEGDAT-class method
			if nargin < 2
				h = plotdiff(obj.segmentData);
			else
				h = plotdiff(obj.segmentData,fh);
			end%if
        
        end%fcn
        
		
		function h = plottangent(obj,ind,varargin)
        %PLOTTANGENT    Plot the street segment and specified tangents.
        %
        %   For the documentation see class SEGDAT.
        %
        %   See also SEGDAT/PLOTTANGENT.
        
        
            % apply plot options if unspecified
            if nargin < 3
                plotColor = {'r'};
            else
                plotColor = varargin;
            end%if
            
            h = plottangent(obj.segmentData,ind,plotColor{:});
            
        end%fcn
		
    end%methods
    
    
    
    %%% GET-Methods
    methods
        
		function value = get.deltaAct(obj)
			
			value = get_deltaAct_(obj);
			
			if any((value < 0) | (value > obj.deltaSet))
				warning('MATLAB:lkaSegment:get:deltaAct',...
					'Actual delta has an invalid value!!!')
			end%if
			
		end%fcn
        
        
        function value = get.nbrOfPoints(obj)
            
            % call abstract method
            value = uint64(getNbrOfPoints_abstract(obj));
            
        end%fcn
        
        
        function value = get.xyStop(obj)
            
            % call abstract method to get xyStop dependent on subclass
            value = getEndPoint_abstract(obj);
            
        end%fcn
        
        
        function segdat = get.segmentData(obj)
            
			% show calculation notes only for other segments than
			% 'connected'
			if strcmp(obj.segmentType,'clothoid')
				showCalcNote = true;
			else
				showCalcNote = false;
			end%if
			
			if showCalcNote
				disp(['*** Calculating segment ''',obj.segmentType,''''])
			end%if
			
            % call abstract method
            segdat = getSegmentData_abstract(obj);
			
			if showCalcNote
				disp('*** Done :)')
				fprintf('\n')
			end%if
            
        end%fcn
        
    end%GET-Methods
    
    
    
    %%% SET-Methods
    methods
        
        function obj = set.deltaSet(obj,value)
            
            % property name
            pN = 'deltaSet';
            
            % dimension check
            obj.check_dimension_isSVM(pN,value,{'scalar'});
            
            % check data of value
            if ~isreal(value)
                obj.errorMsg_xyz('data set',pN,'real valued','complex valued');
            end%if
            
            if ~isnumeric(value)
                obj.errorMsg_xyz('data type',pN,'numeric','non-numeric');
            end%if
            
            if (value <= 0) % value range check
                obj.errorMsg_xyz('value',pN,'> 0','<= 0');
            end%if
            
            % set the property value
            obj.deltaSet = value;
            
        end%fcn
        
        
        function obj = set.xyStart(obj,value)
            
            % property name
            pN = 'xyStart';
            
            % dimension check
            obj.check_dimension_isSVM(pN,value,{'vector'});
            if numel(value) > 2 % < 2 alredy checked
                obj.errorMsg_xyz('dimension',pN,'two-dimensional','N-dimensional (N>2)');
            end%if
            
            % check data of value
            if ~isreal(value)
                obj.errorMsg_xyz('data set',pN,'real valued','complex valued');
            end%if
            
            if ~isnumeric(value)
                obj.errorMsg_xyz('data type',pN,'numeric','non-numeric');
            end%if
            
            % set the property value
            % ensure row-vector (so the values will be displayed)
            obj.xyStart = obj.col(value)';
            
        end%fcn
        
    end%SET-Methods
    
    
    
    %%% Abstract-Methods
    methods (Abstract, Access = protected)
        % abstract methods to be implemented by subclasses
		
		obj		= rotate_abstract(obj,phi);
		
        nbr		= getNbrOfPoints_abstract(obj);
        endP	= getEndPoint_abstract(obj)      
        segdat	= getSegmentData_abstract(obj)
        
    end%Abstract-Methods
    
	
    methods (Access = protected)
		
		% must not be private for overloading in subclasses
		function value = get_deltaAct_(obj)
			value = obj.length/(double(obj.nbrOfPoints)-1);
		end%fcn
		
	end
    
	
    methods (Access = protected)
        
        %%% error message
        function errorMsg_SetDependent(obj,prop)
            
            value = obj.(prop);
            
            switch class(value)
                case 'double'
                    formatString = ['%s',repmat('%.3E\t',1,length(value)),'\n'];
                    
                case 'char'
                    formatString = '%s\n';
                    
                otherwise
                    
            end
            
            % print current value of prop
            fprintf(formatString,['Property ',prop,' is: '],value);
            
            % the number of adjustable design properties
            nbrOfDesignProps = length(obj.designProperties);
            
            % create messages
            msg1 = ['You cannot set ''',prop,''' explicitly. '];
            msg2 = 'This object owns the following adjustable properties: ';
            msg3 = ['Try setting one or more of these to get the desired ''',prop,'''.'];
            
            if nbrOfDesignProps > 1
                error([msg1,msg2,...
                    repmat('%s, ',1,nbrOfDesignProps-1),'%s. ',...
                    msg3],...
                    obj.designProperties{:});
            elseif nbrOfDesignProps > 0
                error([msg1,msg2,...
                    '%s. ',...
                    msg3],...
                    obj.designProperties{:});
            else
                error([msg1,...
                    'There are no adjustable properties for class ',class(obj),'.']);
            end%if
            
        end%fcn
        
    end
    
    
    
    methods (Static, Hidden)
        
        %%% ensure column orientation of input vector
        function in = col(in)
        % COL column vector of input
        % 
        % checks if input is of dimension 1*x or x*1 (vector, not matrix)
        % and returns column vector of input.
            
            % get dimension of input argument
            [a,b] = size(in);
            
            % error if input is matrix
            if (a>1) && (b>1)
                error('Input Argument should be vector but is matrix.')
            end%if
            
            % transpose into column vector if input is row vector
            if (b > 1)
                in = in(:);
            end%if
            
        end%fcn
        
        
        %%% check dimension
        function check_dimension_isSVM(pN,value,accepted_type)
        % check if value is scalar/vector/matrix
        
        
            % create string of accepted dimension types
            string_accepted = accepted_type{1};
            for i = 2:length(accepted_type)
                string_accepted = [string_accepted,' or ',accepted_type{i}];
            end%for
            
            
            % value has more than two dimensions?
            if numel(size(value)) > 2
                errorMsg_dimension(pN,string_accepted,'three dimensional array');
            end%if
            
            
            %%% get dimension type of value
            if isscalar(value) % value is scalar?
                actual_type = 'scalar';
            elseif isvector(value)% && ~isscalar(value) % value is vector?
                actual_type = 'vector';
            elseif isempty(value)
                actual_type = 'empty';
            elseif ismatrix(value)% && ~isscalar(value) && ~isvector(value) % value is matrix?
                actual_type = 'matrix';
            end%if
            
            
            %%% compare actual type to accepted types
            flag = false(length(accepted_type),1);
            for i = 1:length(accepted_type)
                
                % if actual type matches accepted type set flag to true
                if strcmp(accepted_type{i},actual_type)
                    flag(i) = true;
                end%if
                
            end%for
            
            
            %%% print error message if actual type is no accepted type
            flag_full = any(flag);
            if ~flag_full
                error(['Invalid dimension of parameter ''',pN,'''. ',...
                    'Should be ',string_accepted,' but seems to be ',actual_type,'!']);
            end%if
            
        end%fcn
        
        
        %%% error mesage
        function errorMsg_xyz(xyz,pN,string_desired,string_is)
            error(['Invalid ',xyz,' of parameter ''',pN,'''. ',...
                'Should be ',string_desired,' but seems to be ',string_is,'!']);
        end%fcn
        
    end
    
    
end%classdef
