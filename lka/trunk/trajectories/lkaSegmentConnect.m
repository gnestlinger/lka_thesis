classdef lkaSegmentConnect < lkaSegment
%LKASEGMENTCONNECT	Create connected segments.
%	
%	--- (just used by subclasse-methods) ---------------------------------
%	SEG = LKASEGMENTCONNECT(SEGDAT) stores the segment data SEGDAT object.
%	----------------------------------------------------------------------
%	
%	See also LKASEGMENT, LKASEGMENTSTRAIGHT, LKASEGMENTCIRCLE,
%	LKASEGMENTCLOTHOID.

% Subject: lka
% Author: $Author$
% Date: $LastChangedDate$
% Revision: $Revision$



	properties (Constant, Hidden = false)
		
		% designProperties - User adjustable properties.
		%	Since LKASEGMENTCONNECT just holds the street segment data of
		%	connected street segments, there are no adjustable design
		%	properties.
		designProperties = {};
		
    end
    
    
    properties (Dependent)
        
        length
        
    end
    
    
    properties (Hidden, SetAccess = private)
        
        nbrOfPoints_stored;
        segmentData_stored;
        
	end
    
	
	
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	%%% METHODS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %%% CONSTRUCTOR
    methods
        
        function obj = lkaSegmentConnect(segmentData)
            
            % call superclass constructor
            obj = obj@lkaSegment('connected',NaN,[0,0]);
            
            obj.nbrOfPoints_stored = length(segmentData.x);
			obj.segmentData_stored = segmentData;
			
        end%Constructor
		
    end%CONSTRUCTOR-methods
    
    
    %%% User-facing methods
    methods
         
        %%% redefine superclass-method
        function obj = shift(obj,point)
        %SHIFT  Shift the street segment.
		%	See also LKASEGMENT/SHIFT.
        
        
            if nargin < 2
                point = [0,0];
            end%if
			
			% call superlass method to set property xyStart
			obj = shift@lkaSegment(obj,point);
			
			% shift segDat object manually
			sd = shiftTo(obj.segmentData,point);
			
			% recreate LKASEGMENTCONNECT object by calling its constructor
			obj = lkaSegmentConnect(sd);
            
        end%fcn
        
        
        %%% redefine superclass-method
        function obj = resample(obj,deltaNew)
        %RESAMPLE   Apply a set distance between points.
        %   For objects of class LKASEGMENTCONNECT, resampling the
        %   connected street segment is not possible.
        
           
            msg = ['You are trying to resample a connected street segment ',...
                'with a current average delta of %.2f ',...
                'with a new delta of %.2f.\n'];
            msg = sprintf(msg,obj.deltaAct,deltaNew);
            
            errmsg = ['For objects of class lkaSegmentConnect, ',...
                'resampling the street segment is not allowed!'];
            error([msg,errmsg]);
            
        end%fcn
         
    end%methods
	
	
	%%% GET-Methods
    methods
        
        function value = get.length(obj)
%             disp('getting length...')
            value = obj.segmentData.s;
			
			% avoid multiple calls of superclass method for dependent
			% property 'segmentData' by getting the last element from
			% buffered data 'value'
			value = value(end);
			
% 			value = obj.segmentDataConnected.s(end);
%             disp('...done')
        end%fcn

    end%GET-Methods
    
    
    %%% SET-Methods
    methods
        
        %%% error at attempt to set dependent property length
        function obj = set.length(obj,~)
            
            errorMsg_SetDependent(obj,'length'); 
            
        end%fcn
        
        
    end%SET-Methods
    
    
	%%% Implementation of abstract methods
    methods (Access = protected)
        
        function value = getNbrOfPoints_abstract(obj)
% 			disp('getting NbrOfPoints...')
			value = obj.nbrOfPoints_stored;
% 			disp('...done')
        end%fcn
        
        
        function value = getEndPoint_abstract(obj)
%             disp('getting endPoint...')
			sd = obj.segmentData;
			
			% avoid multiple calls of superclass method for dependent
			% property 'segmentData' by getting the last element from
			% buffered data 'sd'
            value = [sd.x(end), sd.y(end)];
%             disp('...done')
        end%fcn
        
            
        function segdat = getSegmentData_abstract(obj)
            segdat = obj.segmentData_stored;
        end%fcn
        
    end%methods
    
        
end%classdef
