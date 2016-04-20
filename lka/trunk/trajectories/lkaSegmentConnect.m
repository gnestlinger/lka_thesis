classdef lkaSegmentConnect < lkaSegment
%LKASEGMENTCONNECT	Create connected segments.
%	
%	--- (just used by subclasse-methods) ---------------------------------
%	SEG = LKASEGMENTCONNECT(SEGDAT) stores the segment data SEGDAT object.
%	----------------------------------------------------------------------
%	
%   See also LKASEGMENT, LKASEGMENTSTRAIGHT, LKASEGMENTCIRCLE,
%   LKASEGMENTCLOTHOID.

% Subject: lka
% Author: $Author$
% Date: $LastChangedDate$
% Revision: $Revision$



    properties (Constant, Hidden)
        
        designProperties = {};
        
    end
    
    
    properties(Dependent)
        
        length
        
    end
    
    
    properties (Hidden, SetAccess = private)
        
        nbrOfPointsConnected
        
    end
    
    
    properties (Hidden, SetAccess = private)
        
        %%% design data: straight segment
        segmentDataConnected; % see superclass
        
    end
    
    
    %%% CONSTRUCTOR
    methods
        
        function obj = lkaSegmentConnect(segmentData)
            
            % call superclass constructor
            obj = obj@lkaSegment('connected',NaN,[0,0]);
            
            obj.segmentDataConnected = segmentData;
            obj.nbrOfPointsConnected = length(obj.segmentDataConnected.x);
            
        end%Constructor
		
    end%CONSTRUCTOR-methods
    
    
    %%% User-facing methods
    methods
         
        %%% redefine superclass-method
        function obj = shift(obj,point)
        %SHIFT  Shift the street segment.
        %   For objects of class LKASEGMENTCONNECT, shifting the starting
        %   point of the street segment is not possible.
        %
        
        
            if nargin < 2
                point = [0,0];
            end%if
            
            msg = ['You are trying to shift initial point ',...
                'from (%.1f,%.1f) to (%.1f,%.1f).\n'];
            msg = sprintf(msg,obj.xyStart,point);
            
            errmsg = ['For objects of class lkaSegmentConnect, ',...
                'shifting the starting point is not allowed!'];
            error([msg,errmsg]);
            
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
            
            value = obj.segmentData.s(end);
            
        end%fcn

    end%GET-Methods
    
    
    %%% SET-Methods
    methods
        
        %%% error at attempt to set dependent property length
        function obj = set.length(obj,~)
            
            errorMsg_SetDependent(obj,'length'); 
            
        end%fcn
        
        
    end%SET-Methods
    
    
    methods (Access = protected)
        
        function value = getNbrOfPoints(obj)
            
            value = obj.nbrOfPointsConnected;
            
        end%fcn
        
        
        function value = getEndPoint(obj)
            
            value = [obj.segmentData.x(end), obj.segmentData.y(end)];
            
        end%fcn
        
            
        function segdat = getSegmentData(obj)
            
            segdat = obj.segmentDataConnected;
            
        end%fcn
        
    end%methods
    
        
end%classdef
