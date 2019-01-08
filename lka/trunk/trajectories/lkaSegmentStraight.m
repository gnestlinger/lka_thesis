classdef lkaSegmentStraight < lkaSegment
%LKASEGMENTSTRAIGHT		Create straight street segment.
%
%	OBJ = LKASEGMENTSTRAIGHT(DELTASET,LENGTH,ANGLE) creates a straight
%	street segment OBJ of length LENGTH and orientation ANGLE. Each two
%	nearby street segment points have a distance of at most DELTASET.
%	
%	OBJ = LKASEGMENTSTRAIGHT([],LENGTH,ANGLE) applies the default value
%	for DELTASET (see superclass LKASEGMENT).
%
%	See also LKASEGMENT.
% 

% Subject: lka
% Author: $Author$
% Date: $LastChangedDate$
% Revision: $Revision$



	properties (Constant, Hidden = false)
		
		% designProperties - User adjustable properties.
		%	Design the street segment LKASEGMENTSTRAIGHT by adjusting its
		%	properties LENGTH and ANGLE.
		designProperties = {'length','angle'};
		
	end%properties
    
    
    properties (SetAccess = private)% design data: straight segment
        
        length; % length of the segment [m]
        angle; % angle (counterclockwise where positive x-axis is 0) [rad]
        
    end%properties
    
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% METHODS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %%% CONSTRUCTOR
    methods
        
        %%% Constructor
        function obj = lkaSegmentStraight(deltaSet,length,angle)
            
            % call superclass constructor
            obj = obj@lkaSegment('straight',deltaSet,[0;0]);
            
            % set length and angle
            obj.length = length;
            obj.angle = angle;
            
        end%Constructor
        
        
%         function obj = laneSep(obj,width,varargin)
%             
%             xyStart = obj.xyStart' - width/2*[cos(obj.angle+pi/2);sin(obj.angle+pi/2)];
%             rSep = lkaSegmentStraight(obj.delta,xyStart,obj.length,obj.angle);
%             
%             xyStart = obj.xyStart' + width/2*[cos(obj.angle+pi/2);sin(obj.angle+pi/2)];
%             lSep = lkaSegmentStraight(obj.delta,xyStart,obj.length,obj.angle);    
%             
%         end%fcn
        
        
    end%CONSTRUCTOR-methods
    
    
    %%% GET-Methods
    methods
    end%GET-Methods
    
    
    %%% SET-Methods
    methods	
    end%SET-Methods
    
    
	%%% Implementation of abstract methods
    methods (Access = protected)
        
		function obj = rotate_abstract(obj,phi)
			
			obj	= shift(obj, obj.rotMat(phi)*obj.xyStart' );
			obj.angle = obj.angle + phi;
			
		end%fcn
		
		
        function value = getNbrOfPoints_abstract(obj)
            
            % calc the number of points of segment to match 'deltaSet'
            value = ceil(obj.length/obj.deltaSet) + 1;
            
        end%fcn
        
        
        function value = getEndPoint_abstract(obj)
        % get endpoint
		
            value = obj.xyStart + obj.length*[cos(obj.angle) sin(obj.angle)];
            
        end%fcn
        
        
        function segdat = getSegmentData_abstract(obj)
		% create straight segment based on object data
            
            % get dependent property 'xyStop'
            xyStopDEP = obj.xyStop;
            
            % get dependent property 'nbrOfPoints'
            nbrOfPointsDEP = obj.nbrOfPoints;
            
            % the segment data
            x = linspace(obj.xyStart(1),xyStopDEP(1),nbrOfPointsDEP);
            y = linspace(obj.xyStart(2),xyStopDEP(2),nbrOfPointsDEP);
            s = sqrt((x-obj.xyStart(1)).^2 + (y-obj.xyStart(2)).^2);
            k = zeros(1,nbrOfPointsDEP);
            phi = obj.angle*ones(1,nbrOfPointsDEP);
%             type = zeros(1,nbrOfPointsDEP);
%             nbr = ones(1,nbrOfPointsDEP);
            
            % store data in segDat class
            segdat = segDat(x,y,s,k,phi,0,1);
            
        end%fcn
        
    end%methods
    
    
end%classdef
