classdef lkaSegmentCircle < lkaSegment
%LKASEGMENTCIRCLE 	Create circular street segment.
%	
%	SEG = LKASEGMENTCIRCLE(DELTASET,ANGLESTART,ANGLESTOP,RADIUS) creates a
%	circular street segment of radius RADIUS>0 starting at the angle
%	ANGLESTART and finishing at the angle ANGLESTOP
%	 .) clockwise if ANGLESTART > ANGLESTOP with curvature < 0,
%	 .) counter-clockwise if ANGLESTART < ANGLESTOP with curvature > 0.
%	
%	SEG = LKASEGMENTCIRCLE([],ANGLESTART,ANGLESTOP,RADIUS) applies the
%	default value for DELTASET (see superclass LKASEGMENT).
% 
%	See also LKASEGMENT.
% 

% Subject: lka
% Author: $Author$
% Date: $LastChangedDate$
% Revision: $Revision$



    properties (Constant, Hidden)
        
        designProperties = {'angleStart','angleStop','radius'};
        
    end
    
    
    properties
        
        %%% design data: circular segment
        angleStart;
        angleStop;
        radius;
        
    end%properties
    
    
    properties (Dependent)
        
        %%% info data
        length
        
    end
    
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% METHODS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %%% CONSTRUCTOR
    methods
        
        function obj = lkaSegmentCircle(deltaSet,angleStart,angleStop,radius)
            
            % call superclass constructor
            obj = obj@lkaSegment('circle',deltaSet,[0;0]);
            
            if (angleStart == angleStop)
                warning('MATLAB:lkaSegmentCircle:angleStartEQUALangleStop',...
                    ['Equality of angles ''angleStart'' and ''angleStop'' ',...
                    'is not allowed. Angle of full circumference has ',...
                    'been added to ''angleStop''.']);
                angleStop = angleStop + 2*pi;
            end%if
            
            % set design data of circular segment
            obj.angleStart = angleStart;
            obj.angleStop = angleStop;
            obj.radius = radius;
            
        end%Constructor
        
    end%CONSTRUCTOR-methods
    
    
    %%% GET-Methods
    methods
        
        function value = get.length(obj)
            
            value = obj.radius*abs(obj.angleStop-obj.angleStart);
            
        end%fcn
        
    end%GET-Methods
    
    
    %%% SET-Methods
    methods
        
        function obj = set.angleStart(obj,value)
            
            if (numel(value) ~= 1)
                error('numel(angleStart) ~= 1');
            end%if
            
            % set value
            obj.angleStart = value;
            
        end%fcn
        
        
        function obj = set.angleStop(obj,value)
            
            if (numel(value) ~= 1)
                error('numel(angleStop) ~= 1');
            end%if
            
            % set value
            obj.angleStop = value;
            
        end%fcn
        
        
        function obj = set.radius(obj,value)
            
            if (value <= 0)
                error('radius <= 0');
            elseif (numel(value) ~= 1)
                error('numel(radius) ~= 1');
            end%if
            
            % set value
            obj.radius = value;
            
        end%fcn
        
        
        %%% error at attempt to set dependent property length
        function obj = set.length(obj,~)
            
            errorMsg_SetDependent(obj,'length'); 
            
        end%fcn
          
    end%SET-Methods
    
    
	%%% Implementation of abstract methods
    methods (Access = protected)
        
        function value = getNbrOfPoints_abstract(obj)
            % calc the number of elements of segment to match 'deltaSet'
            
            % the number of required points to match 'deltaSet'
            value = ceil(obj.length/obj.deltaSet) + 1;
         
        end%fcn
        
        
        function value = getEndPoint_abstract(obj)
            
            value = obj.xyStart + ...
                obj.radius*[cos(obj.angleStop) sin(obj.angleStop)] - ...
                obj.radius*[cos(obj.angleStart) sin(obj.angleStart)];
            
        end%fcn
        
        
        function segdat = getSegmentData_abstract(obj)
        % create circular segment based on object data    
            
%             % ensure column-vectors
%             obj.xyStart = lkaSegment.col(obj.xyStart);
            
            % get sign of curvature k !!!!!!! liefert vlt. FALSCHE Ergebnisse
            if obj.angleStop < obj.angleStart; 
                signk = -1;
            else
                signk = 1;
            end%if
            
            % get dependent property 'nbrOfPoints'
            nbrOfPointsDEP = obj.nbrOfPoints;
            
            % angle of discrete points on path relative to angleStart
            phi = linspace(obj.angleStart,obj.angleStop,nbrOfPointsDEP)';
            
            % create circular path
            x = obj.radius*cos(phi);
            y = obj.radius*sin(phi);
            
%             % shift whole path ([x(1);y(1)] matches xyStart)
%             xShift = obj.xyStart(1) - x(1);
%             yShift = obj.xyStart(2) - y(1);
            
%             % output arguments x/y
%             seg.x = x + xShift;
%             seg.y = y + yShift;
            
            % shift whole path ([x(1);y(1)] matches xyStart)
            x = x - x(1) + obj.xyStart(1);
            y = y - y(1) + obj.xyStart(2);
            
            % output argument s/kappa/type
%             seg.s = sort(abs((phi-phi(1))*obj.radius),'ascend');
            s = abs(phi-phi(1))*obj.radius;
            k = signk/obj.radius*ones(nbrOfPointsDEP,1);
            phi = phi + pi/2;
            type = ones(nbrOfPointsDEP,1);
            nbr = ones(nbrOfPointsDEP,1);
            
            % store data in segDat class
            segdat = segDat(x,y,s,k,phi,type,nbr);
            
        end%fcn
        
    end%methods
    
    
end%classdef
