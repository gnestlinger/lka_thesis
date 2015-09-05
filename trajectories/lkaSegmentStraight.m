classdef lkaSegmentStraight < lkaSegment
% LKASEGMENTSTRAIGHT    Create straight street section.
%   _______
%   Syntax:
%   out = LKASEGMENTSTRAIGHT(deltaSet,length,angle)
%   ________________
%   Input arguments:
%   deltaSet ... (opt.) see superclass
%   length ..... length of segment [m]
%   angle ...... angle of segment [rad]
%   _________________
%   Output arguments:
%   see superclass
% 
%   See also LKASEGMENT.
% 

% Subject: lka
% Author: $Author: georgnoname@gmail.com $
% Date: $LastChangedDate: 2015-04-10 09:40:31 +0200 (Fr, 10 Apr 2015) $
% Revision: $Revision: 157 $



    properties (Constant, Hidden)
        
        designProperties = {'length','angle'};
        
    end
    
    
    properties
        
        %%% design data: straight segment
        length; % length of the segment [m]
        angle; % angle (counterclockwise where positive x-axis is 0) [rad]
        
    end%properties
    
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% METHODS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %%% CONSTRUCTOR & Co
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
        
        
    end%methods
    
    
    %%% GET-Methods
    methods
         
    end%GET-Methods
    
    
    %%% SET-Methods
    methods
        
        function obj = set.length(obj,value)
            
            if (value <= 0)
                error('length <= 0');
            end%if
            
            % set value
            obj.length = value;
            
        end%fcn
        
        
        function obj = set.angle(obj,value)
            
            if (numel(value) ~= 1)
                error('numel(angle) ~= 1');
            end%if
            
            % set value
            obj.angle = value;
            
        end%fcn
        
    end%SET-Methods
    
    
    methods (Access = protected)
        
        %%% get the number of segment points
        function value = getNbrOfPoints(obj)
            
            % calc the number of points of segment to match 'deltaSet'
            value = ceil(obj.length/obj.deltaSet) + 1;
            
        end%fcn
        
        
        %%% get endpoint
        function value = getEndPoint(obj)
            
            value = obj.xyStart + obj.length*[cos(obj.angle) sin(obj.angle)];
            
        end%fcn
        
        
        %%% create straight segment based on object data
        function segdat = getSegmentData(obj)
            
            disp('*** straight calculation ***')
            
            % get dependent property 'xyStop'
            xyStopDEP = obj.xyStop;
            
            % get dependent property 'nbrOfPoints'
            nbrOfPointsDEP = obj.nbrOfPoints;
            
            % the segment data
            x = linspace(obj.xyStart(1),xyStopDEP(1),nbrOfPointsDEP)';
            y = linspace(obj.xyStart(2),xyStopDEP(2),nbrOfPointsDEP)';
            s = sqrt((x-obj.xyStart(1)).^2 + (y-obj.xyStart(2)).^2);
            k = zeros(nbrOfPointsDEP,1);
            phi = obj.angle*ones(nbrOfPointsDEP,1);
            type = zeros(nbrOfPointsDEP,1);
            nbr = ones(nbrOfPointsDEP,1);
            
            % store data in segDat class
            segdat = segDat(x,y,s,k,phi,type,nbr);
            
        end%fcn
        
    end%methods
    
    
end%classdef
