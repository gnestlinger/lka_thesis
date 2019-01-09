classdef lkaSegmentClothoid < lkaSegment
%LKASEGMENTCLOTHOID		Create clothoidal street segment.
%	
%	OBJ = LKASEGMENTCLOTHOID(DELTASET,CURVSTART,CURVSTOP,SLOPESTART,A)
%	creates a clothoidal shaped street segment OBJ with an initial/end
%	curvature CURVSTART/CURVSTOP, an initial slope SLOPESTART and the
%	clothoid parameter A>0.
%	
%	OBJ = LKASEGMENTCLOTHOID([],CURVSTART,CURVSTOP,SLOPESTART,A) applies
%	the default value for DELTASET (see superclass LKASEGMENT).
%	
%	Interpretation of the clothoid parameter A: a clothiods curvature k is
%	proportional to its length s by k = s/A^2.
%	
%	A positive/negative curvature yields a path turning
%	counter-clockwise/clockwise.
%	
%	See also LKASEGMENT.
% 

% Subject: lka
% Author: $Author$
% Date: $LastChangedDate$
% Revision: $Revision$



	properties (Constant, Hidden = true)
		
		% designProperties - User adjustable properties.
		%	Design the street segment LKASEGMENTCLOTHOID by adjusting its
		%	properties CURVSTART, CURVSTOP, SLOPESTART and A.
		designProperties = {'curvStart','curvStop','slopeStart','A'};
		
	end%properties
	
	properties (Constant, Hidden = true)
        % clothoid relation of curve length as a function of curvature
        sOfCurvature = @(A,curvature) A^2*curvature;
        
        % clothoid relation of tangent angle as a function of curvature
        phiOfCurvature = @(A,curvature) A^2*curvature^2/2;
        % k = s/A^2 and phi = s^2/(2*A^2) -> phi = A^2*k^2/2
        
		% clothoid relation of tangent angle as a function of length
        phiOfLength = @(A,s) s^2/(2*A^2);
		
		% squared root of pi
		sqrtPi = 1.7724538509055160272981674833411451827975494561223871282;
		
		% asymptotic points
		P_asymptotic = @(A) A*lkaSegmentClothoid.sqrtPi/2*[1;1];
		
	end%properties
    
    
    properties (SetAccess = private)% design data: clothoid segment
        
        curvStart(1,1) double {mustBeFinite};	% curvature at starting point [1/m]
        curvStop(1,1) double {mustBeFinite};	% curvature at endpoint [1/m]
        slopeStart(1,1) double {mustBeFinite}; % slope of clothoide at starting point [rad]
        A(1,1) double {mustBeFinite,mustBePositive} = 1; % clothoid parameter [m]
        
    end%properties
    
    
    properties (Dependent, SetAccess = private)
        
        %%% info data
        length
        
    end%properties
    
    
    properties (Dependent, Hidden, SetAccess = private)
    
        sStart % length of segment at 'curvStart'
        sStop % length of segment at 'curvStop'
        
	end%properties
    
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% METHODS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %%% CONSTRUCTOR
    methods
        
        function obj = lkaSegmentClothoid(deltaSet,curvStart,curvStop,slopeStart,A)
            
            % call superclass constructor
            obj = obj@lkaSegment("clothoid",deltaSet,[0;0]);
            
            % set design data of clothoid segment
            obj.curvStart = curvStart;
            obj.curvStop = curvStop;
            obj.slopeStart = slopeStart;
            obj.A = A;
            
        end%Constructor
        
    end%CONSTRUCTOR-methods
    
    
    %%% GET-Methods
    methods
            
        function value = get.length(obj)
        % get the length of street segment
            
            value = abs(obj.sStop - obj.sStart);
					
        end%fcn
        
        
        function value = get.sStart(obj)
        % get the length s at starting curvature
            
%             value = obj.A^2*obj.curvStart;
            value = obj.sOfCurvature(obj.A,obj.curvStart);
         
        end%fcn
        
        
        function value = get.sStop(obj)
        % get the length s at starting curvature
            
%             value = obj.A^2*obj.curvStop;
            value = obj.sOfCurvature(obj.A,obj.curvStop);
         
        end%fcn
        
    end%GET-Methods
    
    
    %%% SET-Methods
    methods
    end%SET-Methods
    
	
    %%% User-facing methods
    methods
        
        function plotasypoints(obj)
            
            plot(obj);
			Pasy = obj.P_asymptotic(obj.A);
			
			if obj.curvStop < 0
				Pasy(2) = -Pasy(2);
			end%if
			
            hold on
            plot(Pasy(1),Pasy(2),'rx');
            hold off 
        
        end
        
    end%methods
    
    
	%%% Implementation of abstract methods
    methods (Access = protected)
        
		function obj = rotate_abstract(obj,phi)
			
			obj	= shift(obj, obj.rotMat(phi)*obj.xyStart' );
			obj.slopeStart = obj.slopeStart + phi;
			
		end%fcn
		
        function value = getNbrOfPoints_abstract(obj)
            % calc the number of points of segment to match 'deltaSet'
            
            value = ceil(obj.length/obj.deltaSet) + 1;
         
        end%fcn
        
        
        function value = getEndPoint_abstract(obj)
            
            value = 'currently unsupported to get that value';
            
        end%fcn
        
        
		function segdat = getSegmentData_abstract(obj)
			
			% get the sign of clothoid curvature
			signk = sign(obj.curvStop - obj.curvStart);
			
			% get dependent property 'nbrOfPoints'
			nbrOfPointsDEP = obj.nbrOfPoints;
			
			% pre-calculation
			s = linspace(obj.sStart,obj.sStop,nbrOfPointsDEP)*signk;
			
			% numerical integration of clothoid coordinates
			[x,y] = lkaSegmentClothoid.clothoid_numInt(s,obj.A);
			
			% create segDat object
			segdat = segDat(...
				x,...					% x coordinate
				y,...					% y coordinate
				s-s(1),...				% curve length
				s/obj.A^2,...			% k = s/A^2
				s.^2/(2*obj.A^2),...	% phi = A^2*k^2/2 = s^2/(2*A^2)
				2,...					% segment type
				1);						% segment number
			
			if signk < 0
				segdat = changeSignOfCurvature(segdat);
			end%if
			
			% rotate clothoid
			segdat = rotate(segdat,obj.slopeStart-segdat.phi(1));
			
			% shift whole trajectory so [x(1);y(1)] matches xyStart
			segdat = shiftTo(segdat,obj.xyStart);
			
			% just for debugging
% 			plottangent(segdat,[1,length(segdat.x)]);
			
		end%fcn
        
    end%methods
	
	
	
	methods (Static)%, Access = private)
		
		function [x,y] = clothoid_numInt(s,A,doPlot)
		% CLOTHOID_NUMINT	Numerical integration of clothoid.
		%	[X,Y] = CLOTHOID_NUMINT(S,A) calculates the clothoid points
		%	Y(X) for the vector of clothoid length S using the chlothoid
		%	parameter A>0.
		% 
		%	The input argument S must be vector sized, must not contain
		%	equal elements and must be strictly monotonically increasing or
		%	decreasing!
		% 
		% 
		%	[X,Y] = CLOTHOID_NUMINT(...,DOPLOT) enables some plot output by
		%	setting DOPLOT to logical 1.
		% 
		%	[x;y] = A*sqrt(pi)*Int([cos(pi*t^2/2);sin(pi*t^2/2)],...
		%	S(1)/(sqrt(pi)*A),S(end)/(sqrt(pi)*A))
			
			
			% check input arguments
			if ~isvector(s)
				error('Input argument S has to be vector-sized!');
			elseif any(diff(s) == 0)
				error('Input argument S must not contain equal elements!');
			elseif abs(sum(sign(diff(s)))) ~= length(s)-1
				error('Input argument S must be strictly monotonically decreasing/increasing!');
			else
				% ensure row vector orientation
				s = s(:)';
			end%if
			
			if A <= 0
				error('Input argument A has to be positive!')
			end%if
			
			if nargin < 3
				doPlot = false;
			end%if
			
			
			% clothoid integrands
			integrandX = @(t) cos(pi.*t.^2./2);
			integrandY = @(t) sin(pi.*t.^2./2);
			
			clothX = @(A,sStart,sStop) A*lkaSegmentClothoid.sqrtPi*...
				quad(integrandX,...
				sStart/(lkaSegmentClothoid.sqrtPi*A),...
				sStop/(lkaSegmentClothoid.sqrtPi*A));
			clothY = @(A,sStart,sStop) A*lkaSegmentClothoid.sqrtPi*...
				quad(integrandY,...
				sStart/(lkaSegmentClothoid.sqrtPi*A),...
				sStop/(lkaSegmentClothoid.sqrtPi*A));
			
			% number of points 
			nbrOfPoints = length(s);
			
			% pre-allocation
			x(1,nbrOfPoints) = 0;
			y(1,nbrOfPoints) = 0;
			
			% num. integration
			x(1) = clothX(A,0,s(1));
			y(1) = clothY(A,0,s(1));
			done_percent_old = 0;
			fprintf('    ');
			for i = 2:nbrOfPoints
				x(i) = x(i-1) + clothX(A,s(i-1),s(i));
				y(i) = y(i-1) + clothY(A,s(i-1),s(i));
				done_percent = i/nbrOfPoints*100;
				delta = 10;
				if (done_percent-done_percent_old) > delta
					fprintf('%d%% ',round(done_percent));
					done_percent_old = done_percent;
				end%if
			end%for
			fprintf('%d%% \n',100);
			
			% curvature and tangent angle
			k = s./A^2;
			phi = s.^2/(2*A^2);
			
% 			% derivative to compute tangent vector
% 			% http://mathworld.wolfram.com/TangentVector.html
% 			xD = A*sqrt(pi)*integrandX(s/(sqrt(pi)*A));
% 			yD = A*sqrt(pi)*integrandY(s/(sqrt(pi)*A));
% 			x_tang = xD./sqrt(xD.^2+yD.^2);
% 			y_tang = yD./sqrt(xD.^2+yD.^2);
% 			phi_check = atan2(y_tang,x_tang);
% 			
			% compare phi/phi_check
			phi_atan2 = atan2(sin(phi),cos(phi));
% 			disp(['    phi==phi_check: ',...
% 				num2str(all(abs(phi_atan2 - phi_check) < 1e-12))]);
			
			if doPlot
				
				sd = segDat(x,y,s,k,phi,2,1);
				
				subplot(2,1,1)
				plottangent(sd,[1,length(sd.x)]);
				hold on; plot(x(1),y(1),'o'); hold off;
				grid on; 
				axis equal;
				
				subplot(2,1,2)
				ph = plot(...
					s,phi,...
					s,phi_atan2);
% 					s,phi_check,'.');
				grid on
				set(ph(1),'LineWidth',2);
				set(ph(2),'color','g');
				legend('phi','phi\_atan2','location','Best');
				
			end%if
			
		end%fcn
		
	end%methods
	
	
	%%% Test-Methods
	methods (Static, Hidden)
		
		function test_curvatureSettings()
			
			% signk = sign(obj.curvStop - obj.curvStart)
			
			% signk > 0
			c1p = lkaSegmentClothoid(2,-0.2,0.2,0,25);
			
			% signk < 0
			c1n = lkaSegmentClothoid(2,0.2,-0.2,0,25);
			
			lkaSegmentClothoid.test_curvatureSettings_sub(c1p,1);
			lkaSegmentClothoid.test_curvatureSettings_sub(c1n,2);
			
			
		end%fcn
		
		function test_curvatureSettings_sub(obj,fign)
			
			figure(fign)
			subplot(3,2,[1 3 5])
			plot(obj);
			
			subplot(3,2,2)
			plotdiff_(obj.segmentData,[],'s');
			title('')
			
			subplot(3,2,4)
			plotdiff_(obj.segmentData,[],'k');
			title('')
			
			subplot(3,2,6)
			plotdiff_(obj.segmentData,[],'phi');
			title('')
			
		end%fcn
		
	end%methods
	
end%classdef
