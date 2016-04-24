classdef lkaSegmentClothoid < lkaSegment
%LKASEGMENTCLOTHOID		Create clothoidal street segment.
%	
%	SEG = LKASEGMENTCLOTHOID(DELTASET,CURVSTART,CURVSTOP,SLOPESTART,A)
%	creates a clothoidal shaped street segment with an initial/end
%	curvature CURVSTART/CURVSTOP, an initial slope SLOPESTART and the
%	clothoid parameter A.
%	
%	Interpretation of the clothoid parameter A: a clothiods curvature k is
%	proportional to its length s by k = s/A^2.
%		A > 0: Curvature counter-clockwise
%		A < 0: Curvature clockwise
%	
%	SEG = LKASEGMENTCLOTHOID([],CURVSTART,CURVSTOP,SLOPESTART,A) applies
%	the default value for DELTASET (see superclass LKASEGMENT).
% 
%   See also LKASEGMENT.
% 

% Subject: lka
% Author: $Author$
% Date: $LastChangedDate$
% Revision: $Revision$



    properties (Constant, Hidden)
        
        % user-adjustable design properties for clothoidal street segment
        designProperties = {'curvStart','curvStop','slopeStart','A'};
        
        % clothoid relation of curve length as function of curvature
        sOfCurvature = @(A,curvature) A^2*curvature;
        
        % clothoid relation of tangent angle as function of curvature
        phiOfCurvature = @(A,curvature) A^2*curvature^2/2;
        % k = s/A^2 and phi = s^2/(2*A^2) -> phi = A^2*k^2/2
        
        phiOfLength = @(A,s) s^2/(2*A^2);
        
    end
    
    
    properties
        
        %%% design data: clothoid segment
        curvStart; % curvature at starting point [1/m]
        curvStop; % curvature at endpoint [1/m]
        slopeStart; % slope of clothoide at starting point [rad]
        A; % clothoid parameter [m]
        
    end%properties
    
    
    properties (Dependent)
        
        %%% info data
        length
        
    end
    
    
    properties (Dependent, Hidden, SetAccess = private)
    
        sStart % length of segment at 'curvStart'
        sStop % length of segment at 'curvStop'
        
    end
    
    
    properties (Constant, Hidden)
		
		% squared root of pi
        sqrtPi = 1.7724538509055160272981674833411451827975494561223871282;
        
		% asymptotic points
		P_asymptotic = @(A) A*lkaSegmentClothoid.sqrtPi/2*[1;1];
		
        % clothoid integrand (parameterized form)
        intx = @(t) cos(pi.*t.^2./2);
        inty = @(t) sin(pi.*t.^2./2);
        
        intx_taylor = @(t) 1 ...
            - t^2/factorial(2) + t^4/factorial(4) ...
            - t^6/factorial(6) + t^8/factorial(8) ...
            - t^10/factorial(10) + t^12/factorial(12) ...
            - t^14/factorial(14) + t^16/factorial(16) ...
            - t^18/factorial(18) + t^20/factorial(20);
        
		clothx = @(A,k,t0,t1) abs(A)*lkaSegmentClothoid.sqrtPi*...
			quad(lkaSegmentClothoid.intx,t0,t1)*sign(k);
		clothy = @(A,k,t0,t1) abs(A)*lkaSegmentClothoid.sqrtPi*...
			quad(lkaSegmentClothoid.inty,t0,t1)*sign(A);
        
    end
    
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% METHODS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %%% CONSTRUCTOR & Co
    methods
        
        function obj = lkaSegmentClothoid(deltaSet,curvStart,curvStop,slopeStart,A)
            
            % call superclass constructor
            obj = obj@lkaSegment('clothoid',deltaSet,[0;0]);
            
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
        
        function obj = set.curvStart(obj,value)
            
            if (numel(value) ~= 1)
                error('numel(curvStart) ~= 1');
            end%if
            
            % set value
            obj.curvStart = value;
            
        end%fcn
        
        
        function obj = set.curvStop(obj,value)
            
            if (numel(value) ~= 1)
                error('numel(curvStop) ~= 1');
            end%if
            
            % set value
            obj.curvStop = value;
            
        end%fcn
        
        
        function obj = set.slopeStart(obj,value)
            
            if (numel(value) ~= 1)
                error('numel(slopeStart) ~= 1');
            end%if
            
            % set value
            obj.slopeStart = value;
            
        end%fcn
        
        
        function obj = set.A(obj,value)
            
            if (numel(value) ~= 1)
                error('numel(A) ~= 1');
            end%if
            
			% limit to A > 0
			if value <= 0
				warning('lkaSegmentClothoid:A:valueRange',...
					['Clothoid Parameter ''A'' has to be positive.',...
					'Using the absolute value of the provided value.']);
				value = -value;
			end%if
			
            % set value
            obj.A = value;
            
        end%fcn
        
        
        %%% error at attempt to set dependent property length
        function obj = set.length(obj,~)
            
            errorMsg_SetDependent(obj,'length');
            
        end%fcn
        
    end%SET-Methods
    
    
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


    methods (Access = protected)
        
        %%% get the number of segment points
        function value = getNbrOfPoints(obj)
            % calc the number of points of segment to match 'deltaSet'
            
            value = ceil(obj.length/obj.deltaSet) + 1;
         
        end%fcn
        
        
        %%% get endpoint
        function value = getEndPoint(obj)
            
            value = 'currently unsupported to get that value';
            
        end%fcn
        
        
        %%% create clothoid segment based on object data
        function segdat = getSegmentData(obj)            
            
			% TODO: can be removed since A is limited to A>0
            % get sign of clothoid parameter A
            signA = sign(obj.A);
            absA = abs(obj.A);
            
            % get the sign of clothoid curvature
 			signk = sign(obj.curvStop - obj.curvStart);
            
            % get dependent property 'nbrOfPoints'
            nbrOfPointsDEP = obj.nbrOfPoints;
            
            % pre-calculation
            s = linspace(obj.sStart,obj.sStop,nbrOfPointsDEP)';
            l = s/(absA*obj.sqrtPi)*signk;
            
            % pre-allocation
            cloth.x(nbrOfPointsDEP,1) = 0;
            cloth.y(nbrOfPointsDEP,1) = 0;
            
            % num. integration
			cloth.x(1) = obj.clothx(obj.A,signk,0,l(1));
			cloth.y(1) = obj.clothy(obj.A,signk,0,l(1));
			done_percent_old = 0;
			fprintf('    ');
            for i = 2:nbrOfPointsDEP
				cloth.x(i) = cloth.x(i-1) + obj.clothx(obj.A,signk,l(i-1),l(i));
				cloth.y(i) = cloth.y(i-1) + obj.clothy(obj.A,signk,l(i-1),l(i));
				done_percent = i/nbrOfPointsDEP*100;
				delta = 2;
				if round(done_percent-done_percent_old) > delta
					fprintf('%d%% ',round(done_percent));
					done_percent_old = done_percent;
				end%if
            end%for
			fprintf('%d%% \n',100);
            
            % derivative to compute tangent vector
            % http://mathworld.wolfram.com/TangentVector.html
            xD = absA*obj.sqrtPi*obj.intx(l)*signk;
            yD = absA*obj.sqrtPi*obj.inty(l)*signA;
            tang.x = xD./sqrt(xD.^2+yD.^2);
            tang.y = yD./sqrt(xD.^2+yD.^2);
            
            % rotation angle, curvStart ~= 0 results in a tangent vector
            % with slope ~= 0 at start point
            slopeStartDue2curvStart = angle(tang.x(1) + 1i*tang.y(1)); % falls curvStart ~= 0
            slopeStartDue2curvStart2 = mod(obj.phiOfCurvature(obj.A,obj.curvStart),2*pi);
            if slopeStartDue2curvStart == slopeStartDue2curvStart2
                disp('    is equal')
            else 
                disp('    is NOT equal')
            end%if
            alph = (obj.slopeStart - slopeStartDue2curvStart);
            
            % rotate clothoid
            cloth.rot.x = [cloth.x,cloth.y]*obj.rotMatX(alph);
            cloth.rot.y = [cloth.x,cloth.y]*obj.rotMatY(alph);
            
            % rotate tangent
            tang.rot.x = [tang.x,tang.y]*obj.rotMatX(alph);
            tang.rot.y = [tang.x,tang.y]*obj.rotMatY(alph);
            
            % shift whole trajectory ([x(1);y(1)] matches xyStart)
            xShift = obj.xyStart(1) - cloth.rot.x(1);
            yShift = obj.xyStart(2) - cloth.rot.y(1);
            
            % output arguments
            x = cloth.rot.x + xShift;
            y = cloth.rot.y + yShift;
            sCloth = sort(abs(s));
            sOut = sCloth - sCloth(1);
            k = s/absA^2*signA;
            phi = unwrap(angle(tang.rot.x + 1i*tang.rot.y));%(s).^2/(2*obj.A^2);
            type = 2*ones(nbrOfPointsDEP,1);
            nbr = ones(nbrOfPointsDEP,1);
            
            % store data in segDat class
            segdat = segDat(x,y,sOut,k,phi,type,nbr);
            
        end%fcn
        
    end%methods
    
    
end%classdef
