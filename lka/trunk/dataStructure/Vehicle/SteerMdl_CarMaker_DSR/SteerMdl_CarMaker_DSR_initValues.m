classdef SteerMdl_CarMaker_DSR_initValues
	%UNTITLED3 Summary of this class goes here
	%   Detailed explanation goes here
	
% Subject: lka
% $Author$
% $LastChangedDate$
% $Revision$

	
	properties (SetAccess = public)
		x1@double scalar = 0
		x2@double scalar = 0
	end
	
	
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% METHODS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	
	
	%%% CONSTRUCTOR
	methods (Hidden)
		
		%%% Constructor
		function obj = SteerMdl_CarMaker_DSR_initValues(x10,x20)
			
			if nargin == 0
				% use default values
				
			else
				obj.x1 = x10;
				obj.x2 = x20;
			end%if
			
		end%fcn
		
	end%CONSTRUCTOR-method
	
end%class

