classdef VhclMdl_SingleTrack_initValues
	%UNTITLED3 Summary of this class goes here
	%   Detailed explanation goes here
	
% Subject: lka
% $Author$
% $LastChangedDate$
% $Revision$


	properties (SetAccess = public)
		x1@double scalar = 0
		x2@double scalar = 0
		x3@double scalar = 0
		x4@double scalar = 0
	end
	
	
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% METHODS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	
	
	%%% CONSTRUCTOR
	methods (Hidden)
		
		%%% Constructor
		function obj = VhclMdl_SingleTrack_initValues(x10,x20,x30,x40)
			
			if nargin == 0
				% use default values
			else
				obj.x1 = x10;
				obj.x2 = x20;
				obj.x3 = x30;
				obj.x4 = x40;
			end%if
			
		end%fcn
		
	end%CONSTRUCTOR-method
	
end%class

