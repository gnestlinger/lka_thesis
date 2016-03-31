classdef SteerMdl_CarMaker_DSR_stateNames
	%UNTITLED3 Summary of this class goes here
	%   Detailed explanation goes here
	
% Subject: lka
% $Author$
% $LastChangedDate$
% $Revision$
	
	
	properties (SetAccess = private)
		x1@char vector = 'steering wheel angle'
		x2@char vector = 'steering wheel angular velocity'
	end
	
	
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% METHODS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	
	
	%%% CONSTRUCTOR
	methods (Hidden)
		
		%%% Constructor
		function obj = SteerMdl_CarMaker_DSR_stateNames()
			% nothing to do, use the default values!
		end%fcn
		
	end%CONSTRUCTOR-method
	
end%classdef

