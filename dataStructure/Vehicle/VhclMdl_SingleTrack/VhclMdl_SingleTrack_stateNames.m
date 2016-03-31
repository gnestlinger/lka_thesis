classdef VhclMdl_SingleTrack_stateNames
	%UNTITLED3 Summary of this class goes here
	%   Detailed explanation goes here
	
% Subject: lka
% $Author$
% $LastChangedDate$
% $Revision$


	properties (SetAccess = private)
		x1@char vector = 'sy'
		x2@char vector = 'vy'
		x3@char vector = 'yaw rate'
		x4@char vector = 'yar acceleration'
	end
	
	
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% METHODS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	
	
	%%% CONSTRUCTOR
	methods (Hidden)
		
		%%% Constructor
		function obj = VhclMdl_SingleTrack_stateNames()
			% nothing to do, use the default values!
		end%fcn
		
	end%CONSTRUCTOR-method
	
end%classdef

