classdef SteerMdl_CarMaker_DSR < Vehicle_SubMdl
	%UNTITLED3 Summary of this class goes here
	%   Detailed explanation goes here
	
% Subject: lka
% $Author$
% $LastChangedDate$
% $Revision$

	
	properties (SetAccess = private)
		Parameters%	= VhclMdl_SingleTrack_parameters();
		InitValues%	= VhclMdl_SingleTrack_initValues();
		StateNames%	= VhclMdl_SingleTrack_stateNames();
	end
	
	
	properties (Constant,GetAccess = protected)
		validClass_Parameters = 'SteerMdl_CarMaker_DSR_parameters';
		validClass_InitValues = 'SteerMdl_CarMaker_DSR_initValues';
		validClass_StateNames = 'SteerMdl_CarMaker_DSR_stateNames';
		nbrOfParameters = 10;
		nbrOfStates = 2;
	end
	
	
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% METHODS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	
	%%% CONSTRUCTOR
	methods (Hidden)
		
		%%% Constructor
		function obj = SteerMdl_CarMaker_DSR(param,initValues,stateNames,info)
			
			% handle input arguments
			if nargin < 4; info = ''; end%if
			
			% call superclass constructor
            obj = obj@Vehicle_SubMdl('Steering','CarMaker_DSR',info);
			
			if nargin < 1
				
			else
				
				% set abstract properties of superclass
				obj.Parameters = param;
				obj.InitValues = initValues;
				obj.StateNames = stateNames;
				
			end%if
			
		end%fcn
		
	end%CONSTRUCTOR-methods
	
	
	%%% SET-Methods
	methods
		
	end% SET-Methods
	
	
end%classdef

