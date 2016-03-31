classdef VhclMdl_SingleTrack < Vehicle_SubMdl
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
	
	
	properties (Constant, GetAccess = protected)
		validClass_Parameters = 'VhclMdl_SingleTrack_parameters';
		validClass_InitValues = 'VhclMdl_SingleTrack_initValues';
		validClass_StateNames = 'VhclMdl_SingleTrack_stateNames';
		nbrOfParameters = 6;
		nbrOfStates = 4;
	end
	
	
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% METHODS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	
	
	%%% CONSTRUCTOR
	methods (Hidden)
		
		%%% Constructor
		function obj = VhclMdl_SingleTrack(param,initValues,stateNames,info)
			
			% handle input arguments
			if nargin < 4; info = ''; end%if

			% call superclass constructor
			obj = obj@Vehicle_SubMdl('Vehicle','SingleTrackModel',info);
				
			if nargin < 1
				
			else
				
				% set abstract properties of superclass
				obj.Parameters = param;
				obj.InitValues = initValues;
				obj.StateNames = stateNames;
			
			end%if
			
		end%fcn
		
	end%CONSTRUCTOR-method
	
	
	%%% SET-Methods
	methods
		
		function obj = set.Parameters(obj,value)
			
			if isa(value,obj.validClass_Parameters)
				obj.Parameters = value;
			else
				obj.error_validClass(obj.validClass_Parameters,class(value));
			end%if
			
		end%fcn
		
		function obj = set.InitValues(obj,value)
			
			if isa(value,obj.validClass_InitValues)
				obj.InitValues = value;
			else
				obj.error_validClass(obj.validClass_Parameters,class(value));
			end%if
			
		end%fcn
		
		function obj = set.StateNames(obj,value)
			
			if isa(value,obj.validClass_StateNames)
				obj.StateNames = value;
			else
				obj.error_validClass(obj.validClass_Parameters,class(value));
			end%if
			
		end%fcn
		
	end% SET-Methods
	
	
end%class

