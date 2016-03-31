classdef Vehicle_SubMdl
%Vehicle_SubMdl		Class defining the interface for each vehicle submodel.
%
%	Each property of a Vehicle object contains the properties listed below.
%	
%	
%	Vehicle_SubMdl Properties:
%	Kind		- String identifying the model.
%	Info		- Some additional info string.
%	Parameters	- The parameters of the model.
%	InitValues	- The initial values of the model.
%	StateNames	- The names of the states of the model.
%	
%	
%	Vehicle_SubMdl Methods:
%	Vehicle_SubMdl	- Constructor.
%	
%	See also Vehicle.
	
% Subject: lka
% $Author$
% $LastChangedDate$
% $Revision$

	
	properties (SetAccess = private)
		Kind@char vector = '';
		Type@char vector = '';
	end

	
	properties (SetAccess = public)
		Info@char vector = '';
	end
	
	
	properties (Abstract,SetAccess = private)
		Parameters
		InitValues
		StateNames
	end
	
	
	properties (Abstract,Constant,GetAccess = protected)
		validClass_Parameters
		validClass_InitValues
		validClass_StateNames
		nbrOfStates
		nbrOfParameters
	end
	
	
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% METHODS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	
	%%% CONSTRUCTOR
	methods
		
		function obj = Vehicle_SubMdl(kind,type,info)
			
			obj.Kind = kind;
			obj.Type = type;
			obj.Info = info;
		
		end%fcn
		
	end%CONSTRUCTOR-methods
	
	
	methods (Static,Hidden)
		
		function Parameters = setParameterFromStruct(Parameters,S)
			
			if ~isa(S,'struct')
				error('Class error!');
			else
				try
					fn = fieldnames(Parameters);
					for i = 1:numel(fn)
						Parameters.(fn{i}) = S.(fn{i});
					end%for
				catch exception
					error(exception.message);
				end%try
			end%if

		end%fcn
		
		
		function error_validClass(validClass,currentClass)
			
			error('Current class ''%s'' does not match valid class ''%s''',...
				validClass,currentClass);
			
		end%fcn
		
	end%methods
	
end%class

