classdef Vehicle
%Vehicle	Class defining the Vehicle parameter set.
%
%	The Vehicle parameter set contains the subsets listed below, each of
%	class Vehicle_SubMdl.
%
%	Vehicle Properties:
%	VhclMdl		- Vehicle Model of class Vehicle_SubMdl.
%	SteeringMdl - Steering Model of class Vehicle_SubMdl.
%
%	
%	Vehicle Methods:
%	Vehicle		- Constructor
%
%	See also Vehicle_SubMdl.
	
	properties (SetAccess = private)
		VhclMdl@Vehicle_SubMdl scalar		= SubMdl_empty('Vehicle');
		SteeringMdl@Vehicle_SubMdl scalar	= SubMdl_empty('Steering');
	end
	
% 	properties (Hidden,Transient)
% 		kind2Set = matlab.system.StringSet({'SingleTrackModel','west'});
% 	end
	
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% METHODS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	
	%%% CONSTRUCTOR
	methods
		
		function obj = Vehicle(varargin)
			
			for i = 1:length(varargin)
				obj = addModel(obj,varargin{i});
			end%for
			
		end%fcn
		
	end%CONSTRUCTOR-method
	
	
	methods
		
		function obj = addModel(obj,mdl)
			
			switch mdl.Kind
				
				case 'Vehicle'
					obj.VhclMdl = mdl;
					
				case 'Steering'
					obj.SteeringMdl = mdl;
					
				otherwise
					warning(todo,'An unsupported model kind was provided!');
					
			end%switch
			
		end%fcn
		
	end%methods
	
end%class
