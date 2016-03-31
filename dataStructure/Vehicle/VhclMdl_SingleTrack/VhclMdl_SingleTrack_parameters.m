classdef VhclMdl_SingleTrack_parameters
	%UNTITLED3 Summary of this class goes here
	%   Detailed explanation goes here
	
% Subject: lka
% $Author$
% $LastChangedDate$
% $Revision$

	
	properties (SetAccess = public)
		csv		= []
		csh		= []
		lv		= []
		lh		= []
		m		= []
		Iz		= []
	end
	
	
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% METHODS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	
	
	%%% CONSTRUCTOR
	methods (Hidden)
		
		%%% Constructor
		function obj = VhclMdl_SingleTrack_parameters(csv,csh,lv,lh,m,Izz)
			
			if nargin == 0
				% use default values
				
			elseif nargin == 1
				% input argument has to be struct
				
				if ~isa(csv,'struct')
					error('Class');
				else					
					obj = Vehicle_SubMdl.setParameterFromStruct(obj,csv);
				end%if
				
			else
				% all parameter have to be provided
				
				obj.csv		= csv;
				obj.csh		= csh;
				obj.lv		= lv;
				obj.lh		= lh;
				obj.m		= m;
				obj.Iz		= Izz;
			end%if
			
		end%fcn
		
	end%CONSTRUCTOR-method
	
end%classdef

