classdef SteerMdl_CarMaker_DSR_parameters
	%UNTITLED3 Summary of this class goes here
	%   Detailed explanation goes here
	
% Subject: lka
% $Author$
% $LastChangedDate$
% $Revision$

	
	properties
		J		= [];
		V		= [];
		d_rot	= [];
		d_rack	= [];
		i		= [];
		m_right	= [];
		m_left	= [];
		m_rack	= [];
		alph	= [];
		xi		= [];
	end
	
	
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% METHODS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	
	%%% CONSTRUCTOR
	methods (Hidden)
		
		%%% Constructor
		function obj = SteerMdl_CarMaker_DSR_parameters(J,V,d_rot,d_rack,i,m_left,m_right,m_rack,alph,xi)
			
			if nargin == 0
				% use default values
				
			elseif nargin == 1
				% input argument has to be struct
				
				if ~isa(J,'struct')
					error('Class');
				else
					obj = Vehicle_SubMdl.setParameterFromStruct(obj,J);
				end%if
				
			else
				% all parameter have to be provided
				
				obj.J		= J;
				obj.V		= V;
				obj.d_rot	= d_rot;
				obj.d_rack	= d_rack;
				obj.i		= i;
				obj.m_left	= m_left;
				obj.m_right	= m_right;
				obj.m_rack	= m_rack;
				obj.alph	= alph;
				obj.xi		= xi;
				
			end%if
			
		end%fcn
		
	end%CONSTRUCTOR-method
	
end%classdef

