classdef SimIn
	%UNTITLED Summary of this class goes here
	%   Detailed explanation goes here
	

% Subject: lka
% $Author$
% $LastChangedDate$
% $Revision$
	

	properties
		
		Vehicle
		
		Road
		
		Traffic
		
	end
	
	
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% METHODS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	
	%%% CONSTRUCTOR
	methods
		
		%%% Constructor
		function obj = SimIn(vhcl)
			
			obj.Vehicle = vhcl;
			
			obj.Road = [];
			
			obj.Traffic = [];
		
		end%fcn
		
	end%CONSTRUCTOR-methods
	
end

