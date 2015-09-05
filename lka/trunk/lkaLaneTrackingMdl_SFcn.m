function lkaLaneTracking_SFcn(block)
%MSFUNTMPL A Template for a MATLAB S-Function
%   The MATLAB S-function is written as a MATLAB function with the
%   same name as the S-function. Replace 'msfuntmpl' with the name
%   of your S-function.  
%
%   It should be noted that the MATLAB S-function is very similar
%   to Level-2 C-Mex S-functions. You should be able to get more 
%   information for each of the block methods by referring to the
%   documentation for C-Mex S-functions.
%  
%   Copyright 2003-2010 The MathWorks, Inc.
%   $Revision: 1.1.6.22 $  
  
%
% The setup method is used to setup the basic attributes of the
% S-function such as ports, parameters, etc. Do not add any other
% calls to the main body of the function.  
%   
setup(block);

end%fcn
  
%endfunction

% Function: setup ===================================================
% Abstract:
%   Set up the S-function block's basic characteristics such as:
%   - Input ports
%   - Output ports
%   - Dialog parameters
%   - Options
% 
%   Required         : Yes
%   C-Mex counterpart: mdlInitializeSizes
%
function setup(block)

	% Register number of input and output ports
    block.NumInputPorts  = 8;
	block.NumOutputPorts = 2;
    
	% Setup functional port properties to dynamically inherited.
	block.SetPreCompInpPortInfoToDynamic;
	block.SetPreCompOutPortInfoToDynamic;
    
    
	% Override the input port properties.
	block.InputPort(1).DatatypeID  = 0; % double
	block.InputPort(1).Complexity  = 'Real';
	block.InputPort(1).Dimensions = 1;
	block.InputPort(1).DirectFeedthrough = false;
	
	block.InputPort(2).DatatypeID  = 0; % double
    block.InputPort(2).Complexity  = 'Real';
	block.InputPort(2).Dimensions = 1;
    block.InputPort(2).DirectFeedthrough = false;
    
	block.InputPort(3).DatatypeID  = 0; % double
    block.InputPort(3).Complexity  = 'Real';
	block.InputPort(3).Dimensions = 1;
    block.InputPort(3).DirectFeedthrough = false;
    
	block.InputPort(4).DatatypeID  = 0; % double
    block.InputPort(4).Complexity  = 'Real';
	block.InputPort(4).Dimensions = 1;
    block.InputPort(4).DirectFeedthrough = false;
    
	block.InputPort(5).DatatypeID  = 0; % double
    block.InputPort(5).Complexity  = 'Real';
	block.InputPort(5).Dimensions = 1;
    block.InputPort(5).DirectFeedthrough = false;
    
	block.InputPort(6).DatatypeID  = 0; % double
    block.InputPort(6).Complexity  = 'Real';
	block.InputPort(6).Dimensions = 1;
    block.InputPort(6).DirectFeedthrough = false;
    
	block.InputPort(7).DatatypeID  = 0; % double
    block.InputPort(7).Complexity  = 'Real';
	block.InputPort(7).Dimensions = 1;
    block.InputPort(7).DirectFeedthrough = false;
  
	block.InputPort(8).DatatypeID  = 0; % double
    block.InputPort(8).Complexity  = 'Real';
	block.InputPort(8).Dimensions = 1;
    block.InputPort(8).DirectFeedthrough = false;
    
    
	% Override the output port properties.
	block.OutputPort(1).DatatypeID  = 0; % double
	block.OutputPort(1).Complexity  = 'Real';
	block.OutputPort(1).Dimensions = 1;
    
	% Override the output port properties.
	block.OutputPort(2).DatatypeID  = 0; % double
	block.OutputPort(2).Complexity  = 'Real';
	block.OutputPort(2).Dimensions = 1;
	
	% Set up the continuous states.
    block.NumContStates = 2;
	
	
	block.NumDialogPrms = 0;  
	
	% Set block sample time to inherited
	block.SampleTimes = [-1 0];
	
	% Set the block simStateCompliance to default (i.e., same as a built-in block)
	block.SimStateCompliance = 'DefaultSimState';

    % Run accelerator on TLC
    block.SetAccelRunOnTLC(true);

    % Register methods
    block.RegBlockMethod('Derivatives',@Derivatives);
    block.RegBlockMethod('InitializeConditions',@InitializeConditions);  
    block.RegBlockMethod('Outputs',@Output);  
    block.RegBlockMethod('SetInputPortSamplingMode',@SetInputPortSamplingMode);
	
end%function

% -------------------------------------------------------------------
% The local functions below are provided to illustrate how you may implement
% the various block methods listed above.
% -------------------------------------------------------------------

function Derivatives(block)

% inputs
vx = block.InputPort(1).Data;
vy = block.InputPort(2).Data;
psiDot = block.InputPort(3).Data;
kapL = block.InputPort(4).Data;
L = block.InputPort(5).Data;
yLExact = block.InputPort(8).Data;

% current state
yL = block.ContStates.Data(1);
epsL = block.ContStates.Data(2);

% derivatives
yLDot = vx*tan(epsL) - vy - L*psiDot;
epsLDot = (vx - yLExact*psiDot)*kapL/cos(epsL) - psiDot;

% assign to derivative property
% block.Derivatives(1).Data = block.ContStates.Data(2);
block.Derivatives.Data = [yLDot; epsLDot];

end%function

function InitializeConditions(block)

block.ContStates.Data = [block.InputPort(6).Data;block.InputPort(7).Data];

end%fcn

function Output(block)

%     block.OutputPort(1).Data = [1;2];
    block.OutputPort(1).Data = block.ContStates.Data(1);
    block.OutputPort(2).Data = block.ContStates.Data(2);
  
end%function

function SetInputPortSamplingMode(block, idx, fd)

block.InputPort(idx).SamplingMode = fd;
block.OutputPort(1).SamplingMode = fd;
block.OutputPort(2).SamplingMode = fd;

end%fcn
