function S = paramFile2Struct(parameterFile)
% PARAMFILE2STRUCT  Loads variables from m-file to structure
% 
%   S = PARAMFILE2STRUCT(PARAMETERFILE) creates the structure S with
%   fieldnames and values taken from variables in m-file PARAMETERFILE.
% 
%
%   See also .
% 

% Subject: general purpose
% $Author$
% $LastChangedDate$
% $Revision$



%%% check input arguments
% check if input argument 'parameterFile' is of class char
if ~ischar(parameterFile); error('1st input argument not of class char'); end



%%% load parameter
eval(parameterFile);
S = param2struct();


end%fcn



function S = param2struct(~)
%PARAM2STRUCT   
% 
%   

% get the list of variables from caller workspace
callerspaceVarStrings = evalin('caller','who');

% % remove 'parameterFile' from that list
% relevantVarStrings = setdiff(callerspaceVarStrings,'parameterFile');
relevantVarStrings = callerspaceVarStrings;

% create structure of parameters
for i = 1:length(relevantVarStrings)

    varString_i = relevantVarStrings{i};
    S.(varString_i) = evalin('caller',varString_i);
    
end%for

end%fcn

