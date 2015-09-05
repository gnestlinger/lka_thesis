function out = loadParameter(parameterFile,varargin)
% loadParameter     load parameter from file to structure
%   _______
%   Syntax:
%   out = loadParameter(parameterFile,varargin)
%   ________________
%   Input arguments:
%   parameterFile ... m-file string
%   varargin ........ string(s) of variable used in 'parameterFile' (opt.)
%   _________________
%   Output arguments:
%   out ... structure containing variables from 'parameterFile' stored in
%           out.values.x1
%           out.values.x2
%           out.values. ...
%           
%   out = loadParameter(parameterFile) stores all variables that the m-file
%   'parameterFile' contains in the structure 'ret.values'.
% 
%   out = loadParameter(parameterFile,varargin) does the same as
%   loadParameter(parameterFile) except that the variables with names given
%   in 'varargin' are stored in a substruct like
%       out.varargin{1}
%       out.varargin{2}
%       out. ...
% 
%   All input arguments have to be strings.
% 

% Subject: general purpose
% Author: georgnoname
% Date: 08.10.2012 - 17.02.2015


% check input arguments
if ~ischar(parameterFile); error('1st input argument not of type char'); end
for iLoop = 1:size(varargin,2)
    if ~ischar(varargin{iLoop}); 
        error(['Input argument ',varargin{iLoop},' not of type char']); 
    end
end%for
clear iLoop;


% load parameter
eval(parameterFile);


% get list of variables in workspace
initVarList = who;


% check if variable 'iLoop' (will be used as for-loop counter variable)
% is already used as variable name in 'parameterFile'
if any(strcmp('iLoop',initVarList))
% if ~isempty(strmatch('iLoop',initVarList,'exact')) % pre R2012a
    error(['This function is using a loop counter variable ''iLoop'' ',...
        'which is also used in ''',parameterFile,'.m''.',...
        10,...
        'Change the variable name in ''',parameterFile,'.m'' or ',...
        'the counter loop variable name to make the function working!'])
end%if


% strings of input argument filenames
inpArgsStrings = {'parameterFile','varargin'};


% get indices of position of strings 'inpArgsStrings' in 'initVarList'
removeIndex(2,1) = 0;
for iLoop = 1:length(inpArgsStrings)
    removeIndex(iLoop,1) = ...
        find(strcmp(inpArgsStrings{iLoop},initVarList));
%         strmatch(inpArgsStrings{iLoop},initVarList,'exact');
end%for


% if varargin is present
if ~isempty(varargin)
    
    % get indices of position of strings from 'varargin' in 'initVarList'
    infoIndex(length(varargin),1) = 0;
    try % 
        for iLoop = 1:length(varargin)
            infoIndex(iLoop,1) = find(strcmp(varargin{iLoop},initVarList));
%             infoIndex(iLoop,1) = strmatch(varargin{iLoop},initVarList,'exact'); % pre R2012 a
        end%for
    catch exception % error if varargin{iLoop} is not a valid variable name
%         disp(exception.message);
        error(['Variable ''',varargin{iLoop},''' ',...
            'undefined in parameter-file ',parameterFile,'.m.']);
    end%try

    % create structure of info
    for iLoop = 1:length(infoIndex)
        out.(initVarList{infoIndex(iLoop)}) = ...
            eval(initVarList{infoIndex(iLoop)});
    end%for

    % union of 'removeIndex' and 'infoIndex' (contains indices to be
    % removed from 'initVarList')
    removeIndex = [removeIndex;infoIndex];
    
end%if


% remove strings in 'inpArgsStrings' and 'varargin' from 'initVarList'
keepIndex = true(size(initVarList));
keepIndex(removeIndex) = 0;
postVarList = initVarList(keepIndex);


% create structure of values
for iLoop = 1:length(postVarList)
    out.values.(postVarList{iLoop}) = eval(postVarList{iLoop});
end%for


end%fcn
