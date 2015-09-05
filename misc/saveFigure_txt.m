function [] = saveFigure_txt(filename,rowSoll,varargin)
% saveFigure_txt    saves data to ascii file
%   _______
%   Syntax:
%   saveFigure_txt(filename,rowSoll,varargin)
%   ________________
%   Input arguments:
%   filename ... string of desired filename (will be extended by '.dat')
%   rowSoll .... max. number of rows to save to file
%   varargin ... list with even number of arguments in the style
%                   'columnName1',data1,'columnName2',data2,...
%   _________________
%   Output arguments:
%   none
%
% Subject: general purpose
% Author: georgnoname
% Date: 14.12.2012 - 10.05.2013


% check input arguments
% number of varargin input arguments (must be even)
nVarArgIn = length(varargin);
if mod(nVarArgIn,2) ~= 0; error('not even'); end

% filename
filename = [filename,'.dat'];

% get all strings (every second element of varargin starting with the
% first)
strings = varargin(1:2:nVarArgIn-1);

% get all data (every second element of varargin starting with the
% second)
data = varargin(2:2:nVarArgIn);

% get length of data-vectors (number of rows)
rowIst(nVarArgIn/2,1) = 0;
for i = 1:nVarArgIn/2
    rowIst(i) = length(data{i});
end%for

% set maximum number of rows 'rowSoll' if unspecified (rowSoll=[])
if isempty(rowSoll)
    rowSoll = max(rowIst);
end%if

% berechne jedes wievielte Element der Vektoren in 'data{i}' in file
% gespeichert wird
% Erstes ('data{i}(1)') und letztes ('data{i}(end)') Element müssen
% erhalten bleiben. Erstes bleibt sowieso erhalten. -> es müssen
% 'rowSoll-1' Schritte der Weite 'delta' in '1,..,rowIst-1' untergebracht
% werden
div = (rowIst-1)/(rowSoll-1);
delta = ceil(div);

% gibt an ob 'data{i}(end)' explizit mitübernommen werden muss oder ob
% durch 'delta' bereits eingeschlossen
takeEnd = ~(delta==div);

% limit data-vectors to length <='rowSoll'
rowIstNeu(nVarArgIn/2,1) = 0;
for i = 1:(nVarArgIn/2)
    if rowIst(i) > rowSoll
        if takeEnd(i)
            data{i} = data{i}([1:delta(i):rowIst(i),end]);
        else
            data{i} = data{i}([1:delta(i):rowIst(i)]);
        end%if  
%     else     
%         dataNeu{i} = data{i};  
    end%if
    
    rowIstNeu(i) = length(data{i});
        
end%for

% pre-allocation
M(max(rowIstNeu),nVarArgIn/2) = 0;
M = M + NaN;

% store data-vectors in matrix
for i = 1:(nVarArgIn/2)
    M(1:length(data{i}),i) = col(data{i});
end%for

% create and open file
fid = fopen(filename,'w');

% write strings to file
fprintf(fid,'%13s\t',strings{:});

% write data to file
dlmwrite(filename,...
    M,...
    'delimiter', '\t',...
    'precision','%+13.5e',...
    'roffset',1,...
    '-append');

% close file
fclose(fid);

% info output
disp(' ');
disp(['Dateiname: ',filename]);
disp(['gewünschte max. Zeilenanzahl: ',num2str(rowSoll)]); 
disp(['erreichte Zeilenanzahl: ',num2str(max(rowIstNeu))])
for i = 1:nVarArgIn/2
    disp([num2str(i,'%2.2i'),'. Vektor: Reduktion von ',num2str(rowIst(i)),' auf ',num2str(rowIstNeu(i))]); 
end%for
disp(' ');

end


function ret = col(in)
% COL column vector
% check if 'in' is of dimension 1*x or x*1 (vector, not matrix)

[a,b] = size(in);

if (a>1) && (b>1)
    error('Input Argument should be vector but is matrix.')
end

ret = in(:);

end%fcn