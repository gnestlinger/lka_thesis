function nsys = ssNormForm(inpArg,type)
% ssNormForm    conversion to state-space normal/canonical forms
%   _______
%   Syntax:
%   nsys = ssNormForm(inpArg,type)
%   ________________
%   Input arguments:
%   inpArg ... transfer function model (tf) or state-space model (ss)
%   type ..... select type of normal form (con/obs)
%   _________________
%   Output arguments:
%   nsys ... state-space model in normal form specified by input argument 
%            'type'
%   _______
%   'con': Controller normal form
%    A = [0     1     0 ... 0     0
%         0     0     1 ... 0     0
%         0                       0
%         0     0     0 ... 0     1
%        -a0   -a1   -a2...   -a_{n-1}]
%   B = [0;0; ... ;0;1];
%   C = [b0 - bn*a0, b1 - bn*a1, ... , b{n-1} - bn*a_{n-1}];
%   D = bn;
% 
%   'obs': Observer normal form
%    A = [0 0 ... 0   -a0
%         1 0 ... 0     .
%         0 1     .     .
%         .       0     .
%         0 ... 0 1 -a_{n-1}]
%   B = [b0 - bn*a0; b1 - bn*a1; ... ; b{n-1} - bn*a_{n-1}];
%   C = [0, ... , 0, 1];
%   D = bn;
%   
% 
% Subject: general purpose
% Author: georgnoname
% Date: 08.02.2013 - 15.03.2013


% check input arguments
if ~ischar(type); error('2nd input argument not of type char'); end%if

% handle input argument
switch class(inpArg)
    case {'tf'}
        trFcn = inpArg;
        
    case {'ss'}
        trFcn = tf(inpArg);
        
    otherwise
        error('Unsupported input argument. Only tf- and ss-models are supported');
end%switch

% get numerator/denumerator
b = trFcn.num{1}; b = b(:);
a = trFcn.den{1}; a = a(:);

% normalize, denumerator polynomial has to be 'monisch'
if a(1) ~= 1
    b = b/a(1);
    a = a/a(1);
end%if

% call subfunction
switch type
    case {'con'} % Controller normal form
        ret = observerNormForm(b,a);
        nsys = ss(ret.a',ret.c',ret.b',ret.d);
        
    case {'obs'} % Observer normal form
        nsys = observerNormForm(b,a);
end%switch

end%fcn


function sys = observerNormForm(b,a)
% Observer normal form

% init state space modell
sys = ss([],[],[],[]);

% degree of system
n = length(a)-1;

% set state space values
sys.A = [[zeros(1,n-1);eye(n-1)], -flipud(a(2:end))];
sys.B = flipud( b(2:end) - b(1)*a(2:end) );
sys.C = [zeros(1,n-1),1];
sys.D = b(1);

end%fcn

