function ret = col(in)
% COL column vector
% check if 'in' is of dimension 1*x or x*1 (vector, not matrix)

[a,b] = size(in);

if (a>1) && (b>1)
    error('Input Argument should be vector but is matrix.')
end

ret = in(:);

end%fcn