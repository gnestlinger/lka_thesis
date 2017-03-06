
% number of multiplications
num = 1e5;

A1 = rand(2,1000);
M = rand(2,2);
tic
for i = 1:num
	M*A1;
end
toc


A2 = A1';
tic
for i = 1:num
	A2*M;
end
toc


fprintf('\n')