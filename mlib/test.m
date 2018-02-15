function tests = test()
%
% This unittest can be executed as
%   runtests('test')
% This matlab unittest framework is available from release 2013a
%
    addpath ../mlib
    addpath ../utilities

    rng(1234);
    
    tests = functiontests(localfunctions);
end    

function test_DR(testCase)

m = 10;
n = 30;
A = randn(m, n);
b = randn(m, 1);

cvx_begin
cvx_quiet(true)
variable x_cvx(n)
minimize norm(x_cvx, 1)
A*x_cvx == b
cvx_end

proxth = @(x) S(x, 1);
f = @(a) norm(a, 1);
PQ = @(x) x + A'*((A*A')\(b - A*x));

[x, fxk, nxk] = dr(proxth, PQ, randn(n, 1), 1e5, 1.0, 1e-32, f);

testCase.assertEqual(x, x_cvx, 'absTol', 1e-7)
testCase.assertEqual(f(x), f(x_cvx), 'absTol', 1e-8)

end

function test_DR_SLP(testCase)
% Solves a small instances of the 
% sparse linear prediction problem
% 
% Used for debugging/verfication 
% of the Douglas-Rachford splitting method
% applied to the sparse linear prediction problem

% Tobias Lindstroem Jensen
% tlj@es.aau.dk, Aalborg University, September 14, 2014



frame = 200; k = 20; n = 30; order = 40;
xin = filter(randn(n, 1), 1, randn(frame*100, 1));

s=[0 xin((1+(k-1)*frame):k*frame)' zeros(1,order-1)];
x=[s(2:end) 0]';
X=toeplitz(s, zeros(1, order));
X(1:order,1:order)=tril(X(1:order,1:order));
gamma = 20.0;

cvx_begin
cvx_quiet(true)
variables a_cvx(order)
minimize norm(x-X*a_cvx,1)+gamma*norm(a_cvx,1)
cvx_end

f = @(a) norm(x-X*a,1)+gamma*norm(a,1);

[a_rd1, opt] = dr_slp(x, size(X), gamma, 1e-9, false, true);

[a_rd2, opt] = dr_slp(x, X, gamma, 1e-9, false, false);


testCase.assertEqual(a_rd1, a_rd2, 'relTol', 1e-10)
testCase.assertEqual(a_rd1, a_cvx, 'absTol', 5e-3)

testCase.assertEqual(f(a_cvx), f(a_rd2), 'relTol', 1e-3)


end



function test_ADMM_SLP(testCase)
% Solves a small instances of the 
% sparse linear prediction problem
% 
% Used for debugging/verfication 
% of the ADMM method
% applied to the sparse linear prediction problem

% Tobias Lindstroem Jensen
% tlj@es.aau.dk, Aalborg University, September 14, 2014


frame = 200; k = 20; n = 30; order = 40;
xin = filter(randn(n, 1), 1, randn(frame*100, 1));

s=[0 xin((1+(k-1)*frame):k*frame)' zeros(1,order-1)];
x=[s(2:end) 0]';
X=toeplitz(s, zeros(1, order));
X(1:order,1:order)=tril(X(1:order,1:order));
gamma = 20.0;

cvx_begin
cvx_quiet(true)
variables a_cvx(order)
minimize norm(x-X*a_cvx,1)+gamma*norm(a_cvx,1)
cvx_end

f = @(a) norm(x-X*a,1)+gamma*norm(a,1);

[a, opt] = admm_slp(x, order, gamma, 1e5, 1e-32);

testCase.assertEqual(a, a_cvx, 'absTol', 1e-5)

testCase.assertEqual(f(a), f(a_cvx), 'relTol', 2e-9)

end

