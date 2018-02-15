function tests = test()
%
% This unittest can be executed as
%   runtests('test')
% This matlab unittest framework is available from release 2013a
%

    rng(1234)
    tests = functiontests(localfunctions);

end    

function test_S_simple(testCase)
% Testing of softhresholding

n = 100;
x = randn(n, 1);
t = 1.1;

%% Case for prox_t||u||_1 (x)
cvx_begin
cvx_quiet(true)
variable u(n)
minimize t*norm(u, 1) + 0.5*sum_square(u - x)
cvx_end

uh = S(x, t);
testCase.assertEqual(u, uh, 'absTol', 1e-4)

uh2 = Sm(x, t);
testCase.assertEqual(uh, uh2, 'absTol', 1e-8)


end


function test_S_shift(testCase)
n = 100;
x = randn(n, 1);
t = 1.1;

%% Case for prox_t||v - u||_1 (x) = v-prox_t||u||(v-x)
x = randn(n, 1);
v = randn(n, 1);

cvx_begin
cvx_quiet(true)
cvx_precision high
variable u(n)
minimize t*norm(x - u, 1) + 0.5*sum_square(u - v)
cvx_end

uh = x - S(x - v, t);
testCase.assertEqual(u, uh, 'absTol', 1e-5)

end

function test_levinson(testCase)
%% test for solving using symmetric Toeplitz with general right-hand side

frame = 200; k = 20; n = 12; order = 40;
xin = filter(randn(n, 1), 1, randn(frame*100, 1));

s=[0 xin((1+(k-1)*frame):k*frame)' zeros(1,order-1)];
x=[s(2:end) 0]';
X=toeplitz(s, zeros(1, order));
X(1:order,1:order)=tril(X(1:order,1:order));

R = X'*X;


r = xcorr(x, x, order);
r = r(order+1:end-1);
s = (1/r(1));
r = r*s;

testCase.assertEqual(R(:,1)*s, r, 'relTol', 1e-8)

b = randn(order, 1);
sol = R\b;

solh1 = levinson(order, r', s*b')';
testCase.assertEqual(sol, solh1, 'relTol', 1e-8)
end


function test_levinson_durbin(testCase)
%% test for solving the Yule-Walker equations using Levinson-Durbin

n = 10;

r = randn(n, 1);
R = toeplitz(r(1:end-1));

a = -R\r(2:end);

ah = levinsondurbin(r);

testCase.assertEqual(a, ah(2:end), 'relTol', 1e-12)
testCase.assertEqual(a, ah(2:end), 'absTol', 1e-12)

end