% Solves a small instances of the 
% sparse linear prediction problem
% 
% Mainly used for illustration
% of ADMM and Douglas-Rachford splitting methods

% Tobias Lindstroem Jensen
% tlj@es.aau.dk, Aalborg University, August 21, 2014


addpath ../mlib
addpath ../utilities

close all
clear all

frame = 200; k = 20; n = 10; order = 40;
xin = filter(randn(n, 1), 1, randn(frame*100, 1));

s = [0 xin((1+(k-1)*frame):k*frame)' zeros(1,order-1)];
x = [s(2:end) 0]';
X = toeplitz(s, zeros(1, order));
gamma = 20.0;

cvx_begin
cvx_quiet(true)
variables a_cvx(order)
minimize norm(x-X*a_cvx,1) + gamma*norm(a_cvx,1)
cvx_end
fprintf('\n');

f = @(a) norm(x-X*a,1)+gamma*norm(a,1);

settings.kmax = 100;
[a_dr, opt_dr] = dr_slp(x, [frame+order, order], gamma, 1e-6, true, true, settings);

[a_admm, opt_admm] = admm_slp(x(1:frame), order, gamma, settings.kmax, 1e-6);


figure()
semilogy((opt_dr.fxk-f(a_cvx))/abs(f(a_cvx)), 'r')
hold on
semilogy((opt_admm.fxk-f(a_cvx))/abs(f(a_cvx)), 'k')
xlabel('Iteration k')
ylabel('Relative gap')
legend('DR', 'ADMM')

figure()
semilogy(opt_dr.nxk, 'r')
hold on
semilogy(opt_admm.norm_r, 'g')
semilogy(opt_admm.norm_s, 'k')
xlabel('Iteration k')
ylabel('Residuals')
legend('DR (1/N)||z-z+||_2^2','ADMM (1/N)||z-y||_2^2', 'ADMM (1/N)||y-y^-||_2^2')

