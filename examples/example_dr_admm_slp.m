% Solves a small instances of the 
% sparse linear prediction problem
% 
% Generates an example figure used in the paper
%
% Tobias Lindstroem Jensen
% tlj@es.aau.dk, Aalborg University, August 21, 2014
% Modified Dec 23, 2014 to include tikz output


addpath ../mlib
addpath ../utilities

clear all

[y fs nbits] = wavread('timit16k.wav');
y = y(7000:8200);
frame=320; %%frame length (multiple of 4)
nframes=1;
y=y(1:nframes*frame); %resizing input speech

%% removal of low frequencies components

fc=50; %cut-off frequency (Hz)
ord=2; %order
rp=20; %ripple amplitude
[b,a] = cheby2(ord,rp,fc/fs,'high');
xin=filter(b,a,y);

order = 250;
settings.kmax = 100;
gamma = .12;

s=[0 xin' zeros(1,order-1)];
x=[s(2:end) 0]';
X=toeplitz(s, zeros(1, order));

cvx_begin
cvx_quiet(true)
variables a_cvx(order)
minimize norm(x-X*a_cvx, 1)+gamma*norm(a_cvx, 1)
cvx_end
fprintf('\n');

f = @(a) norm(x-X*a, 1) + gamma*norm(a, 1);

settings.kmax = 150;
[a_dr, opt_dr] = dr_slp(x, [frame+order, order], gamma, 1e-6, true, true, settings);

[a_admm, opt_admm] = admm_slp(x(1:frame), order, gamma, settings.kmax, 1e-6);


figure(1)
clf
semilogy((opt_dr.fxk-f(a_cvx))/abs(f(a_cvx)), 'r')
hold on
semilogy((opt_admm.fxk-f(a_cvx))/abs(f(a_cvx)), 'k')
xlabel('\xlabel')
ylabel('\ylabel')
legend('\DR', '\ADMM')

try
    matlab2tikz('example_dr_admm.tikz', 'height', '\figureheight', 'width', ...
    '\figurewidth');
end

figure(2)
clf
semilogy(opt_dr.nxk, 'r')
hold on
semilogy(opt_admm.norm_r, 'g')
semilogy(opt_admm.norm_s, 'k')
xlabel('Iteration k')
ylabel('Residuals')
legend('DR ||z-z+||_2^2', 'ADMM ||z-y||_2', 'ADMM ||y-y^-||_2')

