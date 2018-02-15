function [a opt] = dr_slp(x, X, gamma, epsilon, verbose, flag, varargin)
% function [a opt] = dr_slp(x, X, gamma, epsilon, verbose, flag, varargin)
%
% Solves
%
% minimize   ||x-Xa||_1 + gamma ||a||_1
%    a
%
% using the Douglas-Rachford splitting method.
%
% X is a convolution matrix such that X'*X is Toeplitz
% and symmetric. If flag is true, the algorithm will use
% the Levinson algorithm and argument X is then the dimension [m, n] of X
% and it is implicit assumed that
% X = Toeplitz([0, x], zeros(1,length(alpha))).
%
% Input:
%    x:
%    X:   An m x n matrix if flag is false. Otherwise a two array [m, n]
%          containing the dimension of X
%    gamma: Regularization parameter (positive)
%    epsilon: Accuracy parameter. Stops if this is reached before 
%           kmax iterations
%    verbose: Controls level of output
%    flag: Controls the contend of input X
%    settings: (Optional) a struct with one or more of the following fields
%      t:    stepsize (default 0.1)
%      rho:  over/under relaxation parameter 0 < rho < 2 (default 1.8)
%      kmax: maximum number of iterations (default 2000)
%
%  -------------------------------------------------------------------------
%
%   Copyright 2014-2015 Tobias L. Jensen <tlj@es.aau.dk>
%   Department of Electronic Systems, Aalborg University, Denmark
%
%  This file is part of slp_sm.
%    
%  slp_sm is free software: you can redistribute it and/or modify
%  it under the terms of the GNU General Public License as published by
%  the Free Software Foundation, either version 3 of the License, or
%  (at your option) any later version.
%  
%  slp_sm is distributed in the hope that it will be useful,
%  but WITHOUT ANY WARRANTY; without even the implied warranty of
%  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%  GNU General Public License for more details.
%  
%  You should have received a copy of the GNU General Public License
%  along with slp_sm.  If not, see <http://www.gnu.org/licenses/>.

if flag 
  m = X(1); n = X(2);
else
  [m, n] = size(X);
end

t = .1;
rho = 1.8;
kmax = 2000;

if length(varargin) == 1
  settings = varargin{1};
  if isfield(settings, 't')
    t = settings.t;
  end
  if isfield(settings, 'rho')
    rho = settings.rho;
  end
  if isfield(settings, 'kmax')
    kmax = settings.kmax;
  end  
end

proxth = @(xb) [Sm(xb(1:n), t*gamma); x-Sm(x-xb(n+1:end), t*1)];

if flag
  % Using the symmetric Toeplitz structure
  r = xcorr(x, x, n);
  a = r(n+1:end-1); 
  s = 1/(a(1) + 1);
  a(1) = a(1) + 1;
  a = a*s;
  PQ = @(xb) PQToeplitz(a, s, x, X, xb(1:n), xb(n+1:end), n);
  f = @(a) norm(x - Xa(a, x), 1) + gamma*norm(a, 1);
else
  C = chol(eye(n) + X'*X);
  PQ = @(xb) PQimp(C, X, xb(1:n), xb(n+1:end), n);
  f = @(a) norm(x - X*a, 1) + gamma*norm(a, 1);
end


g = @(x) f(x(1:n));

[xb, fxk, nxk] = dr(proxth, PQ, zeros(m + n, 1), kmax, rho, epsilon, g);
a = xb(1:n);
opt.status = '';
opt.fxk = fxk;
opt.nxk = nxk;


% ----------------------------------------------------------------
% These are utility functions
function save(k, l, f)
l(k) = f;


function x = PQimp(C, X, x1, x2, n)
v = C\(C'\(x1 + X'*x2));
x = [v; X*v];


function x = PQToeplitz(a, s, x, X, x1, x2, n)
% The used version of the Levinson algorithm 
% needs the diagonal to be normalized, a contains the scaled
% coefficent, and s is the scaling that is needed for the 
% right hand side T x = s b

%X'*x2
XTx2p = filter(flipud(x2), 1, x);
XTx2 = XTx2p(end-1:-1:end-length(a));
v = levinson(n, a', s*(x1 + XTx2)')';
Xvp = filter(v, 1, x);
x = [v; [0; Xvp(1:end-1)]];


function b = Xa(a, x)
Xap = filter(a, 1, x);
b = [0; Xap(1:end-1)];