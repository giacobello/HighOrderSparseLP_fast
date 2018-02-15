function [alpha,opt] = admm_slp(x0, n, gamma, k_max, eps_rel, varargin)
% ADMM_SLP  Alternating Direction Method of Multipliers (ADMM) algorithm
%           for solving l1-regularized linear regression problem
%           on the form
%
%           minimize   ||x-Xa||_1 + gamma ||a||_1
%              a
%
% INPUT   
%   x0     : input speech frame [m x 1]
%   n      : SLP model order
%   gamma  : l1-regularization parameter
%   k_max  : maximum number of iterations
%   eps_rel: relative accuracy
%  settings: (Optional) a struct with the following field
%          rho:  Scaling parameter. (Default rho=100)
%
% OUTPUT
%   alpha  : SLP coefficient vector estimate [n x 1]
%   opt    : A struct with the following fields
%          e: SLP residual vector [m x 1]
%          k: number of iterations performed
%          norm_r: primal residual norm (all iterations) [k x 1]
%          norm_s: dual residual norm (all iterations) [k x 1]
%          fxk: objective at the last iteration
%
%  -------------------------------------------------------------------------
%
%  Copyright 2014-2015 
%  Toon van Waterschoot (toon.vanwaterschoot@esat.kuleuven.be)
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

rho = 100.0; % augmented Lagrangian parameter

if length(varargin) == 1
  settings = varargin{1};
  if isfield(settings, 'rho')
    rho = settings.rho;
  end
end

%% Initialization %%
    %-- ADMM algorithm parameters --
    m = length(x0); % input frame size
    X = toeplitz([0;x0;zeros(n-1,1)],zeros(1,n)); % input speech Toeplitz matrix
    x = [x0;zeros(n,1)]; % input speech frame with zero padding
    
    %-- initialization of internal ADMM variables --    
    alpha = zeros(n,1); % first primal variable (top, scaled) 
                        %= SLP coefficient vector estimate
    e = zeros(m+n,1); % first primal variable (bottom) = SLP residual vector
    y = zeros(m+2*n,1); % second primal variable 
    u = zeros(m+2*n,1); % dual variable (Lagrange multiplier)
    norm_r = zeros(k_max,1); % primal residual norm (all iterations)
    norm_s = zeros(k_max,1); % primal residual norm (all iterations)
    fxk = zeros(k_max,1);  % objective (all iterations)
    norm_r_k = Inf; % primal residual norm (iteration k)
    norm_s_k = Inf; % dual residual norm (iteration k)

    %% Algorithm %%
    %-- calculation of iteration-independent variables --
    %-- solution to regularized l2-norm LP using Levinson-Durbin --
    %-- autocorrelation function estimation using autocorrelation method --
    r = xcorr(x, x, n);
    r_x = r(n+1:end);
    
    %-- modified autocorrelation function calculation --
    r_x_modified = r_x;
    % modified autocorrelation function for regularized l2-norm LP
    r_x_modified(1) = r_x_modified(1) + gamma^2; 
    
    %-- Levinson-Durbin algorithm --
    % regularized l2-norm LP 
    [alpha_gamma2,kappa,sigma2] = levinsondurbin(r_x_modified); 
    alpha_gamma2=-alpha_gamma2(2:end);
    
    %-- ADMM algorithm --
    k = 1; % initialization of iteration index
        
    % check termination criterion
    while ((~((norm_r_k < eps_rel)&&(norm_s_k < eps_rel))) && k <= k_max) 
       
        % solve first primal problem > SLP coefficient vector
        % estimation
        alpha = alpha_gamma2 - Pinv(y-u, X, gamma, r_x_modified);

        % then > SLP prediction error filtering
        e = filter([1;-alpha],1,x); %x-X*alpha;                                
              
        z = [gamma*alpha;e]; % construct full first primal variable

       % save previous estimate of second primal variable 
       %(for calculating dual residual vector)
       y_old = y;            
       % solve second primal problem: shrinkage operation
       y = Sm(z + u, 1/rho); 
       
       alphay = y(1:n)*(1/gamma);

       fxk(k) = norm(filter([1;-alphay],1,x), 1) + gamma*norm(alphay, 1);
       
       % update dual variable
       u = u + z - y;

       % calculate primal residual norm
       norm_r_k = norm(z - y)^2*(1/length(z)); ...
           
       % calculate dual residual norm
       norm_s_k = (rho/length(y))*norm(y - y_old)^2; ...
           
       % save primal residual norm
       norm_r(k) = norm_r_k;      
       
       % save dual residual norm
       norm_s(k) = norm_s_k;                                          

       % increase iteration index
       k = k+1;                 
       
     end
     
     alpha = y(1:n)/gamma;
    
     opt.norm_r = norm_r(1:k-1); % primal residual norm (all iterations)
     opt.norm_s = norm_s(1:k-1); % dual residual norm (all iterations)
     opt.e = filter([1;-alpha], 1, x); % x - X*alpha;
     opt.fxk = fxk(1:k-1);
     opt.k = k-1;


function y = Pinv(z, X, gamma, a);
% Calculates
% y <- pinv([-gamma*eye(n);X])*z;
%
% Where a represent the coeffients in the symmetric Toeplitz 
% matrix T = X'*X+gamma^2I
%

[m, n] = size(X);
zb = -gamma*z(1:n) + X'*z(n+1:end);

s = 1/a(1);
y = levinson(n, s*a', s*zb')';