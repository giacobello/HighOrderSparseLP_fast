function [a,varargout] = levinsondurbin(correlation_vector,varargin) 

% LEVINSONDURBIN  Levinson-Durbin recursion for AR identification
%	
% A = LEVINSONDURBIN(CORRELATION_VECTOR) performs the Levinson-Durbin
% recursive algorithm for linear prediction on the correlation function
% estimates specified in CORRELATION_VECTOR. The parameter vector A
% containing the linear prediction coefficients is returned (A has the same
% size as CORRELATION_VECTOR).
%
% [A,SIGMA2] = LEVINSONDURBIN(CORRELATION_VECTOR) also returns the residual
% signal variance SIGMA2.
%
% [A,KAPPA,SIGMA2] = LEVINSONDURBIN(CORRELATION_VECTOR) also returns the
% reflection coefficients vector KAPPA.
%
% [A,...] = LEVINSONDURBIN(CORRELATION_VECTOR,ORDER) allows for performing
% a linear prediction of lower order than length(CORRELATION_VECTOR)-1.
%
% Author: Toon van Waterschoot (toon.vanwaterschoot@esat.kuleuven.be)
% Reference: S. Haykin, "Adaptive filter theory", 3rd edition, 
%   Prentice-Hall, pp. 254-264.
% Date: 22/11/06
% Version: 1-0

if nargin == 2
    order = varargin{1};
else
    order = length(correlation_vector) - 1;
end

kappa = zeros(order,1); % initialization of reflection coefficients vector
sigma2 = correlation_vector(1); % initialization of zero-order prediction error variance
a = [1;zeros(order,1)]; % initialization of zero-order predictor coefficient

for m = 1:order,
    Delta = correlation_vector(m+1:-1:2)'*a(1:m); % m-th order auxiliary variable Delta
    kappa(m) = -Delta/sigma2; % m-th reflection coefficient of m-th order predictor
    sigma2 = sigma2*(1-kappa(m)^2); % m-th order prediction error variance
    a(1:m+1) = [a(1:m);0] + kappa(m).*[0;flipud(a(1:m))]; % m-th order predictor coefficients
end
    
if nargout==2
    varargout{1} = sigma2;
elseif nargout==3
    varargout{1} = kappa;
    varargout{2} = sigma2;
end
