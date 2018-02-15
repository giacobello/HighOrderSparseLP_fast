function [x, fxk, nxk] = dr(proxth, PQ, z, kmax, rho, epsilon, g)
% function [x, fxk, nxk] = dr(proxth, PQ, z, kmax, rho, epsilon, g)
% Douglas-Rachford splitting for the problem
%
%  minimize   h(x)
%  subject to x \in Q
%
% Input
%   proxth : the prox operator for h and
%   PQ     : the projection onto Q
%   z      : initialization point
%   kmax   : maximum number of iterations
%   rho    : relaxation parameter (0 < rho < 2)
%   epsilon: accuracy
%   g      : call back function for the iterates
%
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

fxk = zeros(kmax, 1);
nxk = zeros(kmax, 1);

k = 0;
while k < kmax
  x = proxth(z);
  y = PQ(2*x - z);
  zp = z + rho*(y - x);
  
  % save some info
  k = k + 1;
  fxk(k) = g(x);
  nxk(k) = norm(zp-z, 2)^2/length(z);

  % test stopping condition
  if nxk(k) <= epsilon
    fxk = fxk(1:k);
    nxk = nxk(1:k);
    break
  end
  
  % update
  z  = zp;  
end  

