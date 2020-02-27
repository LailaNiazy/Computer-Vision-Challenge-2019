%   POLYFIT Fit polynomial to data.
%
%   [P,S,MU] = POLYFIT(X,Y,N) finds the coefficients of a polynomial in
%   XHAT = (X-MU(1))/MU(2) where MU(1) = MEAN(X) and MU(2) = STD(X). This
%   centering and scaling transformation improves the numerical properties
%   of both the polynomial and the fitting algorithm.

%   Copyright 1984-2005 The MathWorks, Inc.
%   $Revision: 5.17.4.8 $  $Date: 2006/06/20 20:11:56 $

% The regression problem is formulated in matrix format as:
%
%    y = V*p    or
%
%          3  2
%    y = [x  x  x  1] [p3
%                      p2
%                      p1
%                      p0]
%
% where the vector p contains the coefficients to be found.  For a
% 7th order polynomial, matrix V would be:
%
% V = [x.^7 x.^6 x.^5 x.^4 x.^3 x.^2 x ones(size(x))];

function [p,S,mu] = polyfit2(x,y,n)

if ~isequal(size(x),size(y))
    error('MATLAB:polyfit:XYSizeMismatch',...
          'X and Y vectors must be the same size.')
end

x = x(:);
y = y(:);

if nargout > 2
    mu = [mean(x); std(x)];
    x = (x - mu(1))/mu(2);
end

% Construct Vandermonde matrix.
V(:,n+1) = ones(length(x),1,class(x));
for j = n:-1:1
    V(:,j) = x.*V(:,j+1);
end

% Solve least squares problem.
[Q,R] = qr(V,0);
p = R\(Q'*y);     % Same as p = V\y;
r = y - V*p;
p = p.';          % Polynomial coefficients are row vectors by convention.

% S is a structure containing three elements: the triangular factor from a
% QR decomposition of the Vandermonde matrix, the degrees of freedom and
% the norm of the residuals.
S.R = R;
S.df = max(0,length(y) - (n+1));
S.normr = norm(r);