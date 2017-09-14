function sol = MultiGammaLn(nu,K)
% Multivariate gamma function
% FORMAT sol = MultiGammaLn(nu,K)
%
% nu  - Degrees of freedom
% K   - Dimension
% sol - Multivariate gamma
%
% See https://en.wikipedia.org/wiki/Multivariate_gamma_function
%
%__________________________________________________________________________
% Copyright (C) 2017 Wellcome Trust Centre for Neuroimaging

% John Ashburner
% $Id$

sol = K*(K-1)/4*log(pi);
for k=1:K
    sol = sol + gammaln(nu+0.5*(1-k));
end


