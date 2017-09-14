function ld = Elogdet(W,nu)
% Expectation of log-determinant of matrix from a Wishart distribution
% FORMAT ld = Elogdet(W,nu)
% W  - Scale matrix
% nu - Degrees of freedom
% ld - Expectation of log-determinant of matrices drawn from this Wishart
%      distribution.
%
% See: 10.65 on page 478 of Bishop's PRML
%      https://en.wikipedia.org/wiki/Wishart_distribution
%__________________________________________________________________________
% Copyright (C) 2017 Wellcome Trust Centre for Neuroimaging

% John Ashburner
% $Id$

K   = size(W,1);
ld = 2*sum(log(diag(chol(W)))); % log(det(W))
ld = ld + K*log(2);
for k=1:K
    ld = ld + psi(0.5*(nu+1-k));
end

