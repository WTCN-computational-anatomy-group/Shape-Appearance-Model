function ld = LogDet(A)
% Log-determinant of matrix
% FORMAT ld = LogDet(A)
%
% A  - A square matrix
% ld - Logarithm of determinant of A
%
% Cholesky factorisation is used to compute a more stable log-determinant.
%
%__________________________________________________________________________
% Copyright (C) 2017 Wellcome Trust Centre for Neuroimaging

% John Ashburner
% $Id$

[C,p] = chol(A);
if p>0
   %warning('Attempting to compute log determinant of matrix that is not positive definite (p=%d).', p);
end

ld = 2*sum(log(diag(C)));

