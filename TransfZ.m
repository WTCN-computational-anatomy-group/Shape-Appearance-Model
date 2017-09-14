function dat = TransfZ(dat,Rz)
% Matrix multiplication of latent variables
% FORMAT dat = TransfZ(dat,Rz)
%
% dat - data (with fields dat(n).z & dat(n).S)
% Rz  - Transpose of matrix to multiply z by
%
%__________________________________________________________________________
% Copyright (C) 2017 Wellcome Trust Centre for Neuroimaging

% John Ashburner
% $Id$

N = numel(dat);
for n=1:N
    dat(n).z = Rz'*dat(n).z;
    if isfield(dat,'S')
        dat(n).S = Rz'*dat(n).S*Rz;
    end
end

