function dat = AddRandZ(dat,Sig)
% Add random nouse to latent variables
% FORMAT dat = AddRandZ(dat,Sig)
%
% dat - Structure containing various information about each image.
%       Fields for each image n are:
%       dat(n).f - Image data.
%       dat(n).z - Expectations of latent variables.
%       dat(n).S - Covariances of latent variables.
% Sig - Noise covariance
%
%__________________________________________________________________________
% Copyright (C) 2017 Wellcome Trust Centre for Neuroimaging

% John Ashburner
% $Id$

randn('seed',1);
C=chol(Sig)';
for n=1:numel(dat)
    dat(n).z  = dat(n).z + single(C*randn(size(C,1),1));
end

