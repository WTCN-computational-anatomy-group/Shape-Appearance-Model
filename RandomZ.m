function dat = RandomZ(dat,K)
% Initialise latent variables randomly
% FORMAT dat = RandomZ(dat,K)
%
% dat - Structure containing various information about each image.
%       Fields for each image n are:
%       dat(n).f - Image data.
%       dat(n).z - Expectations of latent variables.
%       dat(n).S - Covariances of latent variables.
% K   - Dimensionality of latent variables
%
%__________________________________________________________________________
% Copyright (C) 2017 Wellcome Trust Centre for Neuroimaging

% John Ashburner
% $Id$

randn('seed',1);
for n=1:numel(dat)
    dat(n).z  = single(randn(K,1));
    dat(n).S  = zeros(K,K,'single');
end

