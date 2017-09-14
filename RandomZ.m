function dat = RandomZ(dat,K)
% Initialise latent variables randomly
% FORMAT dat = RandomZ(dat,K)
%
% dat - Data structure
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

