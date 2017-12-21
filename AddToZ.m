function dat = AddToZ(dat,za)
% Addition of vector to latent variables
% FORMAT dat = AddToZ(dat,za)
%
% dat - Structure containing various information about each image.
%       Fields for each image n are:
%       dat(n).f - Image data.
%       dat(n).z - Expectations of latent variables.
%       dat(n).S - Covariances of latent variables.
% za  - vector to add
%
%__________________________________________________________________________
% Copyright (C) 2017 Wellcome Trust Centre for Neuroimaging

% John Ashburner
% $Id$

N = numel(dat);
for n=1:N
    dat(n).z = dat(n).z + za;
end

