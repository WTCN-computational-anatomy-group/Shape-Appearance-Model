function dat = TransfZ(dat,T)
% Matrix pre-multiplication of latent variables
% FORMAT dat = TransfZ(dat,T)
%
% dat - Structure containing various information about each image.
%       Fields for each image n are:
%       dat(n).f - Image data.
%       dat(n).z - Expectations of latent variables.
%       dat(n).S - Covariances of latent variables.
% T   - Matrix to multiply z by
%
% Does the following for each image:
%       dat(n).z = T*dat(n).z;
%       dat(n).S = T*dat(n).S*T';
%__________________________________________________________________________
% Copyright (C) 2017 Wellcome Trust Centre for Neuroimaging

% John Ashburner
% $Id$

N = numel(dat);
for n=1:N
    dat(n).z = T*dat(n).z;
    if isfield(dat,'S')
        dat(n).S = T*dat(n).S*T';
    end
end

