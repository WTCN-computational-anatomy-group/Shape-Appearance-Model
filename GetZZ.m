function varargout = GetZZ(dat)
% Return suffient statistics from latent variables
% FORMAT [N,Z,ZZ,S] = GetZZ(dat)
%
% dat - Structure containing various information about each image.
%       Fields for each image n are:
%       dat(n).f - Image data.
%       dat(n).z - Expectations of latent variables.
%       dat(n).S - Covariances of latent variables.
%
% N   - Number of images.
% Z   - Sum over expectations of latent variables (dat(n).z).
% ZZ  - Sum over dat(n).z*dat(n).z'
% S   - Sum over covariances (dat(n).S)
%
%__________________________________________________________________________
% Copyright (C) 2017 Wellcome Trust Centre for Neuroimaging

% John Ashburner
% $Id$

Z  = single(0);
ZZ = single(0);
S  = single(0);
N  = numel(dat);
for n=1:N
    z  = dat(n).z;
    Z  = Z + z;
    ZZ = ZZ + z*z';
    S  = S + dat(n).S;
end
varargout{1} = N;
varargout{2} = Z;
varargout{3} = ZZ;
varargout{4} = S;

