function [s0,s1,s2,mat] = SuffStats(dat,s)
% Compute sufficient statistics
% FORMAT [s0,s1,s2,mat] = SuffStats(dat)
%
% dat - Structure containing various information about each image.
%       Fields for each image n are:
%       dat(n).f - Image data.
%       dat(n).z - Expectations of latent variables.
%       dat(n).S - Covariances of latent variables.
%
% s0  - Zeroeth moment (number of observations)
% s1  - First moment (sum over observations)
% s2  - Second moment (sum of squares of observations)
% mat - Voxel to world mapping
%
% This function is used as part of initialising the shape and appearance
% model fitting.
%__________________________________________________________________________
% Copyright (C) 2017 Wellcome Trust Centre for Neuroimaging

% John Ashburner
% $Id$

s0 = single(0);
s1 = single(0);
s2 = single(0);
mat = eye(4);
for n=1:numel(dat)
    f        = GetDat(dat(n),s);
    d        = [size(f),1,1];
    rho      = ones(d(1:3),'single');
    f1       = single(f);
    msk      = any(~isfinite(f),4);
    f1(repmat(msk,[1 1 1 d(4)])) = 0;
    rho(msk) = 0;
    s0       = s0 + rho;
    s1       = s1 + f1;
    s2       = s2 + f1.*f1;
end


